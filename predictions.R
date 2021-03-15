## Load libraries
library(tidyverse)
## For models
library(MCMCglmm)
library(postMCMCglmm)
## For graphing
library(cowplot)
library(GGally)

## Read in data
data <- read_csv("Data/glmm_dataset.csv") %>% rename(PA = pres_abs)
metadata <- read_csv("Data/metadata_imputed.csv")
load("Models/Model3_noPhylo.Rdata")
model_noPhylo <- mod; rm(mod)

######## Plot correlations between variables ######## 

## Select variables in use
mdat <- select(metadata, elevation, percentcover, average.soil.temperature,
               average.soil.relative.humidity, ec)
names(mdat) <- c("Elevation", "Vegetation", "Soil Temperature", "Soil Humidity", "Conductance")

## Plot pairs graph
pairsPlot <-  ggpairs(mdat, lower = list(continuous = wrap("smooth",alpha = 0.3, size=0.1))) + theme_bw()

alphabeta_plot <- function(covar, mod, data){
  ## MCMC chain for intercept
  chain_i <- mod$Sol[,1]
  ## MCMC chain for Covariate - fixed effect
  chain_cov <- mod$Sol[,colnames(mod$Sol) == covar]
  ## MCMC chain for covariate - random slope variance
  var_slope_col <- which(str_count(colnames(mod$VCV), covar) == 2)
  chain_var_slope <- mean(mod$VCV[,var_slope_col])
  ## Get data column number for mean-centred covariate
  column_sc <- which(colnames(data) == covar)
  ## And also the original covariate
  unscaled_var <- str_split_fixed(covar, "_", n = 2)[[1]]
  column_unsc <- which(colnames(data) == unscaled_var)
  
  ## Generate predictions using means and quantiles
  dat <- expand.grid(param = mean(chain_cov), # Slope estimate for covariate fixed effeect
                     lparam = quantile(chain_cov, probs = 0.025), # Slope upper CI
                     uparam = quantile(chain_cov, probs = 0.975), # Slope lower CI
                     urslope = qnorm(0.975, mean = 0, sd = sqrt(chain_var_slope)), # random slope quantiles
                     lrslope = qnorm(0.025, mean = 0, sd = sqrt(chain_var_slope)),
                     ## Generate 100 covariate predictions on mc scale
                     covar = seq(min(data[,column_sc]), max(data[,column_sc]), length.out = 100)) %>% 
    ## Calculate predictions
    mutate(pred = plogis(mean(chain_i) + param*covar), # Fixed slope 
           lpred = plogis(mean(chain_i) + lparam*covar), # Fixed slope lci
           upred = plogis(mean(chain_i) + uparam*covar), # Fixed slope uci
           urslope_pred = plogis(mean(chain_i) + (urslope + param)*covar), # Random effect uci
           lrslope_pred = plogis(mean(chain_i) + (lrslope + param)*covar), # Random effect lci
           ## Convert mean-centred covar to original units
           covar_unscaled = (covar * sd(as.data.frame(data)[,column_unsc])) + 
             mean(as.data.frame(data)[,column_unsc]))
  
  ## Create plot for fixed effect
  alpha_plot <- ggplot(dat, aes(x = covar_unscaled, y = pred)) +
    geom_ribbon(aes(ymin = lpred, ymax = upred), fill = "darkgoldenrod1", alpha = 1) +
    geom_line() +
    labs(y = "Probability of ASV occurrence", x = names(covar)) +
    lims(y = c(0, 0.015))+
    theme_bw()
  
  ## Create plot for turnover effect
  beta_plot <- ggplot(dat, aes(x = covar_unscaled)) +
    geom_ribbon(aes(ymin = lrslope_pred, ymax = urslope_pred), colour = "darkgreen", alpha = 0) +
    labs(y = "Probability of ASV occurrence", x = names(covar)) +
    lims(y = c(0, 0.27))+
    theme_bw()
  
  return(list(alpha_plot, beta_plot))
}

## Generate graphs for each var of interest:
variables <- c("Elevation" = "elevation_mc", 
               "Vegetation Cover (%)" = "percentcover_mc", 
               "Soil Temperature" = "average.soil.temperature_mc", 
               "Soil Humidity" = "average.soil.relative.humidity_mc",
               "Conductance" = "ec_mc")

Plots <- list()
for(i in 1:length(variables)){
  Plots[[i]] <- alphabeta_plot(variables[i], model_noPhylo, data = data)
}
Plots <- unlist(Plots, recursive = FALSE)

plot_grid(plotlist = Plots, nrow = 2, byrow = FALSE)

########## Taxonomy ##########

## Read taxonomy data
taxonomy <- qiime2R::read_qza("Data/taxonomy.qza")
taxonomy <- taxonomy$data %>%
  separate(Taxon, sep = "; ",
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  rename(ASV = Feature.ID) #%>% 
  # mutate(Phylum = sub("p__", "", Phylum))
taxonomy[is.na(taxonomy$Phylum),]$Phylum <- "Unknown"
taxonomy[is.na(taxonomy$Family),]$Family <- "Unknown"

## Get ASV posterior means for elevation
ASVmeans <- as.data.frame(apply(model_noPhylo$Sol, 2, mean)) %>% 
  rownames_to_column("ASV") %>% 
  filter(grepl("elevation", ASV)) %>% 
  filter(grepl("ASV", ASV)) %>% 
  mutate(ASV = sub("elevation_mc.ASV.", "", ASV))
names(ASVmeans) <- c("ASV", "ElevSlope")

taxonomy <- left_join(taxonomy, ASVmeans) %>% group_by(Phylum) %>% 
  mutate(mean = median(ElevSlope)) %>% 
  filter(!(Phylum %in% c("RCP2-54", "Sumerlaeota", "Fibrobacterota", "WPS-2", "Desulfobacterota",
                         "Deinococcota", "MBNT-15", "MBNT15")))

taxonomy_plot <- ggplot(taxonomy, aes(y = forcats::fct_reorder(Phylum, -mean), x = ElevSlope, fill = Phylum)) +
  geom_boxplot() +
  geom_vline(xintercept = 0, lty = 2) +
  lims(x = c(-0.5, 0.5)) +
  labs(x = "Species response to elevation", y = "Phylum") +
  theme_bw() +
  theme(legend.position = "none")

## Model checking

# sims <- simulate(model_noPhylo, 1000)
load("Models/sims.Rdata")
sims_nzero <- apply(sims, 2, function(x) { sum(x ==0) })

sims_nzero_plot <- qplot(sims_nzero) + 
  geom_vline(xintercept = sum(data$PA == 0), colour = "red") +
  labs(x = "Number of zeroes", y = "Number of simuations")

simtest <- cbind(data, as.data.frame(sims))

simtest_summ <- simtest %>% group_by(SampleID) %>% 
  summarise(across(.cols = c(PA, V1:V1000), .fns = sum)) %>% 
  pivot_longer(!c(SampleID, PA), names_to = "Sim", values_to = "Prediction")

simtest_quant <- data.frame(x = quantile(simtest_summ$PA, seq(0, 1, 0.0001)),
                            y = quantile(simtest_summ$Prediction, seq(0, 1, 0.0001)))

richness_prediction_plot <- ggplot(simtest_quant, aes(x = x, y = y)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  labs(x = "Observed quantiles", y = "Predicted quantiles",
       title = "Observed vs predicted richness")
