library(tidyverse)
library(MCMCglmm)

## Format ASV table
asv_table <- qiime2R::read_qza("QIIME-final/table.qza")
asv_table <- t(asv_table$data) %>% as.data.frame() %>% filter(!rowSums(.) < 50)
asv_table$TotReads <- rowSums(asv_table)
asv_table <- asv_table %>% rownames_to_column("SampleID") %>%
  pivot_longer(!c(SampleID, TotReads), names_to = "ASV", values_to = "Reads") %>%
  ## Convert to presence/absence
  mutate(pres_abs = ifelse(Reads > 0, 1, 0))

# ## Format metadata table, impute missing values for completeness
metadata <- qiime2R::read_q2metadata("QIIME-intermediate/sample-metadata.tsv") %>%
  select(SampleID, `site-name`, elevation, depth:percentcover)
names(metadata) <- gsub("-", ".", names(metadata))
metadata_impute <- complete(mice(metadata)) # Impute missing values
metadata_impute <- metadata_impute %>%
  mutate(transect = substr(site.name, 1, 3)) # Create factor for each transect

# ## Get taxonomy table and filter to retain only those with taxonomic
# ## features of interest - i.e. all must have a phylum
taxonomy <- qiime2R::read_qza("QIIME-final/taxonomy.qza")
taxonomy <- taxonomy$data %>%
  separate(Taxon, sep = "; ",
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  rename(ASV = Feature.ID)
taxonomy[is.na(taxonomy$Phylum),]$Phylum <- "Unknown"
taxonomy[is.na(taxonomy$Family),]$Family <- "Unknown"

# ## Join data together and scale variables
glmm_dataset <- left_join(asv_table, metadata_impute) %>%
  mutate(elevation_mc = scale(elevation)[,1],
         percentcover_mc = scale(percentcover)[,1],
         average.soil.temperature_mc = scale(average.soil.temperature)[,1],
         average.soil.relative.humidity_mc = scale(average.soil.relative.humidity)[,1],
         ph_mc = scale(ph)[,1],
         toc_mc = scale(toc)[,1],
         ec_mc = scale(ec)[,1],
         transect = substr(site.name, 1, 3)) %>%
  left_join(taxonomy)

write_csv(glmm_dataset, "glmm_dataset.csv")

k <- 1000

mod <- MCMCglmm(pres_abs ~ elevation_mc + percentcover_mc + average.soil.temperature_mc + 
                  average.soil.relative.humidity_mc + ec_mc + TotReads,
                random = ~ us(1 + elevation_mc + percentcover_mc + average.soil.temperature_mc + 
                                average.soil.relative.humidity_mc + ec_mc):ASV,
                family = "categorical",
                data = glmm_dataset, 
                nitt = 500000, burnin = 60000,
                prior = list(R=list(V=0.5,fix=1),
                             G=list(G1=list(V=diag(6),nu=1,alpha.mu=rep(0, 6),alpha.V=diag(6)*k))),
                thin = 100,
                pr = T)

save(mod, file = "Model_noPhylo.Rdata")


mod2 <- MCMCglmm(pres_abs ~ elevation_mc + TotReads, 
                 random = ~ us(1 + elevation_mc):Phylum + 
                   us(1 + elevation_mc):ASV,
                family = "categorical",
                data = glmm_dataset, 
                nitt = 200000, burnin = 60000,
                prior = list(R=list(V=0.5,fix=1),
                             G=list(G1=list(V=diag(2),nu=1,alpha.mu=rep(0, 2),alpha.V=diag(2)*k),
                                    G1=list(V=diag(2),nu=1,alpha.mu=rep(0, 2),alpha.V=diag(2)*k))),
                thin = 100,
                pr = T)
save(mod2, file = "Model_Phylo.Rdata")
