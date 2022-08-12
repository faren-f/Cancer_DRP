rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
library(OmnipathR)
# library(tidyr)
# library(dnet)
# library(gprofiler2)
# library(magrittr)
# library(dplyr)


# import Omnipath protein-protein interactions
interactions = import_omnipath_interactions() %>% as_tibble()
interaction_gene_symb = interactions[,c(3,4)]

saveRDS(interaction_gene_symb, "Processed_data/Step8/ppi_Omnipath.rds")



# Convert to igraph objects:
#OPI_g = interaction_graph(interactions = interactions )

#???interactions = import_pathwayextra_interactions(resources=c("BioGRID","STRING"),
                                                #organism = 10090)

# ppi <-
#   import_omnipath_interactions(
#     datasets = 'omnipath',
#     entity_types = 'protein'
#   )

## ----network------------------------------------------------------------------------------------------------
#gri <- import_transcriptional_interactions()
## ----network-igraph-----------------------------------------------------------------------------------------
#gr_graph <- interaction_graph(gri) 
 

## ----go-----------------------------------------------------------------------------------------------------
# go <- go_ontology_download()
# go$rel_tbl_c2p
# ## ----go-graph-----------------------------------------------------------------------------------------------
# go_graph <- relations_table_to_graph(go$rel_tbl_c2p)
# ## ----go-name------------------------------------------------------------------------------------------------
# ontology_ensure_name('GO:0000022')




