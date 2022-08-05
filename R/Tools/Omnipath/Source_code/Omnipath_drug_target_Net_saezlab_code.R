rm(list=ls())
## ----message=FALSE, warning=FALSE---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(OmnipathR)
library(igraph)
library(ggraph)

## -----------------------------------------------------------------------------------------------------------
# Download protein-protein interactions
interactions = import_omnipath_interactions() %>% as_tibble()

# Convert to igraph objects:
OPI_g = interaction_graph(interactions = interactions )

## ---- eval=FALSE--------------------------------------------------------------------------------------------
#  library(dbparser)
#  library(XML)
#  
#  
#  ## parse data from XML and save it to memory
#  get_xml_db_rows("..path-to-DrugBank/full database.xml")
#  
#  ## load drugs data
 drugs <- parse_drug() %>% select(primary_key, name)
 drugs <- rename(drugs,drug_name = name)

#  ## load drug target data
 drug_targets <- parse_drug_targets() %>%
    select(id, name,organism,parent_key) %>%
    rename(target_name = name)
#  
#  ## load polypeptide data
#  drug_peptides <- parse_drug_targets_polypeptides()  %>%
#     select(id, name, general_function, specific_function,
#            gene_name, parent_id) %>%
#     rename(target_name = name, gene_id = id)
#  
#  # join the 3 datasets
#  drug_targets_full <- inner_join(drug_targets, drug_peptides,
#                                  by=c("id"="parent_id", "target_name")) %>%
#     inner_join(drugs, by=c("parent_key"="primary_key")) %>%
#     select(-other_keys)
#  

## -----------------------------------------------------------------------------------------------------------
drug_names = c("Valproat"      = "Valproic Acid",
               "Diclofenac"    = "Diclofenac",
               "Paracetamol"   = "Acetaminophen",
               "Ciproflaxin"   = "Ciprofloxacin",
               "Nitrofurantoin"= "Nitrofurantoin",
               "Tolcapone",
               "Azathioprine",
               "Troglitazone",
               "Nefazodone",
               "Ketoconazole",
               "Omeprazole",
               "Phenytoin",
               "Amiodarone",
               "Cisplatin",
               "Cyclosporin A"  = "Cyclosporine",
               "Verapamil",
               "Buspirone",
               "Melatonin",
               "N-Acetylcysteine"= "Acetylcysteine",
               "Vitamin C"       = "Ascorbic acid",
               "Famotidine",
               "Vancomycin")

## ---- eval=FALSE--------------------------------------------------------------------------------------------
#  
#  drug_target_data_sample <- drug_targets_full %>%
#     filter(organism == "Humans",drug_name %in% drug_names)
#  

## -----------------------------------------------------------------------------------------------------------
drug_targets <- OmnipathR:::drug_target_data_sample %>%
  filter(organism == "Humans",drug_name %in% drug_names)

## -----------------------------------------------------------------------------------------------------------
drug_targets <-  drug_targets %>%
  select(-target_name, -organism) %>%
  mutate(in_OP = gene_id %in% c(interactions$source))
# not all drug-targets are in OP.
print(all(drug_targets$in_OP))

# But each drug has at least one target in OP.
drug_targets %>% group_by(drug_name) %>% summarise(any(in_OP))


## -----------------------------------------------------------------------------------------------------------
POI = tibble(protein = c("NFE2L2","HMOX1","TP53","CDKN1A","BTG2","NFKB1",
                         "ICAM1","HSPA5", "ATF4","DDIT3","XBP1"))

## -----------------------------------------------------------------------------------------------------------
POI <- POI %>% mutate(in_OP = protein %in% interactions$target_genesymbol)
# all POI is in Omnipath
print(all(POI$in_OP))



## -----------------------------------------------------------------------------------------------------------

source_nodes <- drug_targets %>%
  filter(in_OP, drug_name=="Cisplatin") %>%
  pull(gene_name)
target_nodes <- POI %>% filter(in_OP) %>% pull(protein)

collected_path_nodes = list()

for(i_source in 1:length(source_nodes)){
  
  paths <- shortest_paths(OPI_g, from = source_nodes[[i_source]],
                          to = target_nodes,
                          output = 'vpath')
  path_nodes <- lapply(paths$vpath,names) %>% unlist() %>% unique()
  collected_path_nodes[[i_source]] <- path_nodes
}
collected_path_nodes <- unlist(collected_path_nodes) %>% unique()

## -----------------------------------------------------------------------------------------------------------
cisplatin_nodes <- c(source_nodes,target_nodes, collected_path_nodes) %>%
  unique()
cisplatin_network <- induced_subgraph(graph = OPI_g,vids = cisplatin_nodes)
a = as_data_frame(cisplatin_network)

## -----------------------------------------------------------------------------------------------------------
V(cisplatin_network)$node_type = ifelse(
  V(cisplatin_network)$name %in% source_nodes, "direct drug target",
  ifelse(
    V(cisplatin_network)$name %in% target_nodes,"POI","intermediate node"))

ggraph(
  cisplatin_network,
  layout = "lgl",
  area = vcount(cisplatin_network)^2.3,
  repulserad = vcount(cisplatin_network)^1.2,
  coolexp = 1.1
) +
  geom_edge_link(
    aes(
      start_cap = label_rect(node1.name),
      end_cap = label_rect(node2.name)),
    arrow = arrow(length = unit(4, 'mm')
    ),
    edge_width = .5,
    edge_alpha = .2
  ) +
  geom_node_point() +
  geom_node_label(aes(label = name, color = node_type)) +
  scale_color_discrete(
    guide = guide_legend(title = 'Node type')
  ) +
  theme_bw() +
  xlab("") +
  ylab("") +
  ggtitle("Cisplatin induced network")


## ---- sessionInfo, echo=FALSE-------------------------------------------------------------------------------
sessionInfo()
