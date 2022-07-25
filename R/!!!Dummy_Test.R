rm(list=ls())
setwd("~/Desktop/Programming/R_Root/PhD_Project/DRP_PRISM/")

library(igraph)

N_node = 1000
g = make_ring(N_node, circular = TRUE)
plot(g)

# g = make_star(N_node, mode = "undirected")
# plot(g)

edge_index1 = t(as_edgelist(g))
edge_index1 = edge_index1-1
edge_index = cbind(edge_index1, edge_index1[c(2,1),])


node_attr = cbind(matrix(rnorm(length(E(g))*N_graph, .1,.001),length(E(g)),1),
                   matrix(rnorm(length(E(g))*N_graph, .9,.001),length(E(g)),1))
#edge_attr = rep(1,2*N_node)
graph_label = c(rep(0,N_graph),rep(1,N_graph))

#Save Data
write.table(edge_index, file = "Data",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

write.table(edge_attr, file = "Data_Complete/Saved_Data/Test/edge_attrs.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

write.table(graph_label, file = "Data_Complete/Saved_Data/Test/graph_labels.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(node_attr, file = "Data_Complete/Saved_Data/Test/node_attrs.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")