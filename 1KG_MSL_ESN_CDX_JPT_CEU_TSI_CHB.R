library("admixturegraph")

# Define standard admixture graph
leaves <- c("MSL", "ESN", "TSI","CEU","CDX","JPT","CHB")
inner_nodes <- c("r","q","t","u","s","v")
edges <- parent_edges(c(edge("CEU","v"),
                        edge("TSI","v"),
                        edge("v","q"),
                        edge("CHB","u"),
                        edge("JPT","u"),
                        edge("u","t"),
                        edge("CDX","t"),
                        edge("t","q"),
                        edge("q","r"),
                        edge("s","r"),
                        edge("MSL","s"),
                        edge("ESN","s")))
graph <- agraph(leaves, inner_nodes, edges,NULL)



# Define admixture values
admvalues <- c()
if(length(admvalues) > 0){
  colnames(admvalues) <- c("ratename","rate")
  admvalues <- as.data.frame(admvalues)
  admvalues[,2] <- as.numeric(as.character(admvalues[,2]))
}


edgevalues <- unique(cbind(edges[,1],edges[,2],1))
colnames(edgevalues) <- c("child","parent","value")
edgevalues <- as.data.frame(edgevalues)
edgevalues[,3] <- as.numeric(as.character(edgevalues[,3]))
supergraph <- list(graph,edgevalues,admvalues)

# Order edges topologically
supergraph <- OrderBranches(supergraph)

source("LoadFiles.R")
