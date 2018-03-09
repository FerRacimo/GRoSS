library("admixturegraph")

# Define standard admixture graph
leaves <- c("Mandenka", "Yoruba", "EuropeA", "Oceanian", "EastAsian", "NativeAmerican", "Sardinian")
inner_nodes <- c("r","q","t","w","s","v","x","y")
edges <- parent_edges(c(edge("t","r"),
                        edge("q","r"),
			edge("Yoruba","t"),
			edge("Mandenka","t"),
			edge("s","q"),
			edge("y","q"),
			edge("Sardinian","y"),
			edge("w","y"),
			edge("Oceanian","s"),
			edge("v","s"),
                        edge("EastAsian","v"),
			edge("x","v"),
			edge("NativeAmerican","x"),
			edge("w","x"),
			edge("EuropeA","w"),
                        admixture_edge("w","y","x")))
admixtures <- admixture_proportions(c(
  admix_props("w", "y", "x", "gamma")
))
graph <- agraph(leaves, inner_nodes, edges, admixtures)

# Define admixture values
admvalues <- rbind(c("gamma",0.89))

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