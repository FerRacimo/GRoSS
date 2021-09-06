# Define standard admixture graph
print("Loading graph topology")
if(is.null(dotfile)){ print("Reading *graph file"); pregraph <- readr::read_file(graphfile); graph <- graphparse::read_qpgraph(pregraph)
} else { print("Reading *dot file"); pregraph <- readr::read_file(dotfile); graph <- graphparse::read_dot(pregraph) }
vecadm <- attr(graph, "admixture_proportions")
if(is.null(vecadm)){ admvalues <- c() } else{ 
	admvalues <- cbind(names(vecadm),matrix(vecadm))}
leaves <- graph$leaves
inner_nodes <- graph$inner_nodes
edges <- sapply(c(leaves,inner_nodes), function(x){get_edges(graph, x)})
edges <- matrix(unlist(edges),byrow=TRUE,ncol=2)

if(length(admvalues) > 0){
  colnames(admvalues) <- c("ratename","rate")
  admvalues <- as.data.frame(admvalues)
  admvalues[,2] <- as.numeric(as.character(admvalues[,2]))
}

edgevalues <- unique(cbind(edges,1))
colnames(edgevalues) <- c("child","parent","value")
rownames(edgevalues) <- c()
edgevalues <- as.data.frame(edgevalues)
edgevalues[,3] <- as.numeric(as.character(edgevalues[,3]))
supergraph <- list(graph,edgevalues,admvalues)

# Order edges topologically
supergraph <- OrderBranches(supergraph)

if(exists("neutfile")){
    
  # Load Neutral data
  print("Loading SNP data")
  neutfilename <- neutfile
  neutdata <- LoadCounts(neutfilename, leaves)
  firstfreqcol <- 4
  neut_leaves_counts <- as.data.frame(neutdata[,seq(firstfreqcol,dim(neutdata)[2])])
  fcutoff <- 0.05
  raw_freqs <- ObtainFreqs(neut_leaves_counts,fcutoff)
  neut_leaves_freqs <- raw_freqs[[1]]
  checksegneut <- raw_freqs[[2]]
  leaves_finalcounts <- raw_freqs[[3]]
  snpinfo <- neutdata[,c(1,2,3)]
  snpinfo <- snpinfo[checksegneut,]
  
  # Deconstruct graph
  graphedges <- supergraph[[2]]
  deconsgraph <- DeconstructGraph(supergraph)
  leaves <- supergraph[[1]]$leaves

  # Filter for segregating sites
  print("Computing F matrix...")
  #checksegneut <- which( apply(neut_leaves_freqs,1,sum)/dim(neut_leaves_freqs)[2] < (1 - fcutoff)  & apply(neut_leaves_freqs,1,sum)/dim(neut_leaves_freqs)[2] > fcutoff )
  #print(checksegneut)
  #neut_leaves_freqs <- neut_leaves_freqs[checksegneut,]


  # Compute empirical covariance matrix
  #checkLG <- which(snpinfo[,1] != "LG01" & snpinfo[,1] != "LG02" & snpinfo[,1] != "LG12" & snpinfo[,1] != "LG07")
  #snpinfo_cov <- snpinfo[checkLG,]
  #neut_leaves_freqs_cov <- neut_leaves_freqs[checkLG,]
  snpinfo_cov <- snpinfo
  neut_leaves_freqs_cov <- neut_leaves_freqs
  neut_leaves_freqs_means <- apply(neut_leaves_freqs_cov, 1, mean)	  
  mean_hetero <- neut_leaves_freqs_means*(1-neut_leaves_freqs_means)
  Fmat <- sapply(seq(1,dim(neut_leaves_freqs_cov)[2]),function(x){
       sapply(seq(1,dim(neut_leaves_freqs_cov)[2]),function(y){
	cov(neut_leaves_freqs_cov[,x]/sqrt(mean_hetero), neut_leaves_freqs_cov[,y]/sqrt(mean_hetero))
	})
       })
  colnames(Fmat) <- colnames(neut_leaves_freqs_cov)
  rownames(Fmat) <- colnames(neut_leaves_freqs_cov)
  #print(Fmat)

  invFmat <- solve(Fmat)

  # Compute contributions of each branch to each leaf
  contribmat <- c()
  contribmat_names <- c()
  for(branchidx in seq(1,dim(graphedges)[1])){
    targetbranch <- c(as.character(graphedges[branchidx,1]),as.character(graphedges[branchidx,2]))
    name <- paste(targetbranch,collapse="_")
    branchvec <- CollectBranches(supergraph,targetbranch,deconsgraph,leaves)
    contribmat <- rbind(contribmat,branchvec)
    contribmat_names <- rbind(contribmat_names,name)
  }
  rownames(contribmat) <- contribmat_names
  colnames(contribmat) <- leaves
  branchorder <- apply(graphedges,1,function(x){paste(as.character(x[1]),as.character(x[2]),sep="_")})
  
  # Standardize vectors to be of unit length
  contribmat <- t(apply(contribmat,1, function(x) {
    x <- x - mean(x)
    return(x / norm_vec(x))
  }))


  chrvec <- unique(snpinfo[,1])

  finalstats <- do.call(rbind, lapply(chrvec, function(chr){

    print(chr)
    chrset <- which(snpinfo[,1] == chr)
    setsnpinfo <- snpinfo[chrset,]
    setfreqs <- neut_leaves_freqs[chrset,]
    setcounts <- leaves_finalcounts[chrset,]

    numSNPs <- dim(setfreqs)[1]

    winsize <- 1
    winshift <- 1
    SNPrange <- seq(1,numSNPs-winsize,winshift)

    start <- setsnpinfo[SNPrange,2]
    end <- setsnpinfo[SNPrange+(winsize-1),2]
    chrlist <- setsnpinfo[SNPrange,1]

    allstats <- t(sapply(seq(1,numSNPs-winsize,winshift), function(x) {
    #print(c(x,x+(winsize-1)))
    SNPs <- seq(x,x+(winsize-1),1)
    freqs <- setfreqs[SNPs,]
    totcounts <- setcounts[SNPs,]
    print(freqs)
    print(totcounts)
    stats <- ComputeWinRB(branchorder,contribmat,Fmat,invFmat,freqs,finitesamp,totcounts)
    return(stats)
    } ))

    allstats <- cbind(chrlist,start,end,allstats)
    colnames(allstats) <- c("CHR","START","END",branchorder, paste("Pval_",branchorder,sep=""))
    return(allstats)
  }))

  #print(finalstats)
  write.table(finalstats, file=outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

}
