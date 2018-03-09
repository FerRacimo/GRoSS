
if(exists("neutfile")){
    
  # Load Neutral data
  neutfilename <- neutfile
  neutdata <- LoadCounts(neutfilename, leaves)
  firstfreqcol <- 4
  neut_leaves_counts <- as.data.frame(neutdata[,seq(firstfreqcol,dim(neutdata)[2])])
  neut_leaves_freqs <- ObtainFreqs(neut_leaves_counts)
  snpinfo <- neutdata[,c(1,2,3)]

  # Deconstruct graph
  graphedges <- supergraph[[2]]
  deconsgraph <- DeconstructGraph(supergraph)
  leaves <- supergraph[[1]]$leaves

  # Compute empirical covariance matrix
  print("Computing F matrix...")
  fcutoff <- 0.01
  checksegneut <- which( apply(neut_leaves_freqs,1,sum)/dim(neut_leaves_freqs)[2] < (1 - fcutoff)  & apply(neut_leaves_freqs,1,sum)/dim(neut_leaves_freqs)[2] > fcutoff )
  neut_leaves_freqs <- neut_leaves_freqs[checksegneut,]
  neut_leaves_freqs_means <- apply(neut_leaves_freqs, 1, mean)	
  mean_hetero <- neut_leaves_freqs_means*(1-neut_leaves_freqs_means)
  snpinfo <- snpinfo[checksegneut,]
  Fmat <- sapply(seq(1,dim(neut_leaves_freqs)[2]),function(x){
       sapply(seq(1,dim(neut_leaves_freqs)[2]),function(y){
	cov(neut_leaves_freqs[,x]/sqrt(mean_hetero), neut_leaves_freqs[,y]/sqrt(mean_hetero))
	})
       })
  colnames(Fmat) <- colnames(neut_leaves_freqs)
  rownames(Fmat) <- colnames(neut_leaves_freqs)

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

    numSNPs <- dim(setfreqs)[1]
    #numSNPs <- 200

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
    stats <- ComputeWinRB(branchorder,contribmat,Fmat,freqs)
    return(stats)
    } ))

    allstats <- cbind(chrlist,start,end,allstats)
    colnames(allstats) <- c("CHR","START","END",branchorder, paste("Pval_",branchorder,sep=""))
    return(allstats)
  }))

  #print(finalstats)
  write.table(finalstats, file=outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

}
