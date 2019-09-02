suppressMessages(require("admixturegraph"))
suppressMessages(require("msm"))
suppressMessages(require("reshape2"))
suppressMessages(require("pscl"))
suppressMessages(require("parallel"))
suppressMessages(require("gridExtra"))
suppressMessages(require(qqman))
suppressMessages(require(data.table))
suppressMessages(require("ggplot2"))
suppressMessages(require("parallel"))
suppressMessages(require("readr"))
suppressMessages(require("graphparse"))
suppressMessages(require("minilexer"))

suppressMessages(require(IRanges))
suppressMessages(require(biomaRt))

# Function to melt table
MeltTab <- function(finaltab){
CHR <- finaltab[,1]
START <- as.numeric(as.character(finaltab[,2]))
END <- as.numeric(as.character(finaltab[,3]))
MIDPOINT <- (START+END)/2
newtab <- finaltab[,seq(4,dim(finaltab)[2])]
newtab <- newtab[,seq(dim(newtab)[2]/2+1,dim(newtab)[2])]
POS <- seq(1,dim(newtab)[1])
newtab <- apply(newtab,2,as.numeric)
melttab <- melt(newtab)
melttab$value <- -log10(melttab$value)
melttab <- cbind(melttab[,c(2,3)], POS,CHR,START,END)
melttab[,1] <- factor(melttab[,1])
names(melttab) <- c("Branches","Pvals","SNPID","CHR","START","END")
melttab <- melttab[order(melttab$Pvals,decreasing=TRUE),]
return(melttab)
}

# Function to select and collapse top regions
CompressTab <- function(tab,extend,padding){
finaltab <- c()
branches <- unique(tab$Branches)
for(branch in branches){
branchtab <- tab[which(tab$Branches == branch),]

# Collapse overlapping regions
branch <- as.character(branchtab$Branches[1])
begin <- as.numeric(as.character(branchtab$START))
end <- as.numeric(as.character(branchtab$END))
start <- begin - extend - padding
stop <- end + extend + padding
pos <- round((begin + end)/2)
chrom <- as.character(branchtab$CHR)
score <- as.numeric(as.character(branchtab$Pvals))
my.df <- data.frame(chrom=chrom,start=start,stop=stop,pos=pos,score=score)
DT <- as.data.table(my.df)
DT[,group := { 
      ir <-  IRanges(start, stop);
       subjectHits(findOverlaps(ir, reduce(ir)))
      },by=chrom]
#collapsed <- DT[, list(start=min(start),stop=max(stop),chrom=unique(chrom),score=max(score)),by=list(group,chrom)]
collapsed <- DT[, list( start=min(start),stop=max(stop),chrom=unique(chrom),score=max(score),bestsnp=pos[which( score == max(score) )] ),by=list(group,chrom)]
#print(collapsed)
# Remove padding
#compressedtab <- cbind(rep(branch,length(collapsed$chrom)),as.character(collapsed$chrom),collapsed$start+padding,collapsed$stop-padding,collapsed$score)
compressedtab <- cbind(rep(branch,length(collapsed$chrom)),as.character(collapsed$chrom),collapsed$start+padding,collapsed$stop-padding,collapsed$bestsnp,collapsed$score)
finaltab <- rbind(finaltab,compressedtab) 
}
#colnames(finaltab) <- c("Branches","CHR","START","END","MAXSCORE")
colnames(finaltab) <- c("BRANCHES","CHR","START","END","BESTPOS","MAXSCORE")
return(finaltab)
}

# Function to annotate top windows with gene names
AnnotateRegions <- function(tab,hostname,datasetname,padding=0){
oldtab <- tab
tab[,3] <- as.numeric(tab[,3]) - padding
tab[,4] <- as.numeric(tab[,4]) + padding
tab <- as.data.table(tab)
oldtab <- as.data.table(oldtab)
enscall = useEnsembl(host=hostname,biomart="ensembl", dataset=datasetname)
genevec <- sapply(seq(1,dim(tab)[1]), function(i){
genes <- getBM(attributes=c('ensembl_gene_id','gene_biotype','hgnc_symbol','chromosome_name','start_position','end_position'), filters = c('chromosome_name','start','end'), values = as.list(tab[i,c(2,3,4)]), mart = enscall)
genes <- genes[which(genes$gene_biotype == "protein_coding"),]
ensembl <- genes$ensembl_gene_id
hgnc <- genes$hgnc_symbol
ensembl <- ensembl[which(ensembl != "")]
hgnc <- hgnc[which(hgnc != "")]
final <- c( paste(ensembl,collapse=","), paste(hgnc,collapse=",") )
return(final)
})
newtab <- cbind(oldtab,t(genevec))
colnames(newtab) <- c("BRANCHES","CHR","START","END","BESTPOS","SCORE","ENSEMBL","HGNC")
newtab[,6] <- round(as.numeric(unlist(newtab[,6])),3)
newtab[,7][which(newtab[,7] == "")] <- "N/A"
newtab[,8][which(newtab[,8] == "")] <- "N/A"
newtab[,1] <- sapply(unlist(newtab[,1]),function(x){paste(strsplit(x,"_")[[1]][-1],collapse="-")})
newtab <- newtab[order(newtab$SCORE,decreasing=TRUE),]
return(newtab)
}

# Function to sort genes for ontology analysis
GetSortedGenes <- function(melttab,hostname,datasetname,corenum,maxormean,extension=0){
ensembl = useEnsembl(host=hostname,biomart="ensembl", dataset=datasetname)
allgenes <- getBM(attributes=c('ensembl_gene_id','gene_biotype','hgnc_symbol','chromosome_name','start_position','end_position'), mart = ensembl)
allgenes <- allgenes[which(allgenes[,2] == "protein_coding" & allgenes[,3] != "" &  allgenes[,4] %in% c("X","Y",seq(1,22)) ),]
duplicates <- which(duplicated(allgenes[,3]))
allgenes <- allgenes[-duplicates,]
allsorted <- mclapply(levels(melttab[,1]),function(branch){
print(branch)
branchtab <- melttab[which(melttab[,1] == branch),]
allscores <- apply(allgenes,1,function(gene){
genechr <- as.character(gene[4])
genestart <- as.numeric(gene[5]) - extension
geneend <- as.numeric(gene[6]) + extension
winchr <- branchtab[,4]
winstart <- branchtab[,5] 
winend <- branchtab[,6]
overlapwin <- which( as.character(winchr) == as.character(genechr) &
((winstart <= genestart & winend >= genestart ) |
(winstart >= genestart & winend <= geneend ) |
(winstart <= geneend & winend >= geneend))
)
if(maxormean == "max"){ if(length(overlapwin) == 0){score <- NA } else{ score <- max(branchtab[overlapwin,2],na.rm=TRUE)} } else {
if(length(overlapwin) == 0){score <- NA } else{ score <- mean(branchtab[overlapwin,2],na.rm=TRUE)} }
return(score)
})
sortedgenes <- cbind(allgenes[,3],allscores)
sortedgenes <- sortedgenes[order(as.numeric(sortedgenes[,2]),decreasing=TRUE),]
sortedgenes <- sortedgenes[ which( sortedgenes[,2] != "NA" & !is.na( sortedgenes[,2] ) ), ]
return(sortedgenes)
},mc.cores=corenum)
names(allsorted) <- levels(melttab[,1])
return(allsorted)
}


# Get test sites for gene ontology analysis
GetTestSites <- function(melttab,pvalcutoff,corenum){
alltestsnps <- mclapply(levels(melttab[,1]),function(branch){
print(branch)
branchtab <- melttab[which(melttab[,1] == branch),]
testsnps <- branchtab[which(branchtab[,2] > pvalcutoff),]
testsnps <- cbind(testsnps[,4],round((testsnps[,5]+testsnps[,6])/2))
return(testsnps)
},mc.cores=corenum)
names(alltestsnps) <- levels(melttab[,1])
return(alltestsnps)
}


# Function for breakin graph into component pieces (from admixture graph package)
break_graph <- function(graph) {
  nodes <- graph$nodes
  edges <- list()
  admixtures <- list()
  parents <- graph$parents
  probs <- graph$probs
  for (i in seq(1, NROW(parents))) {
    match <- which(parents[i, ] == TRUE)
    if (length(match) == 0) {
      root <- nodes[i]
    } else if (length(match) == 1) {
      edges[[length(edges) + 1]] <- c(nodes[i], nodes[match[1]])
    } else if (length(match) == 2) {
      if (nchar(probs[i, match[1]]) > nchar(probs[i, match[2]])) {
        admixtures[[length(admixtures) + 1]] <- c(nodes[i], nodes[match[2]], nodes[match[1]], probs[i, match[2]])
      } else {
        admixtures[[length(admixtures) + 1]] <- c(nodes[i], nodes[match[1]], nodes[match[2]], probs[i, match[1]])
      }
    }
  }
  return(list(leaves = graph$leaves, inner_nodes = graph$inner_nodes, edges = edges, admixtures = admixtures, root = root))
}


# Function for choosing a tree given a particular combination of admixture choices
choosetree <- function(supergraph,chosenedges){
  
  if( length(chosenedges) > 2) {chosenedgestab <- t(matrix(chosenedges,ncol=2))
  } else{chosenedgestab <-  t(matrix(chosenedges))}
  
  test <- break_graph(supergraph[[1]])
  newedgevalues <- as.matrix(supergraph[[2]])
  newedgestab <- t(matrix(unlist(test$edges),nrow=2))
  newinner_nodes <- test$inner_nodes
  newleaves <- test$leaves
  
  # For each chosen path in list of chosen paths
  for(idx in seq(1,dim(chosenedgestab)[1])){
    
    firsthalf <- chosenedgestab[idx,]
    firsthalf.edge <- newedgevalues[which(newedgevalues[,1] == firsthalf[1] & newedgevalues[,2] == firsthalf[2]),3]
    secondhalf <- newedgestab[which(newedgestab[,2]==firsthalf[1]),]
    secondhalf.edge <- newedgevalues[which(newedgevalues[,1] == secondhalf[1] & newedgevalues[,2] == secondhalf[2]),3]
    taped <- c(secondhalf[1],firsthalf[2])
    taped.edge <- as.numeric(firsthalf.edge) + as.numeric(secondhalf.edge)
    
    newedgestab <- newedgestab[which(!(newedgestab[,1] == secondhalf[1] & newedgestab[,2] == secondhalf[2])),]
    newedgestab <- rbind(newedgestab,taped)
    rownames(newedgestab) <- c()
    
    nodetoextract <- secondhalf[2]
    newinner_nodes <- newinner_nodes[which(newinner_nodes != nodetoextract)]
    
    newedgevalues <- newedgevalues[which(!(newedgevalues[,1] == firsthalf[1] & newedgevalues[,2] == firsthalf[2])),]
    newedgevalues <- newedgevalues[which(!(newedgevalues[,1] == secondhalf[1] & newedgevalues[,2] == secondhalf[2])),]
    newedgevalues <- rbind(newedgevalues,c(taped,taped.edge))
    
    loneredge <- newedgevalues[which(newedgevalues[,1] == firsthalf[1] & newedgevalues[,2] != firsthalf[2]),]
    lonernode <- loneredge[2]
    
    firstloner <- newedgestab[which(newedgestab[,1] == lonernode),]
    firstloner.edge <- newedgevalues[which(newedgevalues[,1] == firstloner[1] & newedgevalues[,2] == firstloner[2]),3]
    secondloner <- newedgestab[which(newedgestab[,2] == lonernode),]
    secondloner.edge <- newedgevalues[which(newedgevalues[,1] == secondloner[1] & newedgevalues[,2] == secondloner[2]),3]
    tapedloner <- c(secondloner[1],firstloner[2])
    tapedloner.edge <- as.numeric(firstloner.edge) + as.numeric(secondloner.edge)
    
    newedgestab <- newedgestab[which(!(newedgestab[,1] == firstloner[1] & newedgestab[,2] == firstloner[2])),]
    newedgestab <- newedgestab[which(!(newedgestab[,1] == secondloner[1] & newedgestab[,2] == secondloner[2])),]
    newedgestab <- rbind(newedgestab,tapedloner)
    rownames(newedgestab) <- c()
    
    newinner_nodes <- newinner_nodes[which(newinner_nodes != lonernode)]
    
    newedgevalues <- newedgevalues[which(!(newedgevalues[,1] == firstloner[1] & newedgevalues[,2] == firstloner[2])),]
    newedgevalues <- newedgevalues[which(!(newedgevalues[,1] == secondloner[1] & newedgevalues[,2] == secondloner[2])),]
    newedgevalues <- rbind(newedgevalues,c(tapedloner,tapedloner.edge))
    
    # Remove dangling branch
    newedgevalues <- newedgevalues[which(!(newedgevalues[,1] == secondhalf[2] & newedgevalues[,2] == firstloner[1])),]
    
  }
  newedgevalues <- as.data.frame(newedgevalues)
  edgevec <- as.vector(matrix(t(newedgestab),nrow=1))
  edgevec <- t(matrix(edgevec,nrow=2))
  newedges <- cbind(edgevec,NA)

  newgraph <- admixturegraph::agraph(newleaves, newinner_nodes, newedges, NULL)
  
  supertree <- list(newgraph,newedgevalues)
  return(supertree)
  
}

# Function for extracting embedded trees and embedded probabilities of each tree from a graph
# Returns SUPERTREE object: 1) graph, 2) edge lengths, 3) product of chosen admixture rates, 4) adm path chosen, 5) intermediate paths for all edges
extract_the_trees <- function(supergraph){
  
  brokengraph <- break_graph(supergraph[[1]])
  graphedges <- supergraph[[2]]
  
  if(length(brokengraph$admixtures) > 0){
    
    alladm <- t(matrix((unlist(brokengraph$admixtures)),nrow=4))
    colnames(alladm) <- c("comb","left","right","ratename")
    alladm <- merge(supergraph[[3]],alladm)
    
    admA <- alladm[,c(2,3,4)]
    admB <- alladm[,c(2,3,5)]; admB[,1] <- 1 - admB[,1]
    admA <- apply(admA,1,function(x){paste(x,collapse="_")})
    admB <- apply(admB,1,function(x){paste(x,collapse="_")})
    admAB <- cbind(admA,admB)
    allchoices <- expand.grid(split(admAB, seq(nrow(admAB))))
    
    finaledges <- t(apply(allchoices,1,function(line){
      linetab <- t(matrix(unlist(sapply(line, function(x){
        strsplit(x,"_")
      })),nrow=3))
      finalprob <- prod(as.numeric(linetab[,1]))
      finaledges <-as.vector(t(linetab[,c(2,3)]))
      return(finaledges)
    }))
    
    finalprob <- as.vector(t(apply(allchoices,1,function(line){
      linetab <- t(matrix(unlist(sapply(line, function(x){
        strsplit(x,"_")
      })),nrow=3))
      finalprob <- prod(as.numeric(linetab[,1]))
      return(finalprob)
    })))
    
    
    treelist <- list()
    for(i in seq(1,dim(finaledges)[1])){
      treelist[[length(treelist)+1]] <- choosetree(supergraph,finaledges[i,])
      treelist[[length(treelist)]][[3]] <- finalprob[i]
      treelist[[length(treelist)]][[4]] <- finaledges[i,]
      
      intermediatelist <- apply(treelist[[length(treelist)]][[2]],1,function(x){
        list(GetConnection(x[1],x[2],finaledges[i,],graphedges))
        })
      namesinter <- apply(treelist[[length(treelist)]][[2]],1,function(x){
        paste(x[1],x[2],sep="_")
      })
      names(intermediatelist) <- namesinter
      treelist[[length(treelist)]][[5]] <- intermediatelist
    }
  }
  else{
    treelist <- list()
    treelist[[1]] <- list()
    treelist[[1]][[1]] <- supergraph[[1]]
    treelist[[1]][[2]] <- supergraph[[2]]
    treelist[[1]][[3]] <- 1
    treelist[[1]][[4]] <- NULL
    
    # Collect intermediate paths
    intermediatelist <- apply(treelist[[1]][[2]],1,function(x){
      list(GetConnection(x[1],x[2],finaledges[i,],graphedges))
    })
    namesinter <- apply(treelist[[1]][[2]],1,function(x){
      paste(x[1],x[2],sep="_")
    })
    names(intermediatelist) <- namesinter
    treelist[[1]][[5]] <- intermediatelist
    
  }
  return(treelist)
}


ZeroOne <- function(vecnum){
  return(pmax(0.001,pmin(0.999,vecnum)))
}





DeconstructGraph <- function(supergraph){
  
  graphedges <- supergraph[[2]]
  nameadmixpars <- as.character(supergraph[[3]][,1])
  allbranches <- sapply( seq(1,dim(supergraph[[1]]$parents)[1]), function(x){
    
    child <- rownames(supergraph[[1]]$parents)[x]
    idxchild <- x
    parents <- which(supergraph[[1]]$parents[x,] == TRUE)
    if( length(parents) == 0 ){ return(c(NaN,NaN,NaN,NaN,NaN,NaN,NaN,FALSE))
    } else if( length(parents) == 1){
      parent <- names(parents)
      branch <- graphedges[which(graphedges[,1] == child & graphedges[,2] == parent),3]
      return(c(child,parent,branch,1,NaN,NaN,NaN,FALSE))
    } else{
      parentsrates <- supergraph[[1]]$probs[x,][ which(names(supergraph[[1]]$probs[x, ]) %in% names(parents) ) ]
      recrate <- parentsrates[which(parentsrates %in% nameadmixpars)]
      unrecrate <- parentsrates[which(!(parentsrates %in% nameadmixpars))]
      parentA <- names(recrate)
      parentB <- names(unrecrate)
      rateA <- supergraph[[3]][which(nameadmixpars == recrate),2]
      rateB <- 1 - rateA
      branchA <- graphedges[which(graphedges[,1] == child & graphedges[,2] == parentA),3]
      branchB <- graphedges[which(graphedges[,1] == child & graphedges[,2] == parentB),3]
      finalvec <- c(child,parentA,branchA,rateA,parentB,branchB,rateB,TRUE)
      return(finalvec)
    }
    
  })

  allbranches <- t(matrix(unlist(allbranches),nrow=8))
  allbranches <- allbranches[which(allbranches[,1] != "NaN"),]
  colnames(allbranches) <- c("child","parent","branch","admixrate","parentB","branchB","admixrateB","admevent")
  return(allbranches)
  
}



RemoveFixed <- function(nodes){
    nodes[which(nodes==1)] <- 0.999
    nodes[which(nodes==0)] <- 0.001
    return(nodes)
}



# Get child-ancestor connection possibly involving more than one intermediate parent 
GetConnection <- function(child,parent,finaledges,graphedges){
  
  retrieve <- which(graphedges[,1] == child & graphedges[,2] == parent)
  if(length(retrieve) == 1){
    branches <- cbind(child,parent,as.character(graphedges[retrieve,3]))
  } else{
      tempparent <- NaN
      tempchild <- child
      allnodes <- c(child)
      while(parent != tempparent){
        if( tempchild %in% finaledges[which(seq(1,length(finaledges)) %% 2 == 1)] ){
          tempparent <- finaledges[which(finaledges == tempchild)+1]
          allnodes <- c(allnodes,tempparent)
          tempchild <- tempparent
        } else{
          tempparent <- as.character(graphedges[which(graphedges[,1] == tempchild),2])
          allnodes <- c(allnodes,tempparent)
          tempchild <- tempparent
        }
      }
      branches <- c()
      for( i in seq(1,length(allnodes)-1)){
        retrieve <- which(graphedges[,1] == allnodes[i] & graphedges[,2] == allnodes[i+1])
        branches <- rbind(branches, c(allnodes[i],allnodes[i+1],as.character(graphedges[retrieve,3])))
      }
  }
  return(matrix(branches,ncol=3))
}


min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  return(pmin(f1, f2))
}

integ.min.f1f2 <- function(m1,m2,v1,v2){
  s1 <- sqrt(v1)
  s2 <- sqrt(v2)
  return(integrate(min.f1f2, -Inf, Inf, mu1=m1, mu2=m2, sd1=s1, sd2=s2)$value)
}

shift <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}


AncestralAverage <- function(leaf_vector,supergraph){
  
  recorded <- leaf_vector
  notrecorded <- supergraph[[1]]$inner_nodes
  
  while( length(notrecorded) != 0 ){
    
    testinner <- notrecorded[1]
    rowidx <- which(rownames(supergraph[[1]]$children) == testinner)
    colidx <- which(supergraph[[1]]$children[rowidx,] == TRUE)
    children <- names(colidx)
    
    if( 0 %in% ( children %in% names(recorded)) ){
      notrecorded <- shift(notrecorded)
    } else{
      assigned <- mean(recorded[children])
      recorded <- c(recorded,assigned)
      names(recorded)[length(names(recorded))] <- testinner
      notrecorded <- notrecorded[-1]
    }
    
  }
  
  innerfinal <- recorded[supergraph[[1]]$inner_nodes]
  return(innerfinal)
  
}




CollectBranches <- function(supergraph,targetbranch,deconsgraph,leaves){
  
  graphedges <- supergraph[[2]]
  
  # Start from child
  firstparent <- targetbranch[1]
  
  # Depth-first search of all downstream paths
  tovisit <- c(firstparent)
  branchpaths <- c()
  discovered <- c()
  while(length(tovisit) != 0){
    check = tovisit[length(tovisit)]
    tovisit = tovisit[-length(tovisit)]
    
    
    if (!(check %in% discovered)){
      
      discovered <- c(discovered,check)
      daughters <-  as.character(graphedges[graphedges$parent == check,1])
      for(daughter in daughters){
        
        daughtervec <- deconsgraph[which(deconsgraph[,1] == daughter),]
        
        
        if(daughtervec[2] == check){ prob <- as.numeric(daughtervec[4])
        } else if(daughtervec[5] == check){ prob <- as.numeric(daughtervec[7])}
        
        path <- c(daughter,check,prob)
        branchpaths <- rbind(branchpaths,path)
        
      }
      if(length(daughters) > 0){
        tovisit <- c(tovisit,daughters)
      }
    }
  }
  

  
  prevpaths <- matrix(NaN)
  # Connect all possible paths from leaves to the starting parent
  #for(rep in seq(1,10)){
  while(TRUE){  
    for(pathidx in which(branchpaths[,1] %in% leaves)){
      path <-  branchpaths[pathidx,]
      for(conidx in which(branchpaths[,1] == path[2])){
        con <- branchpaths[conidx,]
        branchpaths <- rbind(branchpaths, c( path[1], con[2], as.numeric(path[3]) * as.numeric(con[3]) )  )
      }
    }
    
    branchpaths <- unique(branchpaths)
    currpaths <- branchpaths
    #print(branchpaths)
    #print(dim(currpaths))
    #print(dim(prevpaths))
    if( prod(dim(currpaths) == dim(prevpaths))){break}
    prevpaths <- currpaths
    
  }
  

  # Only keep paths that go from a leaf to the starting parent
  branchpaths <- unique(branchpaths[which(branchpaths[,1] %in% leaves & branchpaths[,2] == firstparent),])
  
  leafcontrib <- c()
  
  # Get first adm probability
  candidate <- which(deconsgraph[,1] == targetbranch[1] & deconsgraph[,2] == targetbranch[2])
  if(length(candidate) == 1){ 
    firstprod <- as.numeric(deconsgraph[candidate,4])
  } else{
    candidate <- which(deconsgraph[,1] == targetbranch[1] & deconsgraph[,5] == targetbranch[2])
    firstprod <- as.numeric(deconsgraph[candidate,7])
  }
  if(is.vector(branchpaths)){ branchpaths <- t(matrix(branchpaths))}
  
  
  # Compute vectors of leaf-contributions
  if(firstparent %in% leaves){
    leafcontrib <- as.numeric(leaves == firstparent)
  } else {
    for(leaf in leaves){
      if( !(leaf %in% branchpaths[,1]) ){
        leafcontrib <- c(leafcontrib,0)
      } else{
        
        tempcontrib <- c()
        i <- 1
        for( currparent in branchpaths[which(branchpaths[,1] == leaf),2]){
          currchild <- leaf
          prod <- firstprod * as.numeric(branchpaths[which(branchpaths[,1] == currchild & branchpaths[,2] == currparent),3][i])
          tempcontrib <- c(tempcontrib,prod)
          i <- i + 1
        }
        
        leafcontrib <- c(leafcontrib,sum(tempcontrib))
      }
    }
  }
  
  return(leafcontrib)
  
}


OrderBranches <- function(supergraph){

  oldedgevalues <- supergraph[[2]]
  root <- as.character(unique(supergraph[[2]]$parent[!(supergraph[[2]]$parent %in% supergraph[[2]]$child)]))
  newedgevalues <- c()
  parents <- c(root)
  
  while( length(parents) < length(supergraph[[1]]$nodes) ){

    children_to_test <- supergraph[[1]]$nodes[which( !(supergraph[[1]]$nodes %in% parents) )]
    
    for( child in children_to_test ){
      child <- as.character(child)
      
      allparents <- names(which(supergraph[[1]]$parents[child,]))
      
      if ( prod(names(which(supergraph[[1]]$parents[child,])) %in% parents) == 1){
        toadd <- oldedgevalues[oldedgevalues$child == child,]
        newedgevalues <- rbind(newedgevalues,toadd)
        parents <- c(parents,child)
      }
      parents <- unique(parents)
    }
  }
  
  newedgevalues <- data.frame(newedgevalues)
  colnames(newedgevalues) <- c("child","parent","value")
  
  #print(newedgevalues)
  
  supergraph <- list(supergraph[[1]],newedgevalues,supergraph[[3]])
  
  return(supergraph)
}


norm_vec <- function(x) sqrt(sum(x^2))





ObtainFreqs <- function(countdat,fcutoff){
  dercounts <- apply(countdat,c(1,2),function(x){splitted <- strsplit(x,",")[[1]]; return( as.numeric(splitted[2]) )})
  totalcounts <- apply(countdat,c(1,2),function(x){splitted <- strsplit(x,",")[[1]]; return( as.numeric(splitted[2])+as.numeric(splitted[1]) )})
  dersum <- apply(dercounts,1,function(x){sum(x)})
  totalsum <- apply(totalcounts,1,function(x){sum(x)})
  totalfreq <- dersum/totalsum
  freqs <- apply(countdat,c(1,2),function(x){splitted <- strsplit(x,",")[[1]]; return( as.numeric(splitted[2]) / (as.numeric(splitted[2])+as.numeric(splitted[1])) )})
  checksegneut <- which(totalfreq > fcutoff & totalfreq < (1-fcutoff) ) 
  freqs <- freqs[checksegneut,]
  return( list( as.matrix(freqs), checksegneut ) )
}

# Trim out white space
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


# Modified for single-locus selection
LoadCounts <- function(filename,pops){
  table <- as.matrix(read.table(filename,header=TRUE,sep="\t",strip.white=TRUE))
  tokeep <- which(colnames(table) %in% pops)
  tokeep <- c(1,2,3,tokeep)
  table <- table[,tokeep]
  table[,1] <- trim(table[,1])
  table[,2] <- trim(table[,2])
  noNA <- which(apply(table, 1, function(line){ !("0,0" %in% line) }))
  table <- table[noNA,]
  return(table)
}



# Correct positions
CorrectPos <- function(inputvals,coords,listvertices,root,minbranch=0.075){
  
  # Use only one parent branch when dealing with admixed pops
  inputvals <- inputvals[!duplicated(inputvals[,1]),]
  
  todo <- as.matrix(inputvals[inputvals[,2] == root,])
  done <- rbind( c(root,coords[listvertices == root,]) )
  
  
  while(length(todo) > 0){
    
    if( !(is.vector(todo)) ){
      elem <- todo[1,1]
      parent <- todo[1,2]
    } else {
      elem <- todo[1]
      parent <- todo[2]
    }
    
    elemx <- coords[which(listvertices == elem),1]
    elemy <- coords[which(listvertices == elem),2]
    parenty <- as.numeric(done[which(done[,1] == parent),3])
    minbranch <- as.numeric(minbranch)
    addy <- max(minbranch,inputvals[which(as.character(inputvals[,1]) == elem & as.character(inputvals[,2]) == parent),3])
    addy <- as.numeric(addy)
    newy <- parenty - addy
    
    done <- rbind(done, c(elem,elemx,newy))
    
    if( !(is.vector(todo)) ){
      todo <- todo[-1,]
    }
    else {
      todo <- c()
    }
    
    if(sum(as.character(inputvals[,2]) == elem) > 0){
      toadd <- inputvals[as.character(inputvals[,2]) == elem,]
      toadd <- sapply(toadd,as.character)
      
      if(is.vector(toadd)){
        todo <- rbind(todo,toadd)
      } else{ 
        for(j in seq(1,dim(toadd)[1])){
          todo <- rbind(todo,toadd[j,])
        }
      }
      
    }
  }
  
  
  done <- done[!duplicated(done[,1]),]
  done <- done[match(listvertices,done[,1]),]
  
  final <-  t(apply(done[,c(2,3)],1,as.numeric))
  return(final)
  
}



ChiSquaredReduced <- function(graphedges,contribmat,Fmat,leaves_freqs,effects,total=FALSE,randomize=FALSE){

  checkseg <- which( apply(leaves_freqs,1,sum)/dim(leaves_freqs)[2] < 0.99  & apply(leaves_freqs,1,sum)/dim(leaves_freqs)[2] > 0.01 )
  leaves_freqs <- leaves_freqs[checkseg,]
  effects <- effects[checkseg]

  # Randomize effects if necessary
  if(randomize == TRUE){effects <- effects * sample(c(-1,1),length(effects),replace=TRUE)}

  # Compute mean genetic values
  meangen <- apply(leaves_freqs * effects, 2, function(x){sum(x)})

  # Scale by average genetic value
  meangen <- (meangen - mean(meangen))

  # Compute the estimated ancestral genetic variance over all populations
  meanfreqs <- apply(leaves_freqs,1,mean)

  varmean <- sum(sapply(seq(1,length(meanfreqs)),function(i){
  score = meanfreqs[i]*(1-meanfreqs[i])*effects[i]^2
  return(score)
  }))

  if(total == FALSE){
  contribmat <- contribmat[,names(meangen)]
  }

  if(total == FALSE){
  # Compute Q_B test statistic
  Qteststat <- apply(contribmat,1, function(x){
    #print(varmean)
    lambda <- varmean * (t(x) %*% Fmat %*% x)
    numerator <- (meangen %*% x)^2
    final <- numerator / (lambda)
    return(final)
  } )

  # Compute q_B test statistic
  qteststat <- apply(contribmat,1, function(x){
    #print(varmean)
    lambda <- varmean * (t(x) %*% Fmat %*% x)
    numerator <- (meangen %*% x)
    final <- numerator / sqrt(lambda)
    return(final)
  } )

  branchorder <- apply(graphedges,1,function(x){paste(as.character(x[1]),as.character(x[2]),sep="_")})
  Qteststat <- Qteststat[branchorder]
  Pval <- max( 1.110223e-16, 1 - pchisq(Qteststat,1))
  qteststat <- qteststat[branchorder]
  allstats <- cbind(Qteststat,qteststat,Pval)

  }

  if(total == TRUE){

    # Compute Q_X statistic
    numerator <- t(meangen) %*% solve(Fmat) %*% meangen
    denominator <- varmean
    Qteststat <- numerator / denominator

    Pval <-  max( 1.110223e-16, 1 - pchisq(Qteststat,qr(Fmat)$rank))
    allstats <- c(Qteststat,NaN,Pval)
  }

  return(allstats)
}


ComputeMeanFreq <- function(freqs,invFmat,mode){
if(mode == "A"){
	# Simple mean
	return(mean(freqs))
} else{
	# F-weighted mean
	#print(freqs)
	ones <- rep(1,length(freqs))
	numerator <- t(ones) %*% invFmat %*% freqs
	#print(numerator)
	denominator <- t(ones) %*% invFmat %*% ones
	#print(denominator)
	#print(numerator / denominator)
	Fmean <- c(numerator / denominator)
	Fmean <- min(max(Fmean,0.05),0.95)
	return(Fmean)
}
}



ComputeWinRB <- function(branchorder,contribmat,Fmat,invFmat,freqs){

  contribmat <- contribmat[,colnames(Fmat)]

  if( is.null(dim(freqs)) ){
    freqs <- freqs[colnames(Fmat)]
    ancfreqs <- ComputeMeanFreq(freqs,invFmat,"A")
    #ancfreqs <- ComputeMeanFreq(freqs,invFmat,"B")
    freqs <- freqs - ancfreqs
    anchet <- ancfreqs*(1-ancfreqs)
    sumfreqs <- freqs

  } else{
    freqs <- freqs[,colnames(Fmat)]
    ancfreqs <- apply(freqs,1,mean)
    freqs <- freqs - ancfreqs
    anchet <- ancfreqs*(1-ancfreqs)
    sumfreqs <- apply( freqs, 2, sum)
  }

  #print(Fmat)
  #print(freqs)
  #print(contribmat)
  #print(ancfreqs)
  #print(anchet)
  #print(sumfreqs)  

  sumanchet <- sum(anchet)


  teststat <- apply(contribmat,1, function(x){
    #print(varmean)
    lambda <- sumanchet * (t(x) %*% Fmat %*% x)
    numerator <- (sumfreqs %*% x)^2
    final <- numerator / (lambda)
    return(final)
  } )

  teststat <- teststat[branchorder]

  Pval <- sapply( teststat, function(y){ max( 1.110223e-16, 1 - pchisq(y,1)) })

  names(Pval) <- paste("Pval_",names(teststat),sep="")

  return(c(teststat,Pval))
}


get_edges <- function(g, node_name) {
incidence_vector <- g$parents[node_name,]
parents <- names(incidence_vector)[which(incidence_vector == TRUE)]
if(length(parents) == 0){ return()
} else if(length(parents) == 1){ return(c(node_name,parents))
} else if(length(parents) > 1){ return( rbind( rep(node_name,length(parents)) ,parents)  ) }
}


