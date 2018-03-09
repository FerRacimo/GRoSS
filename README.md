# GRoSS
Graph-aware Recovery of Selective Sweeps


latexImg = function(latex){

    link = paste0('http://latex.codecogs.com/gif.latex?',
           gsub('\\=','%3D',URLencode(latex)))

    link = gsub("(%..)","\\U\\1",link,perl=TRUE)
    return(paste0('![](',link,')'))
}


We modified the $Q_{B}$ statistics (Racimo et al. 2018) to handle single-locus allele frequencies without effect sizes from an association study. For a single SNP, let \textbf{p} be the vector of allele frequencies across populations. We then make a multivariate normal approximation to obtain a distribution with which we can model these frequencies (Coop...):

\begin{equation}
\textbf{p} \sim MVN( e, e(1-e)\textbf{F} )
\end{equation}
