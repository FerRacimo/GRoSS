# GRoSS
Graph-aware Recovery of Selective Sweeps

We modified the $Q_{B}$ statistics (Racimo et al. 2018) to handle single-locus allele frequencies without effect sizes from an association study. For a single SNP, let \textbf{p} be the vector of allele frequencies across populations. We then make a multivariate normal approximation to obtain a distribution with which we can model these frequencies (Coop...):

\begin{equation}
\textbf{p} \sim MVN( e, e(1-e)\textbf{F} )
\end{equation}
