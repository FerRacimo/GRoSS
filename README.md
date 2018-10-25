# GRoSS
Graph-aware Retrieval of Selective Sweeps

Here, I introduce a method to detect selective sweeps across the genome, when using many populations that are related via a complex admixture graph. I made some slight modifications to the Q<sub>B</sub> statistic from Racimo, Berg and Pickrell (2018) which was originally meant to detect polygenic adaptation using admixture graphs (see https://github.com/FerRacimo/PolyGraph). The new statistic - which I call S<sub>B</sub> - does not need GWAS data and works with allele frequency data alone. It can be used to **both scan the genome for regions under strong single-locus positive selection, and pinpoint where in the graph the selective event most likely took place.** See the file detecting-single-locus.pdf for an explanation of how the statistic works.

If you end up using it, please cite the following preprint:
https://www.biorxiv.org/content/early/2018/10/25/453092

# Required R Libraries

```
R
install.packages("msm")
install.packages("reshape2")
install.packages("pscl")
install.packages("parallel")
install.packages("data.table")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("admixturegraph")
install.packages("devtools")
install.packages("readr")
install.packages("qqman")
install.packages("optparse")
devtools::install_bitbucket('coolbutuseless/minilexer')
devtools::install_github("mailund/graphparse")
devtools::install_github("mailund/matchbox")
source("https://bioconductor.org/biocLite.R"); biocLite("biomaRt")
```

# Examples

As a visual example of what GRoSS can do, here are two graphs showing -log10(p-values) for the S<sub>B</sub> statistic computed on SNPs along the genome. The first graph was made using populations from the 1000 Genomes Project Phase 3. The second graph was made using populations from Lazaridis et al. (2014), after imputation on the 1000 Genomes data.

<img src="https://github.com/FerRacimo/GRoSS/blob/master/Q_b_manhattan_human.png" height="800">


# Running GRoSS

Here is an example line for generating the above results. The main R script is GRoSS.R and it requires the user to specify:
- an input file (\*txt) specified with the -e option
- a graph file describing the topology of the graph. This can be:
  * in the same format as the graph file that is used as input for qpGraph (Patterson et al. 2012) (note that the value of the fitted admixture weights must be correctly stated), specified with the -r option OR
  * in dotfile format (outputted from qpGraph after fitting), specified with the -d option 
- an output file on which to write the results, specified with the -o option

We will use the same 1000 Genomes populations as in the example above, but limiting ourselves to chr22. First, unpack the input file KG_popfile_chr22.txt.gz.

```
gzip -c KG_popfile_chr22.txt.gz > KG_popfile_chr22.txt
```

Then, run the program:

```
Rscript GRoSS.R -e KG_popfile_chr22.txt -r 1KG_MSL_ESN_CDX_JPT_CEU_TSI_CHB.graph -o SNPstat_1KG_MSL_ESN_CDX_JPT_CEU_TSI_CHB_chr22.tsv
```

The output contains the chi-squared statistics and the corresponding p-values for each branch of the admixture graph.

The Plotting_and_Windowing.txt file contains some R commands to merge SNPs into windows and plot the output.
