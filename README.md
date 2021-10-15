# GRoSS
Graph-aware Retrieval of Selective Sweeps

Here, I introduce a method to detect selective sweeps across the genome, when using many populations that are related via a complex admixture graph. I made some slight modifications to the Q<sub>B</sub> statistic from Racimo, Berg and Pickrell (2018) which was originally meant to detect polygenic adaptation using admixture graphs (see https://github.com/FerRacimo/PolyGraph). The new statistic - which I call S<sub>B</sub> - does not need GWAS data and works with allele frequency data alone. It can be used to **both scan the genome for regions under strong single-locus positive selection, and pinpoint where in the graph the selective event most likely took place.** See the file detecting-single-locus.pdf for an explanation of how the statistic works.

If you end up using it, please cite the following paper: https://genome.cshlp.org/content/29/9/1506

Free preprint version: https://www.biorxiv.org/content/10.1101/453092v2



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
install.packages("devtools")
install.packages("readr")
install.packages("qqman")
install.packages("optparse")
install.packages("admixturegraph",repos=unique(c(getOption("repos"),repos="https://cran.microsoft.com/snapshot/2019-04-01/")))
devtools::install_github("mailund/graphparse")
source("https://bioconductor.org/biocLite.R"); biocLite("biomaRt")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("biomaRt"))
```

# Examples

As a visual example of what GRoSS can do, here are two graphs showing -log10(p-values) for the S<sub>B</sub> statistic computed on SNPs along the genome. The first graph was made using populations from the 1000 Genomes Project Phase 3. The second graph was made using populations from Lazaridis et al. (2014), after imputation on the 1000 Genomes data.

<img src="https://github.com/FerRacimo/GRoSS/blob/master/Q_b_manhattan_human.png" height="800">


# Running GRoSS

Here is an example line for generating the above results. The main R script is GRoSS.R and it requires the user to specify:
- an input file (\*txt) specified with the -e option
  * The first column of this file is the chromosome name, the second column is the position, the third column is a SNP ID identifier. The identifier can just be equal to "[chromosome]-[position]" if no SNP id is available, but it must be unique to each SNP. All the other columns contain the number of reference and alternative alleles in each population, separated by a comma (e.g. 5,8 if there are 5 reference alleles and 8 alternative alleles in a given population at a given SNP).
  * NOTE: GRoSS will ignore any line where one or more population panels have missing data ("0,0").
- a graph file describing the topology of the graph. This can be:
  * in the same format as the graph file that is used as input for qpGraph (Patterson et al. 2012) (note that the value of the fitted admixture weights must be correctly stated), specified with the -r option OR
  * in dotfile format (outputted from qpGraph after fitting), specified with the -d option
  * NOTE: do not add comments (lines starting with "#") to either of these files
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

The Plotting_and_Windowing.txt file contains some R commands to optionally merge SNPs into windows and plot the output.

For an example of how to create a graph file with admixture, see: https://github.com/DReichLab/AdmixTools/blob/master/examples.qpGraph/gr1x


# A guide for converting a VCF file into a GRoSS input file

If you have a VCF file and would like to convert it to a GRoSS input file (e.g. as in KG_popfile_chr22.txt), then you can use this handy guide created by Gabriel Renaud: https://github.com/FerRacimo/GRoSS/blob/master/VCFtoGRoSS.md

# New feature: correction for low sample sizes

We've recently implemented a modified of the Q_S statistic for cases of low sample sizes. A full description of it can be found here[https://github.com/FerRacimo/GRoSS/blob/1c28b48adca72396c58c70e9084090772ecfa36f/GRoSS_low_sample_size.pdf]
The modified statistic uses a Normal approximation to the binomial distribution that accounts for the increased variance in sample allele frequencies, relative to population allele frequencies, as a consequence of finite sample sizes. We recommend the use of this feature over the "vanilla" version of GRoSS, especially when the number of (diploid) individuals in at least one of the populations is less than 20. It can be run by using the "-s" option:

```
Rscript GRoSS.R -e KG_popfile_chr22.txt -r 1KG_MSL_ESN_CDX_JPT_CEU_TSI_CHB.graph -o SNPstat_1KG_MSL_ESN_CDX_JPT_CEU_TSI_CHB_chr22.tsv -s
```
