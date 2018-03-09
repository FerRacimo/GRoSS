# GRoSS
Graph-aware Recovery of Selective Sweeps

Here, I introduce a method to detect selective sweeps across the genome, when using many populations that are each related via a complex admixture graph. I made some slight modifications to the Q<sub>B</sub> statistic from Racimo, Berg and Pickrell (2018) which was originally meant to detect polygenic adaptation using admixture graphs. The new statistic - which I call R<sub>B</sub> - does not need GWAS data and works with allele frequency data alone. It can be used to **both to scan the genome for regions under strong single-locus positive selection, and to pinpoint where in the graph the selective event most likely took place.** See the file GRoSS.pdf for an explanation of how the statistic works.

Eventually I'll publish a paper on this statistic, but in case you want to use it and not wait for that to come out, you can just cite the paper below and add a link to this github page.

Racimo, F., Berg, J. J., & Pickrell, J. K. (2018). Detecting polygenic adaptation in admixture graphs. Genetics, genetics-300489.

Special acknowledgments to Jeremy Berg who was the original mastermind behind the Q<sub>B</sub> statistic.

**NOTE**: The R<sub>B</sub> statistic is fast and easy to compute but is not as principled as other approaches for multi-population selection. For once, it doesn't rely on an explicit positive selection model, it just detects significant deviations from neutrality. For a more principled approach, see, for example: Lee, K. M., & Coop, G. (2017). Distinguishing among modes of convergent adaptation using population genomic data. Genetics, 207(4), 1591-1619.

# Examples

As a visual example of what GRoSS can do, here are two graphs and their corresponding R<sub>B</sub> statistics computed across the genome. The first graph was made using populations from the 1000 Genomes Project Phase 3. The second graph was made using populations from Lazaridis et al. (2014), after imputation on the 1000 Genomes data.

<img src="https://github.com/FerRacimo/GRoSS/blob/master/Q_b_manhattan_1000G.png" height="300">

<img src="https://github.com/FerRacimo/GRoSS/blob/master/Q_b_manhattan_LazCombo.png" height="300">


# Running GRoSS

Here are two example lines for generating the above results. The main R script is RunMultiBranch.R and it requires the user to specify:
- an input file (\*txt) specified with the -e option
- another R script that specifies the graph using the admixturegraph package (Leppälä et al. 2017), specified with the -r option
- an output file on which to write the results, specified with the -o option

1000 Genomes scan example:

Rscript RunMultiBranch.R -e KG_popfile.txt -r 1KG_MSL_ESN_CDX_JPT_CEU_TSI_CHB.R -o SNPstat_1KG_MSL_ESN_CDX_JPT_CEU_TSI_CHB.tsv

Lazaridis et al. scan example:

Rscript RunMultiBranch.R -e LazCombo_popfile.txt -r LAZ_Mandenka_Yoruba_Oceanian_EastAsian_NativeAmerican_Sardinian_EuropeA.R -o SNPstat_LAZ_Mandenka_Yoruba_Oceanian_EastAsian_NativeAmerican_Sardinian_EuropeA.tsv

The output contains the chi-squared statistics and the corresponding p-values for each branch of the admixture graph.

The Plotting_and_Windowing.txt file contains some R commands to merge SNPs into windows and plot the output.
