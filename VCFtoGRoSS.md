# Converting a VCF file to a GRoSS input file
By Gabriel Renaud

The goal of this section is to show how to convert a multi-individual VCF file into GRoSS input data using [glactools](https://grenaud.github.io/glactools/).  We will replicate a signal of selection on the lactase genes in Europeans
using the 1000 genomes data. To simplify the analysis, we will use the super population panels: Africans (AFR), Europeans (EUR), East Asians (EAS), South Asians (SAS) and Amerindian (AMR).

First, we need to inform glactools about the names and length of the different chromosomes in the genome reference.  We will download the fasta index from the 1000 genomes data:

     wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai

In order to merge individuals from the same continent we will use the following panel information:

     wget https://personal.broadinstitute.org/armartin/ginger/integrated_call_samples_v3.20130502.ALL.panel.txt
     grep -v ^sample integrated_call_samples_v3.20130502.ALL.panel.txt  |cut -f 1,3  > panel.txt

The panel file contains a column that details the ID of the individual
and a second column with the population of origin:

     head -n 5 panel.txt
     HG00096 EUR
     HG00097 EUR
     HG00099 EUR
     HG00100 EUR
     HG00101 EUR

in our case:

    cut -f 2 panel.txt |sort | uniq
    AFR
    AMR
    EAS
    EUR
    SAS
(AFR=African,  AMR=Amerindian, EAS=East Asian,EUR=Europeans, SAS=South Asian).

We will then convert the VCF file into input for GRoSS:

     tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 2:136486829-136653337  |glactools vcfm2acf --onlyGT --fai human_g1k_v37.fasta.fai  -  |glactools meld -f panel.txt -  | glactools acf2gross --noroot -  |gzip > input.gross.gz

The first command uses tabix to retrieve a specific region with the
lactase gene.  The second command converts the VCF file into ACF.  The
third command merges the different individuals from different
super-populations together as to get an allele count per population.
Finally, the fourth command takes this ACF file and converts into
input for GRoSS.

In order to run gross,  we need a population tree:

    echo -e "digraph G {
    R -> AFR [ ] ;
    R -> EURASIA [ ] ;
    EURASIA -> EUR [ ] ;
    EURASIA -> ASIA [ ] ;
    ASIA -> AMR [ ] ;
    ASIA -> EASIA [ ] ;
    EASIA -> SAS [ ] ;
    EASIA -> EAS [ ] ;
    } " > treeGross.dot;


We simply need to execute GRoSS on the allele counts and the tree we
have previously created:

     prevdir=`pwd`
     cd ~/path to gross/GRoSS/
     Rscript ~/path to gross/GRoSS/GRoSS.R -d $prevdir/tree.dot -e
$prevdir/input.gross.gz -o /dev/stdout |grep -v "^\[" |bgzip -c >
$prevdir/output.gross.gz
     cd -

The file:

     output.gross.gz

...contains the output from GRoSS. We will plot the p-value for the
European branch:

     #!/usr/bin/env Rscript

    data <- read.table("EUR.gz",header=TRUE,stringsAsFactors=FALSE);

    pdf("pvalEUR.pdf")
    plot(data$START,-1*log(data$Pval_EUR_EURASIA),xlab="position",ylab="-log(P-value EUR)",main="P-value EUR per position",pch=18,col="darkblue"); dev.off();

The following shows the peak of p-values around the locus where
selection occurred.
