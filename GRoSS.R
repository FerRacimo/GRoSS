#!/usr/bin/env Rscript
library("optparse")

source("MultiBranchFunc.R")


# Options
option_list = list(
  make_option(c("-e", "--inputfile"), type="character", default=NULL, help="Name of input file with SNP allele counts."),
  make_option(c("-r", "--graphfile"), type="character", default=NULL, help="Graph R file name in qpGraph format; don't use this if using the -d option"),
  make_option(c("-d", "--dotfile"), type="character", default=NULL, help="Dot file name; don't use this if using the -r option."),
  make_option(c("-o", "--outfile"), type="character", default=NULL, help="Output file name."),
  make_option(c("-s", "--finitesamp"), action = "store_true", default = FALSE,help = "Run GRoSS accounting for finite sample sizes. [default=FALSE]")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$inputfile)){
  print_help(opt_parser)
  stop("Input file name must be supplied.n", call.=FALSE)
}
if (is.null(opt$outfile)){
  print_help(opt_parser)
  stop("Output file name must be supplied.n", call.=FALSE)
}
if (is.null(opt$graphfile) & is.null(opt$dotfile)){
  print_help(opt_parser)
  stop("Graph file name must be supplied.n", call.=FALSE)
}


neutfile <- opt$inputfile
graphfile <- opt$graphfile
dotfile <- opt$dotfile
outfile <- opt$outfile
finitesamp <- opt$finitesamp

# Load graph file
source("LoadFiles.R")
