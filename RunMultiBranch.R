#!/usr/bin/env Rscript
library("optparse")

source("MultiBranchFunc.R")


# Options
option_list = list(
  make_option(c("-e", "--neutfile"), type="character", default=NULL, help="Neutral input file name"),
  make_option(c("-r", "--graphfile"), type="character", default=NULL, help="Graph R file name"),
  make_option(c("-o", "--outfile"), type="character", default=NULL, help="Output file")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$neutfile)){
  print_help(opt_parser)
  stop("Neutral file name must be supplied.n", call.=FALSE)
}
if (is.null(opt$graphfile)){
  print_help(opt_parser)
  stop("Graph file name must be supplied.n", call.=FALSE)
}
if (is.null(opt$outfile)){
  print_help(opt_parser)
  stop("Output file name must be supplied.n", call.=FALSE)
}

neutfile <- opt$neutfile
graphfile <- opt$graphfile
outfile <- opt$outfile

# Default (before going through graphfile)
pvaltotal <- 0.05 

# Load graph file
source(graphfile)


