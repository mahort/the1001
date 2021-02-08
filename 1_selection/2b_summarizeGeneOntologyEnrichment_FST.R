# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

# gsea (my window version)
#1.) window the gwas results
#2.) determine the significance of the results 
#3.) iterate over the results, w/ FDR in mind
#a.) open the file
#b.) truncate to in_tail_in_gene_set > 1 & in_gene_set > 2 (maybe make the second one optional)
#c.) calc FDR e.g. at 0.1, and make tables of the results (using a minimum # of records to report)
#

#####################################################################################################################
## HARD CODED VARIABLES.
#####################################################################################################################
# the snp file determines the mac cutoff. 
# prune table to unique set.
rm(list=ls());
require(qvalue);

goColumn <- "goterm";
#tailCutoff <- 0.05;

#####################################################################################################################
## pick the panel
#####################################################################################################################
callMethodId <- 2;
snpImputationThreshold <- 0.1;

#####################################################################################################################
## move there.
#####################################################################################################################
setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/", sep=""));

#####################################################################################################################
## PREP FILES.
#####################################################################################################################
# match files of the format:
# win.scores.win_1000.all_categories.t0.05.goterm.fet
methodMatch <- "fet$";
suffix <- "";

files <- list.files(path=getwd(), pattern=methodMatch, recursive=T);
print(files);

minimumNumberOfRecords <- 10;
minimumNumberOfGenesInTail <- 3;
fdrThreshold <- 0.1;

for( filename_i in files ){
	cat("Analyzing file:", filename_i, "\n");

	file_i <- read.table(filename_i, header=T, sep="\t", as.is=T, stringsAsFactors=FALSE, quote="");
	interesting <- subset(file_i, in_tail_in_gene_set  >= minimumNumberOfGenesInTail );

	cat("Number of categories:", nrow(interesting), "\n");

	# now make a table of these results...
	# sort, and calc FDR correction @ 0.1 (just to be conservative).
	# interesting <- interesting[order(interesting[,"fisher_pvalue"]),];
	qvals <- qvalue(interesting[,"fisher_pvalue"], lambda=0.5, fdr=fdrThreshold);
	interesting <- cbind( interesting, q_value=qvals$qvalues );
	significant <- which( qvals$significant );
	reportThreshold <- max(minimumNumberOfRecords, length(significant));
	interesting <- interesting[order(interesting[,"q_value"], interesting[,"fisher_pvalue"]),];
	
	# I want this to fail if nrow(interesting) somehow lower than the reportThreshold so that I can further investigate. i.e. not catching it right now.
	interesting <- interesting[1:reportThreshold,];
	write.table(interesting, paste(filename_i, ".gte", minimumNumberOfGenesInTail, ".fdr", fdrThreshold, sep=""), quote=FALSE, row.names=FALSE, sep="\t");
}