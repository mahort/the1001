# TODO: Add comment
# 
# Author: matt.horton
###############################################################################
rm(list=ls());

buildQsubFilesForTheLoop <- function(
		baseDirectory,
		resultsDirectory,
		categoryColumn,
		cutoff,
		numberOfCores,
		windowedGwasResults,
		maf,
		windowSize,
		phenotypeSuffix,
		jobName){

	scriptsDirectory = paste( baseDirectory, "/qsub_scripts/", sep="");
	if( !file.exists(scriptsDirectory)){
		system(paste("mkdir ", scriptsDirectory, sep=""));
	}

	filename <- paste(scriptsDirectory, "/qs.loop.", categoryColumn, ".t", cutoff, ".win", windowSize, "_",  phenotypeSuffix, ".", gsub(".txt", ".sh", basename(windowedGwasResults)), sep="");

	qsubString <- paste( 
			"#!/bin/bash\n",
			"#PBS -S /bin/bash\n",
			"#PBS -N ", jobName, "\n",
			"#PBS -P cegs\n",
			"#PBS -V\n",
			"#PBS -l walltime=00:10:00\n",
			"#PBS -l select=1:ncpus=", numberOfCores, ":mem=2GB\n",
			"source /home/GMI/matt.horton/module_intel.sh\n", sep="");
		
	qsubString <- paste( qsubString, 
			"cd ", resultsDirectory, "\n", sep="");

	qsubString <- paste(
			qsubString,
			"R --slave \"--args",
			" categoryColumn='",categoryColumn,
			"' cutoff=",cutoff,
			" numberOfCores=${OMP_NUM_THREADS}",
			" windowedGwasResults='", windowedGwasResults, "'",
			" mafThreshold=", maf, 
			" windowSize=", windowSize, "\" < ~/projects/code/2b_loop.enrichment.intail.fet.win.turbo.R >> ${PBS_JOBNAME}.${PBS_JOBID}.Rout", sep="");

	write.table(qsubString, file=filename, quote=F, row.names=F, col.names=F);
}

#####################################################################################################################
## HARD CODED VARIABLES.
#####################################################################################################################
# the snp file determines the maf cutoff. 
# prune table to unique set.
# 1.) go column
goColumn <- "goterm";

# 2.) the cutoffs
tailCutoff <- 0.05;

# 3.) number of cores
nCores <- 8;

# 5.) you can extend the gene model to search for overlap... not sure if this is important, so test w/ a windowSize of 0, 2500, 5000, etc. (?).
resultsWindowedSize <- 1e4;
geneWindowSize <- 1000;
snpImputationThreshold <- 0.1;

# 6.) test the various snp call methods.
callMethodId <- 2;
resultsPath <- paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/", sep="");
setwd(resultsPath);


methodMatch <- "relict|global|_vs_";

files <- list.files(path=getwd(), pattern=methodMatch, recursive=T);
files <- files[grep(paste(".win", resultsWindowedSize, ".scores.txt", sep=""), files)];
cat("Found:", length(files), "matching files.\n");
print(files);

for( i in 1:length(files)){
	buildQsubFilesForTheLoop(
			baseDirectory = resultsPath,
			resultsDirectory = getwd(),
			categoryColumn = goColumn, 
			jobName = paste("gsea_", "fst", sep=""),
			cutoff = tailCutoff,
			numberOfCores = nCores,
			windowedGwasResults = files[i],
			maf = 0,  
			windowSize = geneWindowSize,
			phenotypeSuffix = "");
}
