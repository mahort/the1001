# TODO: Add comment
# 
# Author: matt.horton
###############################################################################
rm(list=ls());

buildFSTjob <- function(
		baseDirectory,
		callMethod,
		contrastPopulation1,
		contrastPopulation2,
		snpImputationThreshold,
		chunkSize,
		ecotypeMap,
		prefix,
		jobName){

	scriptsDirectory = paste( baseDirectory, "/qsub_scripts/", sep="");
	if( !file.exists(scriptsDirectory)){
		dir.create(scriptsDirectory);
	}

	if( is.null( contrastPopulation1 )){
		filename <- paste(scriptsDirectory, "/qs.global", prefix, ".fst.sh", sep="");

	} else {
		filename <- paste(scriptsDirectory, "/qs.", contrastPopulation1, "_vs_", contrastPopulation2, ".fst.sh", sep="");	
	}

	qsubString <- paste( 
			"#!/bin/bash\n",
			"#PBS -S /bin/bash\n",
			"#PBS -N ", jobName, "\n",
			"#PBS -P the1001genomes\n",
			"#PBS -V\n",
			"#PBS -l walltime=12:00:00\n",
			"#PBS -l select=1:mem=", ifelse( is.null(contrastPopulation1), "64gb", "64gb"), "\n",
			"#PBS -J 1-5\n",
			"source /home/GMI/matt.horton/module_xcms.sh\n", sep="");

	qsubString <- paste( qsubString, 
			"cd ", baseDirectory, "/results/\n", sep="");

	qsubString <- paste(
			qsubString,
			"R --slave \"--args",
			" callMethodId=", callMethod, 
			" nCutoff=10", sep="");

	if( !is.null(contrastPopulation1) & !is.null(contrastPopulation2 )){
		qsubString <- paste(
				qsubString, 
				" pop1='", contrastPopulation1, "'",
				" pop2='", contrastPopulation2, "'",
				sep="");
	}

	qsubString <- paste(qsubString, 
					ifelse( is.null(prefix), " prefix=NULL", paste(" prefix='", prefix, "'", sep="")),
					" chr=${PBS_ARRAY_INDEX} threshold=", snpImputationThreshold, " batchSize=", chunkSize,
					ifelse( is.null(ecotypeMap), 
					" map=NULL", 
					paste( " map='", ecotypeMap, "' ", sep="")), "\" < ~/projects/code/svc_hdf5_calcWeirsTheta.R >> ${PBS_JOBNAME}.${PBS_JOBID}.Rout", sep="");

	write.table(qsubString, file=filename, quote=F, row.names=F, col.names=F);
}

#####################################################################################################################
## HARD CODED VARIABLES.
#####################################################################################################################
callMethodId <- 2;
snpImputationThreshold <- 0;
hdf5chunkSize <- 5e3;

admixedFile <- TRUE;

dataset <- "pbsrelicts";
dataset <- match.arg(dataset, c("relict", "pbsrelicts", "global", "pairwise"));

################################################################################
### use the accession map...
################################################################################
if( dataset == "relict" ){
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/pairwise/admixed/", sep=""));
	mapFileName <- "~/projects/the_1001/selection/fst/2/imp1/global/results/relicts_nonrelict.only.txt";
	accessionMap <- read.table(mapFileName, header=T, sep="\t", as.is=T, stringsAsFactors=F);

} else if( dataset == "pbsrelicts" ){
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/pairwise/admixed/", sep=""));
	mapFileName <- "~/projects/the_1001/selection/fst/2/imp1/global/results/pbs_relicts.txt";
	accessionMap <- read.table(mapFileName, header=T, sep="\t", as.is=T, stringsAsFactors=F);

} else if( admixedFile & ( dataset != "relict")){
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/pairwise/admixed/", sep=""));
	mapFileName <- "~/projects/the_1001/selection/fst/admix_pop_metadata/metadata_structure.txt";
	accessionMap <- read.table(mapFileName, header=T, sep="\t", as.is=T, stringsAsFactors=F);
	colnames(accessionMap)[which(colnames(accessionMap) == "tg_ecotypeid")] <- "ecotype_id";
	colnames(accessionMap)[which(colnames(accessionMap) == "Admixed")] <- "population_name";
	mapFileName <- NULL;			

} else {
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/pairwise/popn/", sep=""));
	mapFileName <- "~/projects/the_accessions/the1001.genomes.txt";
	accessionMap <- read.table(mapFileName, header=T, sep="\t", as.is=T, stringsAsFactors=F, quote="");	
}



if( dataset == "global" ){
	## we should count the samples, and use a compatible memory requirement
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/global/", sep=""));

	buildFSTjob(
			baseDirectory = getwd(),
			callMethod = callMethodId,
			contrastPopulation1 = NULL,
			contrastPopulation2 = NULL,
			snpImputationThreshold = snpImputationThreshold,
			chunkSize = hdf5chunkSize, 
			ecotypeMap = mapFileName,
			prefix = ifelse(admixedFile, "adm_", "popn_"),
			jobName = paste("gl.fst", sep=""));

} else if( dataset == "pairwise" | dataset == "relict" | dataset == "pbsrelicts" ){
	cat("preparing:", dataset, '\n');
	subpops <- unique(accessionMap[,"population_name"]);
	uniquePopulations <- sort(subpops[which(!is.na(subpops))]);
	uniquePopulations <- gsub(" ", "_", uniquePopulations);
	
	for( j in 1:(length(uniquePopulations) - 1)){
		pop1 <- uniquePopulations[j];
		
		for( k in (j+1):length(uniquePopulations)){
			pop2 <- uniquePopulations[k];	
			
			## we should count the samples, and use a compatible memory requirement
			buildFSTjob(
					baseDirectory = getwd(),
					callMethod = callMethodId,
					contrastPopulation1 = pop1,
					contrastPopulation2 = pop2,
					snpImputationThreshold = snpImputationThreshold,
					chunkSize = hdf5chunkSize, 
					ecotypeMap = mapFileName,
					prefix = NULL,
					jobName = paste("fst2", sep=""));
		}
	}
}
