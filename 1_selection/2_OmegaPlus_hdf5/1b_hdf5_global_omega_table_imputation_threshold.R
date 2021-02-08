# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());
require(rhdf5);
numberOfReads <- 1e4;

snpImputationThreshold <- 0.1;
callMethodId <- 2; ## gots to put 'em somewhere.

###############################################################################
## setup the hdf5 connection and retrieve metadata.
###############################################################################
hdf5Path <- paste("/lustre/scratch/projects/the1001genomes/VCF_1135g/", sep="");
filename <- paste(hdf5Path, "1135g_SNP_BIALLELIC.SNPmatrix_6-Oct-2015_uncompressed.hdf5", sep="");

accessions <- as.numeric(h5read(filename, "accessions"));

#attributes <- h5ls(filename, all=T);
positions <- h5read(filename, name="positions", read.attributes=T);
chromosomes <- attributes(positions)$chr_regions;
chromosomes[1,] <- chromosomes[1,] + 1;
chr.ids <- attributes(positions)$chrs;
snpsPerChromosome <- apply(chromosomes, 2, function(x){ x[2] - x[1]; }) + 1;

###############################################################################
## omega+ parameters
###############################################################################
gridInterval <- 500; ## grid size, not selection scan size.
selectionScanMinWindowStart <- 2e3;
selectionScanMaxWindowStop <- 5e4;
range <- paste(selectionScanMinWindowStart, selectionScanMaxWindowStop, sep="-");

###############################################################################
## code to setup pbs jobs - 3 functions, and (below) the code to setup
## each (population specific) job.
###############################################################################
buildOmegaPlusJob <- function(
		callMethodId,	
		jobName,
		sizeOfSample,
		chromosomeNumber,
		numberOfSlices,
		lengthOfSegment,
		gridIntervalDensity,
		minimumSweepWindowSize,
		maximumSweepWindowSize,
		snpImputationThreshold){
	
	filename <- paste(callMethodId, "/imp", snpImputationThreshold, "/omega_", jobName, ".", sizeOfSample, ".chr", chromosomeNumber, ".grd", gridIntervalDensity, ".win", range, ".om.sh", sep="");
	
	qsubString <- paste( 
			"#!/bin/bash\n",
			"#PBS -S /bin/bash\n",
			"#PBS -N op", chromosomeNumber, "\n",
			"#PBS -P the1001genomes\n",
			"#PBS -V\n",
			"#PBS -l walltime=48:00:00\n",
			"#PBS -l select=1:ncpus=1:mem=16gb\n", sep="");
	
	if( numberOfSlices > 1 ){
		qsubString <- paste(qsubString, "#PBS -J 1-", numberOfSlices, "\n", sep="");
	}
	
	qsubString <- paste( qsubString, 
			"SUBPOP=", jobName, "\n",
			"SIZE=", sizeOfSample, "\n",
			"CHR=", chromosomeNumber, "\n", sep="");
	if( numberOfSlices > 1 ){
		qsubString <- paste(qsubString, "JOB=${PBS_ARRAY_INDEX}\n", sep="");
		
	} else {
		qsubString <- paste(qsubString, "JOB=1\n", sep="");
	}
	
	qsubString <- paste(qsubString, 
			"cd ~/projects/the_1001/selection/omega/", callMethodId, "/imp", snpImputationThreshold, "/${SUBPOP}/chr${CHR}/slice${JOB}\n",
			"INPUTFILE=${SUBPOP}_chr${CHR}.n_${SIZE}.slice${JOB}.omega\n", sep="");
	
	qsubString <- paste( qsubString, "OmegaPlus -name ${INPUTFILE}.grd", gridIntervalDensity, ".win", minimumSweepWindowSize, "-", maximumSweepWindowSize, " -input ${INPUTFILE}.txt -grid ", floor(lengthOfSegment/gridIntervalDensity), " -minwin ", minimumSweepWindowSize, " -maxwin ", maximumSweepWindowSize, " -length ", lengthOfSegment, " -all\n", sep="");
	
	write.table(qsubString, file=filename, quote=F, row.names=F, col.names=F);
}

###############################################################################
## distinguish between the lines we want, and those we don't...
###############################################################################
metadata <- "~/projects/the_1001/selection/fst/admix_pop_metadata/metadata_structure.txt";
metadata <- read.table(metadata, header=T, sep="\t", as.is=T, stringsAsFactors=F);
colnames(metadata)[which(colnames(metadata) == "Admixed")] <- "population_name";
colnames(metadata)[which(colnames(metadata) == "tg_ecotypeid")] <- "ecotype_id";
#subpops <- unique(metadata[,"population_name"]);
#metadata <- metadata[which(!is.na(metadata[,"population_name"])),]; # %in% subpops),];

setwd("~/projects/the_1001/selection/omega");

hdf5EcotypeIdIndices <- which(accessions %in% metadata[,"ecotype_id"]);
cat("We were able to recover the indices for the following lines:\n", accessions[hdf5EcotypeIdIndices], "\n");
cat("Total # of accessions:", length(hdf5EcotypeIdIndices), "\n");

for( j in 1:5 ){

	cat("Working on chromosome:", j, "\n");
	indices.chr_j <- chromosomes[1,j]:chromosomes[2,j];
	positions.chr_j <- positions[indices.chr_j];
	totalLength <- positions.chr_j[length(positions.chr_j)] - positions.chr_j[1];

	chr_j <- h5read(filename, name="snps", index=list(hdf5EcotypeIdIndices, indices.chr_j));
	cat("There are:", ncol(chr_j), "SNPs from sequencing/genotyping this population.\n");

	## handle monomorphic SNPs ourselves.
	proportionOfNs <- apply(chr_j, 2, function(x){ sum( x == -1)})/nrow(chr_j);
	dropouts <- which( proportionOfNs > snpImputationThreshold);
	if( length(dropouts) > 0 ){
		cat("There were:", length(dropouts), "dropouts.\n");
		chr_j <- chr_j[,-dropouts];
		positions.chr_j <- positions.chr_j[-dropouts];
	}

	sums <- colSums(chr_j);
	dropouts <- which(sums == 0 | sums == nrow(chr_j));
	if( length(dropouts) > 0 ){
		cat("There were:", length(dropouts), "dropouts.\n");
		chr_j <- chr_j[,-dropouts];
		positions.chr_j <- positions.chr_j[-dropouts];
	}

	chr_j[chr_j == -1] <- "N";

	cat(ncol(chr_j), "of these are biallelic.\n");	
	outputDirectory <- paste(callMethodId, "/imp", snpImputationThreshold, "/global/chr", j, "/slice1", sep="");
	
	if( !file.exists(outputDirectory)){
		cat("Making directory(ies) now.\n");
		dir.create(outputDirectory, recursive=T);
		
		if( !file.exists(outputDirectory)){
			cat("There was an issue creating the directory.\nRetrying now.\n");
			Sys.sleep(10);
			dir.create(outputDirectory, recursive=T);
		}
	}
	
	outputFileName <- paste(callMethodId, "/imp", snpImputationThreshold,  "/global/", "chr", j, "/slice1/global_chr", j, ".n_", nrow(chr_j), ".slice1.omega.txt", sep="");
	
	sink(outputFileName)
	cat("ms", nrow(chr_j), "1 -s", ncol(chr_j), "\n", sep=" ");
	cat(positions.chr_j[1], totalLength, gridInterval, "\n", sep=" ");
	cat("\n//\n");
	cat("segsites:", ncol(chr_j), "\n", sep=" ");
	
	positions.chr_j <- positions.chr_j/totalLength;
	cat("positions:", positions.chr_j, "\n", sep=" ");
	
	for( k in 1:nrow(chr_j)){
		cat(paste( chr_j[k,], collapse="", sep=""), "\n");
	}
	
	sink(NULL);
	
	buildOmegaPlusJob( jobName="global", 
			sizeOfSample=nrow(chr_j), 
			chromosomeNumber=j, 
			numberOfSlices=(1), 
			callMethodId=callMethodId, 
			lengthOfSegment=totalLength, 
			gridIntervalDensity=gridInterval, 
			minimumSweepWindowSize=selectionScanMinWindowStart, 
			maximumSweepWindowSize=selectionScanMaxWindowStop,
			snpImputationThreshold=snpImputationThreshold);
	
	cat("Files for chr:", j, "written.\n\n");
	gc();
}

