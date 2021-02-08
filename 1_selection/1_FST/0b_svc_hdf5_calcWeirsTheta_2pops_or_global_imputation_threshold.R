# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());
require(rhdf5);

###############################################################################
# core function: Weir & Cockerham (1984) theta. this is FST
# this code expects homozygous AA to be coded with the number: 2. Heterozygote Aa: 1. Homozygote aa: 0.
###############################################################################
Fst <- function(n, ns, ps, hs){ 
	
	nbar <- mean(ns);
	pbar <- sum(ns*ps)/(nbar*n);
	
	nc <- (n*nbar - sum(ns^2)/(n*nbar))/(n-1);
	
	s2 <- sum(ns*((ps - pbar)^2))/((n-1)*nbar);
	hbar <- sum(ns*hs)/(n*nbar);
	
	# from pages 1359, 1360 in Weir & Cockerham (1984)
	a <- (nbar/nc)*(s2 - (1/(nbar-1))*(pbar*(1-pbar) - ((n-1)/n)*s2 - (1/4)*hbar));
	b <- (nbar/(nbar-1))*((pbar*(1-pbar)) - ((n-1)/n)*s2 - ((2*nbar-1)/(4*nbar))*hbar);
	c <- (1/2)*hbar;

	return(a/(a+b+c));
}

#####################################################################################################################
## Handle incoming command line arguments.
#####################################################################################################################
args <- commandArgs(TRUE);
for( i in 1:length(args)){
	eval(parse(text=args[i])); #
}

stopifnot(
		exists("callMethodId"),
		exists("map"),
		exists("chr"),
		exists("threshold"),
		exists("batchSize"),
		exists("nCutoff"),
		exists("prefix"));

#
#################################################################################
#### metadata for the accessions, this is for testing...
#################################################################################
#map <- "~/projects/the_accessions/the1001.genomes.txt";
#threshold <- 0.1;
#batchSize <- 3000;
#chr <- 1;
#nCutoff <- 10;

#pop1 <- "Americas";
#pop2 <- "British-Isles";

chromosomeNumber <- chr;
snpImputationThreshold <- threshold;		## e.g. 0.1
numberOfReads <- batchSize;					## e.g. 1e4
sampleSizeCutoff <- nCutoff;				## e.g. 10
=

###############################################################################
## setup the hdf5 connection and retrieve metadata.
###############################################################################
hdf5Path <- paste("/lustre/scratch/projects/the1001genomes/VCF_1135g/", sep="");
filename <- paste(hdf5Path, "1135g_SNP_BIALLELIC.SNPmatrix_6-Oct-2015_uncompressed.hdf5", sep="");

accessions <- as.numeric(h5read(filename, "accessions"));
positions <- h5read(filename, name="positions", read.attributes=T);
chromosomes <- attributes(positions)$chr_regions;
#chromosomes[1,] <- chromosomes[1,] + 1;
chr.ids <- attributes(positions)$chrs;
snpsPerChromosome <- apply(chromosomes, 2, function(x){ x[2] - x[1]; }) + 1;


###############################################################################
## distinguish between the lines we want, and those we don't...
###############################################################################
if( is.null(map)){
	map <- "~/projects/the_1001/selection/fst/admix_pop_metadata/metadata_structure.txt";
	metadata <- read.table(map, header=T, sep="\t", as.is=T, stringsAsFactors=F);
	colnames(metadata)[which(colnames(metadata) == "Admixed")] <- "population_name";
	colnames(metadata)[which(colnames(metadata) == "tg_ecotypeid")] <- "ecotype_id";
	
} else {
	metadata <- read.table(map, header=T, sep="\t", as.is=T, stringsAsFactors=F);
}


if( exists("pop1") & exists("pop2")){
	usableNamePop1 <- pop1;
	if( grepl("_", pop1)){ pop1 <- gsub("_", " ", pop1); }

	usableNamePop2 <- pop2;
	if( grepl("_", pop2)){ pop2 <- gsub("_", " ", pop2); }

	uniquePopulations <- c(pop1, pop2);

} else {
	uniquePopulations <- table(metadata[,"population_name"]);
	uniquePopulations <- names(uniquePopulations[uniquePopulations >= sampleSizeCutoff]);
}


allIndices <- c();
listOfPopulations <- list();

for( j in 1:length(uniquePopulations)){
	pop_j <- uniquePopulations[j];
	
	subset_j <- subset(metadata, population_name == pop_j);
	targetIds <- subset_j[,"ecotype_id"];

	sampleSize <- length(targetIds);
	cat(pop_j, "has", sampleSize,"members.\n");
	indices <- which(accessions %in% targetIds);
	allIndices <- c(allIndices, indices);
	listOfPopulations[[j]] <- list("name"=pop_j, "ecotype_ids"=targetIds, "indices"=indices, "n"=sampleSize);
}

###############################################################################
# begin the analysis
###############################################################################
cat("This software assumes the data are in the form: \n");
cat("----------  2 for the homozygote (AA)\n");
cat("------ and: 1 for the heterozygote (Aa)\n"); ## these don't exist in our dataset.
cat("------ and: 0 for the homozygote (aa)\n");

populations <- list();
populations$n <- length(listOfPopulations);

cat("Retrieving snps for chromosome", chromosomeNumber, "\n");
indices.chr_j <- chromosomes[,chromosomeNumber];
breaks <- unique(c(seq(indices.chr_j[1], indices.chr_j[2], by=numberOfReads), chromosomes[2, chromosomeNumber])); #=ceiling(snpsPerChromosome[chromosomeNumber]/numberOfReads)));

results <- list();
for( k in 1:(length(breaks) - 1 )){
	cat("On slice:", k, "of", length(breaks), "\n");
	sliceStart <- breaks[k] + 1; 
	sliceStop <- breaks[k + 1];

	snps_k <- h5read(filename, name="snps", index=list(allIndices, sliceStart:sliceStop));
	rownames(snps_k) <- accessions[allIndices];
	pos_k <- positions[sliceStart:sliceStop];

	## from here, we need to eliminate the SNPs that don't pass the imputation Threshold
	proportionOfNs <- apply(snps_k, 2, function(x){ sum( x == -1)})/nrow(snps_k);
	dropouts <- which( proportionOfNs > snpImputationThreshold);
	if( length(dropouts) > 0 ){
		cat("There were:", length(dropouts), "dropouts.\n");
		snps_k <- snps_k[,-dropouts];
		pos_k <- pos_k[-dropouts];

		if( length(pos_k) == 1 ){
			dim(snps_k) <- c(length(snps_k), 1);
			rownames(snps_k) <- accessions[allIndices];

		} else if( length(pos_k) == 0 ){
			next;
		}
	}

	## that said, we still have SNPs flagged as -1 here (that is, N). We need to eliminate these before estimating theta for particular SNPs.
	is.na(snps_k[which(snps_k == -1)]) <- TRUE;
	snps_k <- snps_k * 2; ## and rescale.

	ps_k <- matrix(nrow=length(pos_k), ncol=populations$n);
	hs_k <- matrix(nrow=length(pos_k), ncol=populations$n);
	ns_k <- matrix(nrow=length(pos_k), ncol=populations$n);

	for( l in 1:length(listOfPopulations)){
		pop_l <- listOfPopulations[[l]];
		pop_l_snps <- snps_k[which(rownames(snps_k) %in% pop_l[['ecotype_ids']]), 1:ncol(snps_k)];
		if((is.null(dim(pop_l_snps))) & (length(pos_k) == 1)){
			dim(pop_l_snps) <- c(length(pop_l_snps), 1);
		}
		
		ps_k[,l] <- apply(pop_l_snps, 2, mean, na.rm=T) / 2;
		ns_k[,l] <- apply(pop_l_snps, 2, function(x){ length(which(!is.na(x))); });
		hs_k[,l] <- apply(pop_l_snps, 2, function(x){ sum( x == 1, na.rm=T) }) / ( ns_k[,l] );
	}

	## not sure what happens if the population size goes to 0 for a particular SNP?!
	dropouts <- which(ns_k == 0, arr.ind=T);
	if( length(dropouts) > 0 ){
		cat("Found:", length(dropouts), "snps that are (completely) imputed in one of the pops. Dropping it/them now.\n");
		pos_k <- pos_k[-dropouts];

		if( length(pos_k) == 0 ){
			next;
		}

		dropouts <- unique(dropouts[,1]);
		ps_k[-dropouts,];
		ns_k[-dropouts,];
		hs_k[-dropouts,];
	}

	fsts <- numeric(length(pos_k));
	for( m in 1:length(fsts)){
		fsts[m] <- Fst(populations$n, ns=ns_k[m,], ps=ps_k[m,], hs=hs_k[m,]);
	}

	results[[k]] <- cbind(chromosomeNumber, pos_k, round(fsts, 5));
	cat("\n");
}

results <- do.call(rbind, results);
colnames(results) <- c("chr", "pos", "fst");

if( exists("pop1") & exists("pop2")){
	outputFileName <- paste(usableNamePop1, "_vs_", usableNamePop2, ".chr", chromosomeNumber, ".fst.txt", sep="");

} else {
	outputFileName <- paste(prefix, "global.chr", chromosomeNumber, ".fst.txt", sep="");
}

write.table(results, outputFileName, quote=F, sep=",", row.names=F);
cat("Analyses finished and file written.\n");
quit(save="no");
