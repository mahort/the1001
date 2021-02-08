# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());

###############################################################################
## go to the directory, find the contrasts that match, and combine them...
###############################################################################
callMethodId <- 2;
snpImputationThreshold <- 0.1;

admixedFile <- FALSE;
dataset <- "pairwise"; ## when relict, we ignore the admixedFile variable.
dataset <- match.arg(dataset, c("global", "pairwise", "relict", "pbsrelicts"));

################################################################################
###
################################################################################
regexPattern <- NULL;
if( dataset == "relict" | dataset == "pbsrelicts" ){
	prefix <- "geo-nonrelicts_vs_Relicts";

	if( dataset == "pbsrelicts" ){
		prefix <- c(prefix, "french-nonrelicts_vs_Relicts", "french-nonrelicts_vs_geo-nonrelicts"); #"french-nonrelicts_vs_([a-zA-Z0-9\\-])+");
	}

	regexPattern <- unlist(lapply(lapply(prefix, paste, ".chr[1-5].fst.txt$", sep=""), unlist));
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/pairwise/admixed/results/", sep=""));

} else if( dataset == "global" ){
	if( admixedFile ){
		prefix <- "adm_global";

	} else {
		prefix <- "popn_global";
	}

	regexPattern <- paste(prefix, ".chr[1-5].fst.txt$", sep="");
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/global/results/", sep=""));
	
} else if( dataset == "pairwise" ){
	if( admixedFile ){
		accessionMap <- read.table("~/projects/the_1001/selection/fst/admix_pop_metadata/metadata_structure.txt", header=T, sep="\t", as.is=T, stringsAsFactors=F);
		colnames(accessionMap)[which(colnames(accessionMap) == "tg_ecotypeid")] <- "ecotype_id";
		colnames(accessionMap)[which(colnames(accessionMap) == "Admixed")] <- "population_name";
		setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/pairwise/admixed/results/", sep=""));

	} else {
		accessionMap <- read.table("~/projects/the_accessions/the1001.genomes.txt", header=T, sep="\t", as.is=T, stringsAsFactors=F, quote="");
		setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/pairwise/popn/results/", sep=""));
	}
}

################################################################################
### use the accession map...
################################################################################
if( !is.null(regexPattern)){ ## just 5
	for( k in 1:length(regexPattern)){
		files <- list.files(path=".", pattern=regexPattern[k]);
		stopifnot(length(files) == 5);

		cat("Preparing to combine files for:", prefix[k], "\n");
		results <- list();
		for( l in 1:length(files)){
			file_l <- files[l];
			results[[l]] <- read.table(file_l, header=T, sep=",", as.is=T, stringsAsFactors=F);
		}
		
		results <- do.call(rbind, results);
		results <- results[order(results[,"chr"], results[,"pos"]),];
		results <- results[which(!is.na(results[,"fst"])),];
		
		filename <- paste(prefix[k], ".combined.fst.txt", sep="");
		write.table(results, file=filename, quote=F, sep=",", row.names=F);
		cat("Files from:", prefix[k], "combined.\n");
	}

} else { ## we have to iterate and find our samples.

	uniquePopulations <- unique(accessionMap[,"population_name"]);
	uniquePopulations <- sort(uniquePopulations[!is.na(uniquePopulations)]);
	
	for( j in 1:(length(uniquePopulations) - 1)){
		pop1 <- gsub(" ", "_", uniquePopulations[j]);
		
		for( k in (j+1):length(uniquePopulations)){
			pop2 <- gsub(" ", "_", uniquePopulations[k]);
			
			pairwiseRegexPattern <- paste(pop1, "_vs_", pop2, sep="");
			files <- list.files(path=".", pattern=pairwiseRegexPattern);
			files <- files[grep("combined", invert=T, files)];
			stopifnot(length(files) == 5);
			
			cat("combining files from:", pairwiseRegexPattern, "\n");		
			results <- list();
			for( l in 1:length(files)){
				file_l <- files[l];
				results[[l]] <- read.table(file_l, header=T, sep=",", as.is=T, stringsAsFactors=F);
			}

			results <- do.call(rbind, results);
			results <- results[order(results[,"chr"], results[,"pos"]),];
			results <- results[which(!is.na(results[,"fst"])),];

			filename <- paste("fst_scores.", pairwiseRegexPattern, ".combined.txt", sep="");
			write.table(results, file=filename, quote=F, sep=",", row.names=F);
		}
	}
}