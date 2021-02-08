# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());

###############################################################################
## user parameters, which act to identify the files we'll concatenate.
###############################################################################
imputationThreshold <- 0.1;

callMethodId <- 2;
gridInterval <- 500; ## grid size, not selection scan size.
selectionScanMinWindowStart <- 1e4;
selectionScanMaxWindowStop <- 1e5;
range <- paste(selectionScanMinWindowStart, selectionScanMaxWindowStop, sep="-");
populationSize <- 1135;

###############################################################################
## the population def file
###############################################################################
uniquePopulations <- "global";

for( k in 1:length(uniquePopulations)){

	population_k <- uniquePopulations[k];
	usableName_k <- gsub(" " , "_", population_k);
	cat("Finding files for population: ", population_k, "\n");
	setwd(paste("~/projects/the_1001/selection/omega/", callMethodId, "/imp", imputationThreshold, "/", usableName_k, sep=""));

	###############################################################################
	## there are 5 chr[1-5] directories that contain our files.
	## and within those 5 directories are separate 'slices' (the chr-specific files are sliced up to allow for faster processing)
	files <- list.files(path=".", pattern="Report", recursive=TRUE);
	files <- files[grep(populationSize, files)];
	files <- files[grep(paste("grd", gridInterval, "\\b", sep=""), files)];
	files <- files[grep(paste("win", gsub("e+", "e\\+", range, fixed=T), "\\b",  sep=""), files)];
	
	stopifnot(length(files) == 5);
	print(files);
	###############################################################################	

	omega_combined <- list();
	for( chr_i in 1:5 ){
		
		files.chr_i <- files[grep(paste("chr", chr_i, sep=""), files)];
		
		stopifnot( length(files.chr_i) == 1 );
		cat("Opening singlet file:", files.chr_i, "\n");
		chrNumber <- as.numeric(unlist(strsplit(unlist(strsplit(files.chr_i, "chr"))[2], "/"))[1]);
		
		data <- read.table(files.chr_i, sep="\t", as.is=T, stringsAsFactors=FALSE, header=F, skip=2);
		omega_combined[[files.chr_i]] <- cbind(chr=rep(chrNumber, nrow(data)), data);
	}
	
	omega_combined <- do.call(rbind, omega_combined);
	colnames(omega_combined) <- c("chr", "position", "omega", "win_start", "win_end", "valid");
	omega_combined <- omega_combined[order(omega_combined[,"chr"], omega_combined[,"position"]),];
	omega_combined <- omega_combined[omega_combined$valid == 1,];
	
	###############################################################################
	## now mask the centromeres.
	###############################################################################
	centro_start <- c(14364752, 3602775, 12674550, 2919690, 11668616);
	centro_end   <- c(15750321, 3735247, 13674767, 4011692, 12082583);
	###############################################################################
#
#	for( chr_i in 1:5 ){
#		indices <- which(omega_combined[,"chr"] == chr_i & omega_combined[,"position"] >= centro_start[chr_i] & omega_combined[,"position"] <= centro_end[chr_i]);
#		if( length(indices) > 0 ){
#			cat("Masking", length(indices), "SNPs in centromere:", chr_i, "\n");
#			omega_combined[indices, "omega"] <- 3.7e-06;	
#		}
#	}
	
	outputFileName <- paste(usableName_k, ".n_", populationSize, ".cm", callMethodId, unlist(strsplit(basename(files[1]), "omega"))[2], ".omega_scores.txt", sep="");
	write.table(omega_combined, outputFileName, quote=F, sep="\t", row.names=F);
}


