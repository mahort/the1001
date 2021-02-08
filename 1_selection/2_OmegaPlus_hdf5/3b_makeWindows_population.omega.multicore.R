# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());
require(multicore);

###############################################################################
windowSize <- 1e4;
numberOfCores <- 7;

###############################################################################
## user parameters, which act to identify the files we'll concatenate.
###############################################################################
snpImputationThreshold <- 0.1;

callMethodId <- 2;
gridInterval <- 500; ## grid size, not selection scan size.

populationSize <- 1135;
selectionScanMinWindowStart <- 1e4;
selectionScanMaxWindowStop <- 1e5;
range <- paste(selectionScanMinWindowStart, selectionScanMaxWindowStop, sep="-");

#####################################################################################################################
## make the windows for the gwas results.
#####################################################################################################################
makeWindows <- function( arg_set, phenotype, scores){
	
	chromosome <- arg_set$chr;
	windowSize <- arg_set$window;
	
	# quick sort just in case.
	subset <- scores[which(scores[,"chr"]==chromosome),];
	subset <- subset[order(subset[,"position"]),];
	
	# specify a split size and then read these separately before combining the objects.
	windows <- seq(1, max(subset[,"position"]), by=windowSize);
	windows <- c(windows, max(subset[,"position"] + 2e4));
	windowIds <- cut(subset[,"position"], breaks=windows, labels=FALSE);
	
	uniqueWindows <- unique(windowIds);
	results <- list();
	
	for( j in 1:length(uniqueWindows)){
		targetWindow <- uniqueWindows[j];
		indices <- which( windowIds == targetWindow );
		
		if( j %% 100 == 0 ){
			cat("On window:", j, "of", length(uniqueWindows),"\n");
		}
		
		result_id <- paste("results_", j, sep="");
		vals <- subset[indices, "omega"];
		
		topSnp <- which.max(vals);
		topSnp <- subset[indices[topSnp], "position"];
		topScore <- max(vals);
		
		results[[result_id]] <- c(
				j, 
				phenotype, 
				chromosome,
				windows[targetWindow], 
				windows[targetWindow+1]-1,
				topSnp,
				topScore);
	}
	
	output <- do.call(rbind,results);
	colnames(output) <- c("window_id", "phenotype", "chr", "win_start", "win_end", "position", "score");
	return(output);
}

###############################################################################
## identify the populations we're collapsing on...
###############################################################################
metadata <- "~/projects/the_1001/selection/fst/admix_pop_metadata/metadata_structure.txt";
metadata <- read.table(metadata, header=T, sep="\t", as.is=T, stringsAsFactors=F);
colnames(metadata)[which(colnames(metadata) == "Admixed")] <- "population_name";
colnames(metadata)[which(colnames(metadata) == "tg_ecotypeid")] <- "ecotype_id";
#metadata <- subset(metadata, population_name %in% subpops);

subpops <- unique(metadata[,"population_name"]);
uniquePopulations <- subpops[which(!is.na(subpops))];
uniquePopulations <- uniquePopulations[uniquePopulations != "Relicts"];

for( k in 1:length(uniquePopulations)){

	population_k <- uniquePopulations[k];
	usableName_k <- gsub(" " , "_", population_k);
	cat("Finding files for population: ", population_k, "\n");
	setwd(paste("~/projects/the_1001/selection/omega/", callMethodId, "/imp", snpImputationThreshold, "/", usableName_k, sep=""));

	###############################################################################
	## identify the files we're windowing...
	###############################################################################
	files <- list.files(path=".", recursive=TRUE, pattern="omega_scores.txt");
	files <- files[grep(".win.scores.txt", files, invert=T)];
	files <- files[grep(gsub("+", "\\+", range, fixed=T), files)];
	populationSize <- as.numeric(unlist(strsplit(unlist(strsplit(files[1], "\\.n_"))[2], "\\."))[1]);
	files <- files[grep(populationSize, files)];

	stopifnot(length(files) == 1);
	print(files);
	
	###############################################################################	
	argList <- lapply(1:5, function(x){ 
				list( chr = x[1], window = windowSize );}); # make a list out of it.
	
	for( i in 1:length(files)){
		filename <- files[i];
		outputFileName <- gsub(".txt$", paste(".win", windowSize, ".scores.txt", sep=""), filename);
		cat("Preparing to window results for file:", filename, "\n");
		file <- read.table(filename, header=T, sep="\t");
		file[,"position"] <- floor(file[,"position"]);
		
		resultSet <- mclapply(
				argList,
				makeWindows,
				filename,
				file,
				mc.cores=numberOfCores,
				mc.preschedule=FALSE, 
				mc.cleanup=TRUE);
		
		windowedPvalues <- do.call(rbind, resultSet);
		write.table(windowedPvalues, outputFileName, quote=F, sep="\t", row.names=F);
	}
}
