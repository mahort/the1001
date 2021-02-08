# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());
assign("last.warning", NULL, envir = baseenv())

require(multicore);
numberOfCores <- 7;

window <- 1e4;
callMethodId <- 2;
snpImputationThreshold <- 0.1;

setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/global/results/", sep="")); #")
files <- list.files(path=".", pattern="combined.fst.txt$");
print(files);

summaryStatistic <- "min";
summaryStatistic <- match.arg(summaryStatistic, c("mean", "max", "median", "min"));

#####################################################################################################################
## make the windows for the gwas results.
#####################################################################################################################
makeWindows <- function( arg_set, phenotype, scores, method){
	
	chromosome <- arg_set$chr;
	windowSize <- arg_set$window;
	
	# quick sort just in case.
	subset <- scores[which(scores[,"chr"]==chromosome),];
	subset <- subset[order(subset[,"pos"]),];
	
	# specify a split size and then read these separately before combining the objects.
	windows <- seq(1, max(subset[,"pos"]), by=windowSize);
	windows <- c(windows, max(subset[,"pos"] + 2e4));
	windowIds <- cut(subset[,"pos"], breaks=windows, labels=FALSE);
	
	uniqueWindows <- unique(windowIds);
	results <- list();
	
	for( j in 1:length(uniqueWindows)){
		
		targetWindow <- uniqueWindows[j];
		
		# uniqueWindows is called on windowIds, so it should be impossible for windowIds to not have the targetWindow in it.
		indices <- which( windowIds == targetWindow );
		
		if( j %% 100 == 0 ){ cat("On window:", j, "of", length(uniqueWindows),"\n"); }
		
		result_id <- paste("results_", j, sep="");
		vals <- subset[indices, "fst"];
		
		if( sum(is.na(vals)) == length(vals)){
			topSnp <- subset[indices[1], "pos"];
			topScore <- 0;
			
		} else {
			topSnp <- which.max(vals);
			topSnp <- subset[indices[topSnp], "pos"];

			if( method == "min" ){
				topScore <- -log(min(vals, na.rm=T), 10);

			} else if( method == "mean" ){
				topScore <- mean(vals, na.rm=T);	

			} else if( method == "median" ){
				topScore <- median(vals, na.rm=T);

			} else if( method == "max" ){
				topScore <- max(vals, na.rm=T);
			}
			
		}
		
		results[[result_id]] <- c(
				j, 
				phenotype, 
				chromosome,
				windows[targetWindow], 
				windows[targetWindow+1] - 1,
				topSnp,
				topScore);
	}
	
	output <- do.call(rbind,results);
	if( !is.null(warnings())){
		stop("Error for chr:", chromosome, "and phenotype:", phenotype, "\n");
	}
	
	colnames(output) <- c("window_id", "phenotype", "chr", "win_start", "win_end", "position", "score");
	return(output);
}

argList <- lapply(1:5, function(x){ 
			list( chr = x[1], window = window );}); # make a list out of it.


## need to rename the file and read it...
for( j in 1:length(files)){

	filename_j <- files[j];
	cat("Opening:", filename_j, "\n");

	file_j <- read.table(filename_j, header=T, sep=",", as.is=T, stringsAsFactors=F);
	outputFileName <- gsub(".txt", paste(".", summaryStatistic, ".win", window, ".scores.txt", sep=""), filename_j);
	
	resultSet <- mclapply(
			argList,
			makeWindows,
			paste(summaryStatistic, "_fst", sep=""),
			file_j,
			method=summaryStatistic,
			mc.cores=numberOfCores );

	windowedPvalues <- do.call(rbind, resultSet);
	if( !is.null(warnings())){
		stop("Error for file:", filename_j, "\n");
	}

	cat("Writing windowed-results for file:", outputFileName, "\n");
	write.table(windowedPvalues, outputFileName, quote=F, sep="\t", row.names=F);
}






