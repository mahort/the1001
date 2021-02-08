# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());

callMethodId <- 2;
setwd("~/projects/the_1001/selection/");

snpImputationThreshold <- 0.1;
windowSize <- 1e4;

###############################################################################
## 1.) determine the tests
###############################################################################
{	
	chooseFile <- function(selectionScanPath, pattern){
		
		if( pattern[1] != "" ){
			files <- list.files(path=selectionScanPath, pattern=pattern[1]);
			if( length(pattern) > 1 ){
				for( k in 2:length(pattern)){
					files <- files[grep(pattern[k], files)];
				}
			}
			
		} else {
			files <- list.files(path=selectionScanPath, pattern="txt$");
		}
		
		fileChoices <- paste(1:length(files), files, sep=": ");
		
		while( TRUE ){
			cat("Available files:\n");
			print(fileChoices);
			
			fileNumber <- readline("Please choose a file ");
			if( fileNumber == "" ){
				cat("You didn't choose anything.\n");
				
			} else {
				
				fileNumber <- as.numeric(fileNumber);
				if( fileNumber >= 1 | fileNumber <= length(files)){
					cat("You chose file:", files[fileNumber], "\n");
					break;

				} else {
					cat("Please choose a real file.\n");
				}
			}
		}
		
		return(list("base_dir"=selectionScanPath, "filename"=files[fileNumber]));
		
	};
	
	###############################################################################
	## 1.) determine the tests
	###############################################################################
	## combining the results from global scans for FST, omega+, and CLR
	## 1.) FST
	fstResultsPath <- paste(getwd(), "/fst/", callMethodId, "/imp", snpImputationThreshold, "/global/results/", sep="");
	omegaResultsPath <- paste(getwd(), "/omega/", callMethodId, "/imp", snpImputationThreshold, "/global/", sep="");
	clrResultsPath <- paste(getwd(), "/clr/", callMethodId, "/imp", snpImputationThreshold, "/global/", sep="");
	
	results <- list();
	results[["fst"]] <- chooseFile(selectionScanPath=fstResultsPath, pattern=paste("win", windowSize, ".scores.txt$", sep=""));
	results[["omega"]] <- chooseFile(selectionScanPath=omegaResultsPath, pattern=paste("win", windowSize, ".scores.txt$", sep=""));
	results[["clr"]] <- chooseFile(selectionScanPath=clrResultsPath, pattern=paste("win", windowSize, ".scores.txt$", sep=""));

}

###############################################################################
## 2.) join the files
###############################################################################
for( j in 1:length(results)){
	test_j <- names(results)[j];
	setwd(results[[j]][['base_dir']]);
	dataset <- read.table(results[[j]][['filename']], header=T, sep='\t', as.is=T, stringsAsFactors=F);
	dataset <- cbind(dataset, rank_ = rank(-1*dataset[,"score"]));
	dataset <- cbind(dataset, emp_ = dataset[,"rank_"] / (nrow(dataset)));

	colnames(dataset)[which(colnames(dataset) == 'rank_')] <- paste(test_j, "_rank", sep="");
	colnames(dataset)[which(colnames(dataset) == 'emp_')] <- paste(test_j, "_empP", sep="");
	colnames(dataset)[which(colnames(dataset) == 'score')] <- paste(test_j, "_score", sep="");
	colnames(dataset)[which(colnames(dataset) == 'position')] <- paste(test_j, "_position", sep="");
	dataset <- dataset[,which(!colnames(dataset) %in% c("window_id", "phenotype"))];
	results[[j]][['dataset']] <- dataset;
}

tmp <- results[[1]][['dataset']];
for( j in 2:length(results)){
	tmp <- merge(tmp, results[[j]][['dataset']], by=c("chr", "win_start", "win_end"))
}

output <- tmp;
output <- output[order(output[,"chr"], output[,"win_start"]),];

indices <- grep("empP", colnames(output));
chisq <- -2*rowSums(sapply(output[,indices], log));
chisq.p <- 1 - pchisq(chisq, (2*length(indices)));
output <- cbind(output, chisq=chisq, chisq_p=chisq.p);

setwd("~/projects/the_1001/selection/");
write.table(output, "the1001.selection_scan_ranks.txt", quote=F, sep="\t", row.names=F);

