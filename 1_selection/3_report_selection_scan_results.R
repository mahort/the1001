# TODO: Add comment
#
# Author: matt.horton
###############################################################################
# one reasonable way to review the results is to generate a report file per p-values list, using either a p-value or n-count cutoff.
# 1.) specify a regex to match the p-value files
# 2.) specify a cutoff
# 3.) enumerate over the files that we want to generate reports for
# 4.) consider the snps in turn, in new matrix, record score, URL and gene models close by (using gff files).
# 5.) write the report.

rm(list=ls());

selectionScoreColumn <- "fst"; #"Likelihood"; #LR #omega
minimumResults <- 10;
numberOfNeighbors <- 1;
callMethodId <- 2;

snpImputationThreshold <- 0.1;

##########################################################################################################################################################################
## GET THE ANNOTATION INFO
##########################################################################################################################################################################
filename <- "~/projects/resources/tair_stuff/genes_final.txt";
genes <- read.table(filename, header=T, sep="\t", as.is=T, stringsAsFactors=FALSE, quote="");

baseUrl <- "http://www.ncbi.nlm.nih.gov/projects/mapview/maps.cgi?TAXID=3702&CHR=TARGETCHROMOSOME&MAPS=cntg-r,clone,tair_marker,genes[POSITIONSTART%3APOSITIONEND]&amp;ZOOM=0.1";
urlBuffer <- 1000;

functional <- read.table("~/projects/resources/tair_stuff/TAIR10_functional_descriptions",header=T,sep="\t",as.is=T,stringsAsFactors=FALSE,quote="",comment="");

focal <- functional[grep("cold|cbf|arr|myb|lrr|defens|disease|resistan|antimicrob|antifung", functional[,5], ignore.case=T),];
focal.ids <- strsplit(focal[,1], "\\.");
focal.ids <- unique(do.call(rbind, focal.ids)[,1]);

isItInteresting <- function(candidate){
	return(length(which( candidate %in% focal.ids)) != 0);
}

##########################################################################################################################################################################
## OPEN THE FILES.
##########################################################################################################################################################################
delimiter = '\t';

if( selectionScoreColumn == "fst"){
	## we need to choose the file up here, but for now...
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/", sep=""));
	regexPattern <- "^fst_scores|^global|^adm|^popn|relict";
	files <- list.files(path=".", pattern=regexPattern, recursive=T);
	files <- files[grep("txt$", files)];
	files <- files[grep("combined", files)];
	files <- files[grep("win", invert=T, files)];
	delimiter = ',';

} else if( selectionScoreColumn == "LR" ){
	setwd(paste("~/projects/the_1001/selection/clr/", callMethodId, "/", sep=""));
	regexPattern <- "clr.combined.txt$";
	files <- list.files(path=".", pattern=regexPattern, recursive=T);

} else if( selectionScoreColumn == "omega" ){
	regexPattern <- "omega_scores.txt$";

	setwd(paste("~/projects/the_1001/selection/omega/", callMethodId, "/imp", snpImputationThreshold, "/", sep=""));	
	files <- list.files(path=".", pattern=regexPattern, recursive=T);
}

files <- files[grep("eps$|chr|report|goterm$", invert=T, files)];
print(files);

for( k in 1:length(files)){

	file_k <- read.table(files[k], header=T, sep=delimiter, as.is=T, stringsAsFactors=F); ## probably have to make the delimiter flexible for clr

	file_k <- file_k[which(!is.na(file_k[,selectionScoreColumn])),];
	file_k <- file_k[order(file_k[,selectionScoreColumn], decreasing=T),];
	topScores <- file_k[1:minimumResults,];

	##########################################################################################################################################################################
	## ITERATE OVER THE SCORES.
	##########################################################################################################################################################################
	results <- list();
	for( j in 1:nrow(topScores)){
		
		candidate <- topScores[j,];
		result_id <- paste("results_",j,sep="");
		results[[result_id]][["filename"]] <- files[k];
		results[[result_id]][["chr"]] <- chr <-  candidate[1,1];
		results[[result_id]][["pos"]] <- pos <- candidate[1,2];
		results[[result_id]][["score"]] <- candidate[1, selectionScoreColumn];

		# what are the genes in this area???
		genes_chr_i <- genes[which(genes[,"chr"] == chr),];

		category_index <- which(genes_chr_i[,"start"] <= pos & genes_chr_i[,"stop"] >= pos);
		numberOfIndices <- length(category_index);
		if( numberOfIndices != 1 ){
			# different symbols for the same model. take the first one.
			category_index <- category_index[1];
		}

		# determine (2) 5 prime features.
		neighbor_index <- category_index;
		neighbors5prime <- c();
		neighborsDesc <- c();
		while( neighbor_index > 1 ){
			neighbor_index <- neighbor_index - 1;
			if( genes_chr_i[neighbor_index,"feature"] != "intergenic" ){
				# record it.
				interesting <- isItInteresting(genes_chr_i[neighbor_index, "ids"]);
				symbol <- genes_chr_i[neighbor_index, "symbol"];

				if( is.na(symbol)){
					symbol <- paste(genes_chr_i[neighbor_index,"ids"], " (", genes_chr_i[neighbor_index,"feature"], ")", sep="");
				}

				desc <- genes_chr_i[neighbor_index,"full_name"];
				if( interesting ){
					cat("Found something interesting.\n");
					desc <- paste(desc, "; [annotation-support]", sep="");
				}

				neighbors5prime <- c(symbol, neighbors5prime);
				neighborsDesc <- c(desc, neighborsDesc);
				if( length(neighbors5prime) == numberOfNeighbors ){
					break;
				}
			}
		}

		results[[result_id]][["neighbors_5prime"]] <- paste(neighbors5prime, collapse=",");
		results[[result_id]][["neighbors_5prime_desc"]] <- paste(neighborsDesc, collapse=";");
		
		if( genes_chr_i[category_index, "feature"] != "intergenic" ){
			symbol <- genes_chr_i[category_index, "symbol"];
			if( is.na(symbol)){
				symbol <- paste(genes_chr_i[category_index,"ids"], " (", genes_chr_i[category_index,"feature"], ")", sep="");
			}
			
			results[[result_id]][["name"]] <- symbol;
			interesting <- isItInteresting(genes_chr_i[category_index, "ids"]);
			desc <- genes_chr_i[category_index, "full_name"];
			
			if( interesting ){
				cat("Found something interesting.\n");
				desc <- paste( desc, "; [annotation-support]", sep="");			
			}
			
			results[[result_id]][["desc"]] <- desc;
			
		} else {
			results[[result_id]][["name"]] <- "intergenic";
			results[[result_id]][["desc"]] <- NA;
		}
		
		# determine (2) 3 prime features.
		neighbor_index <- category_index + numberOfIndices - 1;
		neighbors3prime <- c();
		neighborsDesc <- c();

		while( neighbor_index < nrow(genes_chr_i)){
			neighbor_index <- neighbor_index + 1;
			if( genes_chr_i[neighbor_index,"feature"] != "intergenic" ){
				# record it.

				interesting <- isItInteresting(genes_chr_i[neighbor_index, "ids"]);
				symbol <- genes_chr_i[neighbor_index,"symbol"];
				
				if( is.na(symbol)){
					symbol <- paste(genes_chr_i[neighbor_index,"ids"], " (", genes_chr_i[neighbor_index,"feature"], ")", sep="");
				}
				
				desc <- genes_chr_i[neighbor_index,"full_name"];
				if( interesting ){
					cat("Found something interesting.\n");
					desc <- paste(desc, "; [annotation-support]", sep="");
				}

				neighbors3prime <- c(neighbors3prime, symbol);
				neighborsDesc <- c(neighborsDesc, desc);

				if( length(neighbors3prime) == numberOfNeighbors ){
					break;
				}
			}
		}

		results[[result_id]][["neighbors_3prime"]] <- paste(neighbors3prime, collapse=",");
		results[[result_id]][["neighbors_3prime_desc"]] <- paste(neighborsDesc, collapse=";");

		# make it easier to use mapview.
		url <- gsub( "TARGETCHROMOSOME", candidate[1], baseUrl);
		url <- gsub( "POSITIONSTART", candidate[2] - urlBuffer, url);
		url <- gsub( "POSITIONEND", candidate[2] + urlBuffer, url);

		results[[result_id]][["url"]] <- url;
	}

	results <- do.call(rbind,results);

	outputFileName <- gsub("txt", "report.txt", files[k]);
	write.table(results, outputFileName, quote=F,sep="\t", row.names=F);
	cat("Wrote report for file:", files[k], "\n");
}