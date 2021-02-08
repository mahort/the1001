# TODO: Update to ggplot
# 
# Author: matt.horton
###############################################################################

rm(list=ls());

callMethodId <- 2;
offsetBetweenChromosomes <- 4e6;
lowerThreshold <- 0;
snpImputationThreshold <- 0.1;
targetWindowSize <- 1e4;
movingWindowAverage <- 1e5;

plotLrrs <- TRUE;
summaryStatistic <- "max";
summaryStatistic <- match.arg(summaryStatistic, c("median", "max", "mean"));

###############################################################################
# from here down, everything else is generic and can be used for whichever contrast.
#scoreColor <- c("dodgerblue3", "mediumpurple", "firebrick1", "darkolivegreen4", "tan1", "hotpink3", "cadetblue2", "darkslateblue", "orange3", "darksalmon", "forestgreen", "darksalmon", "indianred3");
pointColor <- "cadetblue3";
borderColor <- "cadetblue4";

###############################################################################
setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/global/results", sep=""));
files <- list.files(path=".", pattern="^adm");
files <- files[grep(summaryStatistic, files)];
files <- files[grep(paste("win", targetWindowSize, sep=""), files)];
files <- files[grep("txt$", files)];
print(files);



## drop the rug.
if( plotLrrs ){
	## read the file from resources and then put pipes under the global scores.
	lrrs <- read.table("~/projects/resources/tair_stuff/lrr_disease.models.pos.txt", header=T, sep="\t", as.is=T, stringsAsFactors=F);
	lrrs <- cbind(lrrs, midpoint=(lrrs[,"start"] + ((lrrs[,"stop"] - lrrs[,"start"])/2)));
}



averageFileName <- "adm_global.combined.fst.mean.win1e+05.scores.txt";
average <- read.table(averageFileName, header=T, sep="\t", as.is=T, stringsAsFactors=F);
average <- cbind(average, midpoint=(average[,"win_start"] + (average[,"win_end"] - average[,"win_start"])/2))

for( j in 1:length(files)){
	filename_j <- files[j];
	cat("Processing file: ", filename_j, "\n");
	
	scores1 <- read.table(filename_j, header=T, sep="\t", as.is=T, stringsAsFactors=FALSE);
	scores1 <- subset(scores1, score >= lowerThreshold);
	windowSize <- as.numeric(do.call(rbind, lapply(lapply(do.call(rbind, strsplit(filename_j, "win"))[,2], strsplit, "\\."), unlist))[,1]);

	## chromosomal endpoints as markers.
	scores1 <- cbind(scores1, midpoint=(scores1[,"win_start"] + (scores1[,"win_end"] - scores1[,"win_start"])/2));
	endPoints <- as.numeric(c(tapply(scores1[,"midpoint"], scores1[,"chr"], max)));
	maxEndPoint <- max(endPoints);

	## place the centromeres in space according to the offset (b/w chromosomes) and chromosome #
#	centro_start <- c(14364752, 3602775, 12674550, 2919690, 11668616);
#	centro_end   <- c(15750321, 3735247, 13674767, 4011692, 12082583);
	centro_start <- c(13.7, 2.45, 11.3, 1.8, 11)
	centro_end <- c(15.9, 5.5, 14.3, 5.15, 13.35);
	centro_start <- centro_start * 1e6;
	centro_end <- centro_end * 1e6;

	## make the image. 
	outputFileName <- gsub(".txt", paste(".allchr.ll", lowerThreshold, ".eps", sep=""), filename_j);

	setEPS(paper="special", horizontal=TRUE);
	postscript(outputFileName, width=5, height=9);

	layout(matrix(1:5, ncol=1))
	par(mar=c(5, 4, 0, 0));

	###############################################################################
	## plot the 5 separate chromosomes for each population
	for( k in 1:5 ){
		fsub_k <- subset(scores1, chr == k);

		plot.new();
		plot.window(xlim=c(0, max(scores1[,"win_end"])), ylim=c(lowerThreshold, 1.02));		
		axis(2);
	
		title(ylab=expression(F[ST]));
		if( k == 5 ){
			title(xlab="Mb");
		}

		tickMarks <- seq( 0, endPoints[k], by=4e6);
		chrLabels <- tickMarks / 1e6;
		chrLabels[seq(2, length(chrLabels), by=2)] <- "";
		
		rect( centro_start[k], lowerThreshold - 0.02, centro_end[k], 1.02, col="honeydew3", border="honeydew3");
		
		# x-axis offset!
		if( plotLrrs ){
			lrrs_k <- subset(lrrs, chr == k);
			rug( x=lrrs_k[,"midpoint"], col="brown3", lwd=1.4);
		}

		points(fsub_k[,"midpoint"], fsub_k[,"score"], col=borderColor, cex=0.2, pch=21, bg=pointColor);

		## add circles for the completely distinguished regions.
		ones <- subset(fsub_k, score == 1);
		points(ones[,"midpoint"], ones[,"score"], col="red", cex=0.6, pch=21, bg="red");
		axis(1, at=tickMarks, labels=chrLabels, cex.axis=0.85);
		cat("We've made a manhattan plot for the global fst analysis, chr:", k, "\n"); #subset\n");
	}

	dev.off();
}

setEPS(paper="special", horizontal=TRUE);
outputFileName <- gsub(".txt", paste(".allchr.ll", lowerThreshold, ".eps", sep=""), averageFileName);

postscript(outputFileName, width=5, height=9);

layout(matrix(1:5, ncol=1))
par(mar=c(5, 4, 0, 0));

for( k in 1:5 ){
	average_k <- subset(average, chr==k);
	
	plot.new();
	plot.window(xlim=c(0, max(average[,"win_end"])), ylim=c(lowerThreshold, max(average[,"score"], na.rm=T)));		
	axis(2);
	
	title(ylab=expression(F[ST]));
	if( k == 5){
		title(xlab="Mb");	
	}
	
	rect( centro_start[k], lowerThreshold - 0.03, centro_end[k], max(average_k[,"score"], na.rm=T) + 0.03, col="honeydew3", border="honeydew3");
	lines(average_k[,"midpoint"], average_k[,"score"]);
	axis(1, at=tickMarks, labels=chrLabels, cex.axis=0.85);
	
	# x-axis offset!
	if( plotLrrs ){
		lrrs_k <- subset(lrrs, chr == k);
		rug( x=lrrs_k[,"midpoint"], col="brown3", lwd=1.4);
	}
}

dev.off();
