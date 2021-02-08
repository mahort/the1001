# TODO: Update to ggplot
# 
# Author: matt.horton
###############################################################################

rm(list=ls());

callMethodId <- 2;
offsetBetweenChromosomes <- 4e6;
lowerThreshold <- 0.485;
snpImputationThreshold <- 0.1;

plotLrrs <- TRUE;

###############################################################################
# from here down, everything else is generic and can be used for whichever contrast.
#scoreColor <- c("dodgerblue3", "mediumpurple", "firebrick1", "darkolivegreen4", "tan1", "hotpink3", "cadetblue2", "darkslateblue", "orange3", "darksalmon", "forestgreen", "darksalmon", "indianred3");
pointColor <- "cadetblue3";
borderColor <- "cadetblue4";

dataset <- "global";
dataset <- match.arg(dataset, c("global", "relict", "padm", "ppopn"));

###############################################################################
if( dataset == "global" ){
	regexPattern <- "combined.fst.txt$";
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/", dataset, "/results/", sep=""));

} else if( dataset == "relict" ){
	regexPattern <- "combined.fst.txt$";
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/pairwise/admixed/results/", sep=""));

} else if( dataset == "padm" ){
	regexPattern <- "combined.txt$";
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/pairwise/admixed/results/", sep=""));

} else if( dataset == "ppopn" ){
	regexPattern <- "combined.txt$";
	setwd(paste("~/projects/the_1001/selection/fst/", callMethodId, "/imp", snpImputationThreshold, "/pairwise/popn/results/", sep=""));
}


## drop the rug.
if( plotLrrs ){
	## read the file from resources and then put pipes under the global scores.
	lrrs <- read.table("~/projects/resources/tair_stuff/lrr_disease.models.pos.txt", header=T, sep="\t", as.is=T, stringsAsFactors=F);
	lrrs <- cbind(lrrs, midpoint=(lrrs[,"start"] + ((lrrs[,"stop"] - lrrs[,"start"])/2)));
}

files <- list.files(path=".", pattern=regexPattern);
print(files);

for( j in 1:1 ){ #length(files)){
	filename_j <- files[j];
	cat("Processing file: ", filename_j, "\n");

	scores1 <- read.table(filename_j, header=T, sep=",", as.is=T, stringsAsFactors=FALSE);

	if(exists("lowerThreshold")){
		lowerThreshold_j <- lowerThreshold;

	} else {
		lowerThreshold_j <- quantile(scores1[,"fst"], probs=0.75);	
	}
	
	scores1 <- subset(scores1, fst >= lowerThreshold_j);
	
	## chromosomal endpoints as markers.
	endPoints <- as.numeric(c(0, tapply(scores1[,"pos"], scores1[,"chr"], max)));
	
	## place the centromeres in space according to the offset (b/w chromosomes) and chromosome #
	centro_start <- c(14364752, 3602775, 12674550, 2919690, 11668616);
	centro_end   <- c(15750321, 3735247, 13674767, 4011692, 12082583);
	centro_start <- centro_start + cumsum(endPoints[1:length(centro_start)]) + (0:4)*offsetBetweenChromosomes;
	centro_end <- centro_end + cumsum(endPoints[1:length(centro_end)]) + (0:4)*offsetBetweenChromosomes;
	
	## make the image. 
	if( exists("lowerThreshold")){
		outputFileName <- gsub(".txt", paste(".ll", lowerThreshold, ".eps", sep=""), filename_j);

	} else {
		outputFileName <- gsub(".txt", paste(".llq3.eps", sep=""), filename_j);
	}
	

	setEPS(paper="special", horizontal=TRUE);
	postscript(outputFileName, width=9, height=5, bg="white");

	plot.new();
	plot.window(xlim=c(0, sum(endPoints) + (4*offsetBetweenChromosomes)), ylim=c(lowerThreshold_j, max(scores1[,"fst"], na.rm=T)));
	par(mar=c(5, 4, 1, 1));
	
	axis(2);
	title(ylab=expression(F[ST]));
	title(xlab="Mb");

	tickMarks <- c();
	labels <- c();

	###############################################################################
	## plot the 5 separate chromosomes for each population
	for( k in 1:5 ){
		fsub_k <- subset(scores1, chr == k);

		breaks <- seq( 0, endPoints[k + 1], by=4e6);
		chrLabels <- breaks / 1e6;
		chrLabels[seq(2, length(chrLabels), by=2)] <- "";
		
		labels <- c(labels, chrLabels);
		tickMarks <- c(tickMarks, sum(endPoints[1:k]) + breaks + (k-1)*offsetBetweenChromosomes);

		# x-axis offset!
		if( plotLrrs ){
			lrrs_k <- subset(lrrs, chr == k);
			lrrs_k[,"midpoint"] <- lrrs_k[,"midpoint"] + sum(endPoints[1:k]) + (k-1)*offsetBetweenChromosomes;
			rug( x=lrrs_k[,"midpoint"], col="brown3", lwd=1.4);
		}

		fsub_k[,"pos"] <- fsub_k[,"pos"] + sum(endPoints[1:k]) + (k-1)*offsetBetweenChromosomes;
		points(fsub_k[,"pos"], fsub_k[,"fst"], col=borderColor, cex=0.2, pch=21, bg=pointColor);

	}

	axis(1, at=tickMarks, labels=labels, cex.axis=0.85);
	box();

	dev.off();
	cat("We've made a manhattan plot for the global fst analysis.\n"); #subset\n");
}
