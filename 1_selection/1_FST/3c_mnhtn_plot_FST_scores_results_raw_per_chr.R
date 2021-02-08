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
	scores1 <- subset(scores1, fst >= lowerThreshold);
	
	## chromosomal endpoints as markers.
	endPoints <- as.numeric(tapply(scores1[,"pos"], scores1[,"chr"], max));
	maxEndPoint <- max(endPoints);

	## place the centromeres in space according to the offset (b/w chromosomes) and chromosome #
	centro_start <- c(14364752, 3602775, 12674550, 2919690, 11668616);
	centro_end   <- c(15750321, 3735247, 13674767, 4011692, 12082583);

	###############################################################################
	## plot the 5 separate chromosomes for each population
	for( k in 1:5 ){
		fsub_k <- subset(scores1, chr == k);

		## make the image. 
		outputFileName <- gsub(".txt", paste(".chr", k, ".ll", lowerThreshold, ".eps", sep=""), filename_j);

		setEPS(paper="special", horizontal=TRUE);
		postscript(outputFileName, width=9*(endPoints[k]/maxEndPoint), height=5, bg="white");
	
		plot.new();
		plot.window(xlim=c(0, endPoints[k]), ylim=c(lowerThreshold, max(fsub_k[,"fst"], na.rm=T)));
		
		par(mar=c(5, 4, 1, 1));
		
		axis(2);
		title(ylab=expression(F[ST]));
		title(xlab="Mb");
	
		tickMarks <- c();
		labels <- c();

		breaks <- seq( 0, endPoints[k], by=4e6);
		chrLabels <- breaks / 1e6;
		chrLabels[seq(2, length(chrLabels), by=2)] <- "";

		tickMarks <- c(breaks);

		# x-axis offset!
		if( plotLrrs ){
			lrrs_k <- subset(lrrs, chr == k);
			rug( x=lrrs_k[,"midpoint"], col="brown3", lwd=1.4);
		}

		points(fsub_k[,"pos"], fsub_k[,"fst"], col=borderColor, cex=0.2, pch=21, bg=pointColor);
		axis(1, at=tickMarks, labels=chrLabels, cex.axis=0.85);
		box();
		
		dev.off();
		cat("We've made a manhattan plot for the global fst analysis, chr:", k, "\n"); #subset\n");
	}

}
