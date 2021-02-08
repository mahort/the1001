# TODO: Add comment
# 
# Author: matt.horton
###############################################################################


rm(list=ls());

setwd("~/projects/the_1001/selection/omega/");

callMethodId <- 2;
imputationThreshold <- 0.1;
windowSize <- 1e4;
baselineOmega <- 0.3;

###############################################################################
## parameters
###############################################################################
gridInterval <- 500; ## grid size, not selection scan size.
selectionScanMinWindowStart <- 1e4;
selectionScanMaxWindowStop <- 1e5;
range <- paste(selectionScanMinWindowStart, selectionScanMaxWindowStop, sep="-");

legendPositionScalar <- 1;
offsetBetweenChromosomes <- 4e6;

#pop.colors <- c("darksalmon", "forestgreen", "firebrick1", "darkolivegreen4", "tan1", "hotpink3", "cadetblue2", "darkslateblue", "orange3"); # "darksalmon", "forestgreen", "darksalmon" );
pointColor <- "cadetblue3";
borderColor <- "cadetblue4";

###############################################################################
## open the datasets, and plot.
###############################################################################
setwd(paste("~/projects/the_1001/selection/omega/", callMethodId, "/imp", imputationThreshold, "/", sep=""));
regexPattern <- paste(".omega_scores.win", windowSize, ".scores.txt", sep="");

files <- list.files(path=".", pattern=regexPattern, recursive=T);
files <- files[grep(gsub("+", "\\+", range, fixed=T), files)];
print(files);

################################################################################
## from here down, everything else is generic and can be used for whichever contrast.
## determine the population specific colors.
################################################################################

for( l in 1:length(files)){
	file_l <- files[l];
	scores1 <- read.table(file_l, header=T, sep="\t", as.is=T, stringsAsFactors=FALSE);
	scores1 <- cbind(scores1, midpoint=(scores1[,"win_start"] + scores1[,"win_end"]) / 2);
	scores1 <- subset(scores1, score >= baselineOmega);

	## determine the population under survey...
	focalPopulation <- unlist(strsplit(file_l, "/"))[1];
	populationSize <- as.numeric(unlist(strsplit(unlist(strsplit(file_l, "n_"))[2], "\\."))[1]);

	## chromosomal endpoints as markers.
	endPoints <- as.numeric(c(0, tapply(scores1[,"midpoint"], scores1[,"chr"], max)));

	## place the centromeres in space according to the offset (b/w chromosomes) and chromosome #
	centro_start <- c(14364752, 3602775, 12674550, 2919690, 11668616);
	centro_end   <- c(15750321, 3735247, 13674767, 4011692, 12082583);
	centro_start <- centro_start + cumsum(endPoints[1:length(centro_start)]) + (0:4)*offsetBetweenChromosomes;
	centro_end <- centro_end + cumsum(endPoints[1:length(centro_end)]) + (0:4)*offsetBetweenChromosomes;

	outputFileName <- paste("omega_", focalPopulation, ".n_", populationSize, ".grd", gridInterval, ".scwin", range, ".win", windowSize, ".eps", sep="");

	## make the image. 
	setEPS(paper="special", horizontal=TRUE);
	postscript(outputFileName, width=9, height=5, bg="white");

	plot.new();
	plot.window(xlim=c(0, sum(endPoints) + (4*offsetBetweenChromosomes)), ylim=c(baselineOmega, legendPositionScalar*(max(scores1[,"score"]))));
	par(mar=c(5, 4, 1, 1));
	
	axis(2);
	title(ylab=expression(omega));
	title(xlab=paste("Mb, ", focalPopulation, sep=""))

	tickMarks <- c();
	labels <- c();

	###############################################################################
	## plot the 5 separate chromosomes for each population
	for( j in 1:5 ){
		fsub_j <- subset(scores1, chr == j);

		breaks <- seq( 0, endPoints[j + 1], by=4e6);
		chrLabels <- breaks / 1e6;
		chrLabels[seq(2, length(chrLabels), by=2)] <- "";

		labels <- c(labels, chrLabels);
		tickMarks <- c(tickMarks, sum(endPoints[1:j]) + breaks + (j-1)*offsetBetweenChromosomes);

		## x-axis offset!
		fsub_j[,"midpoint"] <- fsub_j[,"midpoint"] + sum(endPoints[1:j]) + (j-1)*offsetBetweenChromosomes;
		points(fsub_j[,"midpoint"], fsub_j[,"score"], col=pointColor, cex=0.2, pch=21, bg=borderColor);
	}

	axis(1, at=tickMarks, labels=labels, cex.axis=0.85);
	box();
#	legend(x=(centro_end[3] + 8e6), y=legendPositionScalar*max(scores1[, "score"]), legend=c(focalPopulation), pch=21, pt.bg=c(focalPopulationColor), pt.cex=1.25, bty="n");
	dev.off();

	cat("We've made a manhattan plot for the global analysis\n");
}