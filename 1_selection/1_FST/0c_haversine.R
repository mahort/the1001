# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());
require(rhdf5);

callMethodId <- 2;
numberOfReads <- 1e4;

###############################################################################
# core function: weir & cockerham (1984) theta. 
# this code expects homozygous AA to be coded with the number: 2. Heterozygote Aa: 1. Homozygote aa: 0.
###############################################################################
Fst <- function(n, ns, ps, hs){ 
	
	nbar <- mean(ns);
	pbar <- sum(ns*ps)/(nbar*n);
	
	nc <- (n*nbar - sum(ns^2)/(n*nbar))/(n-1);
	
	s2 <- sum(ns*((ps - pbar)^2))/((n-1)*nbar);
	hbar <- sum(ns*hs)/(n*nbar);
	
	# from pages 1359, 1360 in Weir & Cockerham (1984)
	a <- (nbar/nc)*(s2 - (1/(nbar-1))*(pbar*(1-pbar) - ((n-1)/n)*s2 - (1/4)*hbar));
	b <- (nbar/(nbar-1))*((pbar*(1-pbar)) - ((n-1)/n)*s2 - ((2*nbar-1)/(4*nbar))*hbar);
	c <- (1/2)*hbar;
	
	return(a/(a+b+c));
}

###############################################################################
# to convert matrix to TWOS format, just multiply it by 2, or, when using a file:
#  sed -ie 's/,1\b/,2/g' tmp.csv
###############################################################################

setwd(paste("/net/gmi.oeaw.ac.at/nordborg/lab/Data/genotype-callmethods/PYGWAS_GENOTYPES/", callMethodId, "/", sep=""));
filename <- "all_chromosomes_binary_uncompressed.hdf5";
accessions <- as.numeric(h5read(filename, "accessions"));

#attributes <- h5ls(filename, all=T);
positions <- h5read(filename, name="positions", read.attributes=T);
chromosomes <- attributes(positions)$chr_regions;
chromosomes[1,] <- chromosomes[1,] + 1;
chr.ids <- attributes(positions)$chrs;
snpsPerChromosome <- apply(chromosomes, 2, function(x){ x[2] - x[1]; }) + 1;

###############################################################################
## distinguish between the lines we want, and those we don't...
###############################################################################
#metadata <- read.table("~/projects/the_accessions/the1001.genomes.txt", header=T, sep="\t", as.is=T, stringsAsFactors=F);
#subpops <- c("Fennoscandia", "Iberia", "Central-Europe", "Americas", "Eastern-Range", "Appenine", "British-Isles", "Southeast-Europe", "France", "Benelux");
metadata <- "~/projects/the_1001/selection/fst/admix_pop_metadata/metadata_structure.txt";
metadata <- read.table(metadata, header=T, sep="\t", as.is=T, stringsAsFactors=F);
colnames(metadata)[which(colnames(metadata) == "Admixed")] <- "population_name";
colnames(metadata)[which(colnames(metadata) == "tg_ecotypeid")] <- "ecotype_id";
metadata <- subset(metadata, (country == "ESP" | country == "POR"));

##############################################################################
#  determine the lat/long for individual lines, and the distance from each line...
###############################################################################
metadata[,"latitude"] <- metadata[,"latitude"] * (pi/180);
metadata[,"longitude"] <- metadata[,"longitude"] * (pi/180);

haversineDistance <- function(lat1, long1, lat2, long2, radiusOfTheEarth=6371){
	
	deltaLongitude <- abs(long2 - long1);
	
	y <- sqrt((cos(lat2) * sin(deltaLongitude))^2 + ((cos(lat1) * sin(lat2)) - (sin(lat1) * cos(lat2) * cos(deltaLongitude)))^2);
	x <- (sin(lat1) * sin(lat2)) + (cos(lat1) * cos(lat2) * cos(deltaLongitude));

	deltaSig <- atan2(y, x);	
	d <- radiusOfTheEarth * deltaSig;
	return(d);
}


## for each relict, identify a matching non-relict.
relicts <- subset(metadata, population_name == "Relicts");
nonrelicts <- subset(metadata, population_name != "Relicts");

pairs <- matrix(nrow=0, ncol=3);
for( j in 1:nrow(relicts)){
	relict_j <- relicts[j, c("latitude", "longitude")];

	for( k in 1:nrow(nonrelicts)){
		nonrelict_k <- nonrelicts[k, c("latitude", "longitude")]; #"
		haversine_jk <- haversineDistance(long1=relict_j$longitude, lat1=relict_j$latitude, long2=nonrelict_k$longitude, lat2=nonrelict_k$latitude);
		pairs <- rbind(pairs, c(relicts[j, "ecotype_id"], nonrelicts[k, "ecotype_id"], haversine_jk));
	}
}

colnames(pairs) <- c("relict", "non_relict", "distance_in_km");
write.table(pairs, paste("~/projects/the_1001/selection/fst/", callMethodId, "/results/matches_relicts_nonrelicts.txt", sep=""), quote=F, row.names=F, sep="\t");
pairs <- read.table(paste("~/projects/the_1001/selection/fst/", callMethodId, "/results/matches_relicts_nonrelicts.txt", sep=""), header=T, sep="\t", as.is=T, stringsAsFactors=F);
#minimums <- stack(tapply(pairs[,"distance_in_km"], pairs[,"relict"], min));
#minimums[,2] <- as.numeric(as.character(minimums[,2]));
#colnames(minimums) <- c("distance", "ecotype")
## determine the minimum for each id.
pairs <- pairs[order(pairs[,"distance_in_km"]),];
relictIds <- unique(pairs[,"relict"]);

matches <- matrix(nrow=0, ncol=3);
for( j in 1:length(relictIds)){
	relict_j <- relictIds[j];
	subset_j <- subset(pairs, relict == relict_j);

	## who is the closest line? have we already used it?
	for( k in 1:nrow(subset_j)){
		match_k <- subset_j[k, "non_relict"];
		if( !match_k %in% matches[,2]){
			matches <- rbind(matches, c(relict_j, match_k, subset_j[k, "distance_in_km"]));
			break;
		}
	}
}

colnames(matches) <- c("relict", "non_relict", "distance_in_km");
write.table(matches, paste("~/projects/the_1001/selection/fst/", callMethodId, "/results/matches_relicts_nonrelicts_paired.txt", sep=""), quote=F, row.names=F, sep="\t");
pairs <- read.table(paste("~/projects/the_1001/selection/fst/", callMethodId, "/global/imp1/results/matches_relicts_nonrelicts_paired.txt", sep=""), header=T, sep="\t", as.is=T, stringsAsFactors=F);

pairs <- subset(pairs, distance_in_km <= 300);
relicts <- subset(metadata, ecotype_id %in% pairs[,"relict"]);
nonrelicts <- subset(metadata, ecotype_id %in% pairs[,"non_relict"]);
nonrelicts[,"population_name"] <- "nonrelict_contrast";
metadata <- rbind(relicts, nonrelicts);

subpops <- unique(metadata[,"population_name"]);

populations <- list();
populations$n <- 2;
ns <- numeric(populations$n);

for( j in 1:ncol(chromosomes)){
	cat("Retrieving snps for chromosome", j, "\n");
	indices.chr_j <- chromosomes[,j];

	breaks <- floor(seq(indices.chr_j[1] - 1, indices.chr_j[2], len=ceiling(snpsPerChromosome[j]/numberOfReads)));

	results <- list();
	for( k in 1:(length(breaks) - 1 )){
		cat("On slice:", k, "of", length(breaks), "\n");
		sliceStart <- breaks[k] + 1; 
		sliceStop <- breaks[k + 1];

		pos.slice_k <- as.numeric(positions[sliceStart:sliceStop]); ## get the start/stop positions for the SNPs.
		ps.slice_k <- matrix(nrow=length(pos.slice_k), ncol=populations$n);
		hs.slice_k <- matrix(nrow=length(pos.slice_k), ncol=populations$n);

		for( l in 1:length(subpops)){

			subset_k <- subpops[l];
			ecotypeIds <- metadata[which(metadata[,"population_name"] == subset_k), "ecotype_id"];
			indices <- which(accessions %in% ecotypeIds);

			snps_k <- t(h5read(filename, name="snps", index=list(indices, sliceStart:sliceStop))) * 2;

			ns[l] <- length(ecotypeIds);
			ps.slice_k[,l] <- rowMeans(snps_k)/2;
			hs.slice_k[,l] <- apply(snps_k, 1, function(x){ sum(x == 1)})/ns[l];
			cat(".");
		}

		fsts <- numeric(length(pos.slice_k));
		for( m in 1:length(fsts)){
			fsts[m] <- Fst(populations$n, ns=ns, ps=ps.slice_k[m,], hs=hs.slice_k[m,]);
		}

		results[[k]] <- cbind(j, pos.slice_k, round(fsts, 5));
		cat("\n");
	}

	results <- do.call(rbind, results);
	colnames(results) <- c("chr", "pos", "fst");
	write.table(results, paste("~/projects/the_1001/selection/fst/", callMethodId, "/results/relict_nonrelict.fst.chr", j, ".txt", sep=""), quote=F, sep=",", row.names=F);
}
