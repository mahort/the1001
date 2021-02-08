# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());

callMethodId <- 2;
snpImputationThreshold <- 0.1;

###############################################################################
## here are the Iberian relicts and non-relicts
###############################################################################
mapFileName <- "~/projects/the_1001/selection/fst/2/imp1/global/results/relicts_nonrelict.only.txt";
relictMap <- read.table(mapFileName, header=T, sep="\t", as.is=T, stringsAsFactors=F);

###############################################################################
## identify non-relicts that aren't already matched to Iberian relicts.
###############################################################################
mapFileName <- "~/projects/the_1001/selection/fst/admix_pop_metadata/metadata_structure.txt";
accessionMap <- read.table(mapFileName, header=T, sep="\t", as.is=T, stringsAsFactors=F);
colnames(accessionMap)[which(colnames(accessionMap) == "tg_ecotypeid")] <- "ecotype_id";
colnames(accessionMap)[which(colnames(accessionMap) == "Admixed")] <- "population_name";

## omit the prev 2 categories.
others <- which(!accessionMap[,"ecotype_id"] %in% relictMap[,"ecotype_id"]) 
accessionMap <- accessionMap[others,];

###############################################################################
## 
###############################################################################
table(accessionMap[,"population_name"]);

outgroupC <- subset(accessionMap, population_name == "Western Europe" & country == "FRA" );
outgroupC <- outgroupC[sample(1:nrow(outgroupC), size=21),];

setwd(paste("~/projects/the_1001/selection/fst/2/imp1/global/results/", sep=""));
write.table(outgroupC, "french_nonrelicts.txt", quote=F, sep="\t", row.names=F); #")
