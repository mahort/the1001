# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());

setwd("~/projects/invasion/public_phenos/korte/");

ft10 <- "ft10_1001.csv";
ft10 <- read.table(ft10, header=T, sep=",", as.is=T, stringsAsFactors=F);

## bjarni's format: 
## phenotype_id, phenotype_name, ecotype_id, value, replicate_id

tmp <- cbind(phenotype_id=rep(1, nrow(ft10)), phenotype_name=rep("ft10", nrow(ft10)), ecotype_id=ft10[,"ecotypeid"], value=ft10[,"mean"], replicate_id=rep(1, nrow(ft10)));

## output this to the public_phenos folder...
write.table(tmp, "ft10_1001.txt", quote=F, sep=",", row.names=F);



# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());

setwd("~/projects/invasion/public_phenos/korte/");
ft16 <- "ft16_1001.csv";
ft16 <- read.table(ft16, header=T, sep=",", as.is=T, stringsAsFactors=F);

## bjarni's format: 
## phenotype_id, phenotype_name, ecotype_id, value, replicate_id

tmp <- cbind(phenotype_id=rep(1, nrow(ft16)), phenotype_name=rep("ft16", nrow(ft16)), ecotype_id=ft16[,"ecotypeid"], value=ft16[,"mean"], replicate_id=rep(1, nrow(ft16)));

## output this to the public_phenos folder...
write.table(tmp, "ft16_1001.txt", quote=F, sep=",", row.names=F);

