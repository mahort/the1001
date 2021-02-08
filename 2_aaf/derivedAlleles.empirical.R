# take 
# 1. all snps
# 2. ancestral allele frequency
# 3. ecotype country map.

# make Derived Allele Frequency spectra per populations:
# Sweden/UK/USA/France-Spain-Portugal/Germany/Eastern Europe-Asia

map <- read.table("~/agdp/newsnps/accession_locale.txt",header=T,sep="\t",stringsAsFactors=FALSE);
lyrata <- read.table("~/agdp/newsnps/alleles_w_anc_021910",header=T,sep="\t",stringsAsFactors=FALSE);
snps.all <- read.table("~/agdp/newsnps/snps.derived.letters.txt",header=T,sep="\t",stringsAsFactors=FALSE);

populations <- c("Americas","Iberia","UK-Isles","France","Central-Europe","South-Central","Scandinavia","Austria-Hungary","Eastern-Range");
all.aafs <- list();

sample <- 1:8
for( i in min(sample):max(sample)){

	pop_i <- populations[i];
	ids <- subset(map,pop==pop_i)[,"ecotype_id"];

	cat("Working on population:",pop_i,"\n");

	snps.pop <- snps.all[,paste("X",ids,sep="")];
	snps.pop <- cbind(1:nrow(snps.pop),snps.pop);
	colnames(snps.pop)[1] <- "num";
	
	# determine AAF for this subset.
	pop.aafs <- apply(snps.pop,1,function(x){ 
				tab <- table(t(x[-1]));
				index <- which(names(tab)==lyrata[as.numeric(x[1]),"lyrata_consensus"]);
				if(length(index)==0){
					NA;

				} else {
					tab[index];
				}
	});

	n <- length(ids);
	pop.aafs <- cbind(rep(pop_i,nrow(snps.pop)),pop.aafs/n);

	fileName <- paste(pop_i,".snps.empirical.n",n,".Robj",sep="");
	save(pop.aafs,file=fileName);
	write.table(pop.aafs,gsub(".Robj",".txt",fileName),quote=F,sep="\t",row.names=F,col.names=F);
}

quit(save="no");

file <- "europe.empirical.txt"; # these are resampled w/o replacement now.

draws <- read.table(file,header=F,sep="\t",as.is=T,stringsAsFactors=FALSE);
colnames(draws) <- c("pop",paste("emp",1,sep=""));

mean.aafs <- unstack(draws, emp1 ~ pop);
populations <- c("Americas","Iberia","UK.Isles","France","Central.Europe","South.Central","Scandinavia","Austria.Hungary","Eastern.Range");



pop.colors <- c("dodgerblue3","mediumpurple","firebrick1","darkolivegreen4","hotpink3","cadetblue2","darkslateblue","darksalmon","mediumseagreen");
# nb: check your table order to make sure that the colors matchup... namely:

aaf.tables <- c();
for( i in 1:length(populations)){
	aaf.tables[[populations[i]]] <- table(cut(mean.aafs[,populations[i]],breaks=seq(0,1,by=0.05),include.lowest=T));
}

aaf.tables <- do.call(rbind,aaf.tables);
aaf.tables <- t(aaf.tables);
aaf.tables <- scale(aaf.tables,center=F,scale=colSums(aaf.tables));

# do the linear model and add the slopes to our plot (below).
range <- seq(0,1,by=0.05);

for( i in 1:5 ){
	start <- range[i];
	end <- 1 - start;

	start_ind <- which(range==start);
	end_ind <- length(range) - start_ind;

	cat("start:",start," on index:",start_ind,"\n");
	cat("end:",end," on index:",end_ind,"\n");

	tables_subset <- aaf.tables[start_ind:end_ind,];

	colors <- matrix(nrow=length(pop.colors),ncol=nrow(tables_subset));
	colors[,c(1:ncol(colors))] <- pop.colors
	colors <- t(colors);
	colors <- as.vector(colors);

	mids <- seq((start+0.025),(end-0.025),by=0.05);
	lm.aafs <- apply(tables_subset,2,function(x){ lm(x ~ mids); });

	slopes <- c();
	for( i in 1:length(lm.aafs)){
		slopes <- append(slopes,lm.aafs[[i]]$coefficients[2]);
		names(slopes)[i] <- names(lm.aafs[i]);
	}

	slopes <- paste("(",base::format(slopes,sci=F,digits=2),")",sep="");

	# per pop
	postscript(paste("pop.aaf.hap.resample.bypop.n17.maf",start,".eps",sep=""),width=7,height=7);
	bp <- barplot(tables_subset,beside=T,xaxt="n",xlab=paste("Ancestral allele frequencies (",start,"-",end,"%)",sep=""),col=colors,ylim=c(0,max(tables_subset)+0.015));
	box();
	axis(side=1,at=c(bp[floor(nrow(tables_subset)/2),]),labels=c("US","IB","UK","FR","CE","SC","FS","AH","ER"));
	mtext(slopes,side=1,at=c(bp[floor(nrow(tables_subset)/2),]),line=1.8,cex=0.7);
	dev.off();
}
