# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());

###############################################################################
## omega+ parameters
###############################################################################
callMethodId <- 2;

gridInterval <- 500; ## grid size, not selection scan size.
selectionScanMinWindowStart <- 1e4;
selectionScanMaxWindowStop <- 1e5;

snpImputationThreshold <- 0.1;

{
	###############################################################################
	## code to setup pbs jobs - 3 functions, and (below) the code to setup
	## each (population specific) job.
	###############################################################################
	buildOmegaPlusJob <- function(
			callMethodId,	
			jobName,
			sizeOfSample,
			chromosomeNumber,
			numberOfSlices,
			lengthOfSegment,
			gridIntervalDensity,
			minimumSweepWindowSize,
			maximumSweepWindowSize,
			snpImputationThreshold){
		
		filename <- paste("omega_", jobName, ".", sizeOfSample, ".chr", chromosomeNumber, ".grd", gridIntervalDensity, ".win", minimumSweepWindowSize, "-", maximumSweepWindowSize, ".om.sh", sep="");
		
		qsubString <- paste( 
				"#!/bin/bash\n",
				"#PBS -S /bin/bash\n",
				"#PBS -N op", chromosomeNumber, "\n",
				"#PBS -P the1001genomes\n",
				"#PBS -V\n",
				"#PBS -l walltime=48:00:00\n",
				"#PBS -l select=1:ncpus=1:mem=64gb\n", sep="");
		
		if( numberOfSlices > 1 ){
			qsubString <- paste(qsubString, "#PBS -J 1-", numberOfSlices, "\n", sep="");
		}
		
		qsubString <- paste( qsubString, 
				"SUBPOP=", jobName, "\n",
				"SIZE=", sizeOfSample, "\n",
				"CHR=", chromosomeNumber, "\n", sep="");
		if( numberOfSlices > 1 ){
			qsubString <- paste(qsubString, "JOB=${PBS_ARRAY_INDEX}\n", sep="");
			
		} else {
			qsubString <- paste(qsubString, "JOB=1\n", sep="");
		}
		
		qsubString <- paste(qsubString, 
				"cd ~/projects/the_1001/selection/omega/", callMethodId, "/imp", snpImputationThreshold, "/${SUBPOP}/chr${CHR}/slice${JOB}\n",
				"INPUTFILE=${SUBPOP}_chr${CHR}.n_${SIZE}.slice${JOB}.omega\n", sep="");
		
		qsubString <- paste( qsubString, "OmegaPlus -name ${INPUTFILE}.grd", gridIntervalDensity, ".win", minimumSweepWindowSize, "-", maximumSweepWindowSize, " -input ${INPUTFILE}.txt -grid ", floor(lengthOfSegment/gridIntervalDensity), " -minwin ", minimumSweepWindowSize, " -maxwin ", sprintf("%d", maximumSweepWindowSize), " -length ", lengthOfSegment, " -all\n", sep="");
		
		write.table(qsubString, file=filename, quote=F, row.names=F, col.names=F);
	}

	setwd(paste("~/projects/the_1001/selection/omega/", callMethodId, "/imp", snpImputationThreshold, "/", sep=""));

	chooseFile <- function(pattern){
		
		if( pattern[1] != "" ){
			files <- list.files(path=".", pattern=pattern[1]);
			if( length(pattern) > 1 ){
				for( k in 2:length(pattern)){
					files <- files[grep(pattern[k], files)];
				}
			}
			
		} else {
			files <- list.files(path=".", pattern="sh$");
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
		
		return(list("base_dir"=dirname(getwd()), "filename"=files[fileNumber]));
		
	}; phenotypeFile <- chooseFile(pattern=""); ## need to write the regex to match on method[1] and method[2] ### sometimes, the latter element isn't specified, so condition it 

	filename <- gsub("omega_", "", phenotypeFile$filename);
#	population <- unlist(strsplit(unlist(strsplit(filename, "_"))[2], "\\."))[1];
	population <- unlist(strsplit(filename, "\\."))[1]
	populationSize <- as.numeric(unlist(strsplit(unlist(strsplit(filename, population))[2], "\\."))[2]);
	gridInterval <- as.numeric(unlist(strsplit(unlist(strsplit(filename, "grd"))[2], "\\."))[1]);
	range <- unlist(strsplit(unlist(strsplit(filename, "win"))[2], "\\."))[1];
	range <- gsub("e+", "e\\+", fixed=T, range);

	files <- list.files(path=getwd(), pattern="sh$");
	files <- files[grep(population, files)];
	files <- files[grep(populationSize, files)];
	files <- files[grep(paste(".grd", gridInterval, sep=""), files)];
	files <- files[grep(paste(".win", range, sep=""), files)];
	stopifnot(length(files) == 5);

	for( j in 1:length(files)){
		## determine the length of the section...
		totalLength <- -1;
		cnxn <- file(files[j], open="r"); #header <- unlist(strsplit(readLines(cnxn, n=1, warn=F), '\t'));
		while( length(line <- readLines(cnxn, n=1)) > 0 ){ ## could do this in a batch, but whatever.
			if( grepl("^OmegaPlus -name", line)){
				## got it, so determine the length of the segment.
				fields <- unlist(strsplit(line, " "));
				fieldOfInterest <- grep("-length", fields);
				totalLength <- as.numeric(fields[fieldOfInterest + 1]);
				cat("Determined the length of the region.\n");
				break;
			}

		}; close(cnxn);

		stopifnot(totalLength > 0);

		buildOmegaPlusJob( jobName=population, 
				sizeOfSample=populationSize, 
				chromosomeNumber=j, 
				numberOfSlices=(1), 
				callMethodId=callMethodId, 
				lengthOfSegment=totalLength, 
				gridIntervalDensity=gridInterval, 
				minimumSweepWindowSize=selectionScanMinWindowStart, 
				maximumSweepWindowSize=selectionScanMaxWindowStop,
				snpImputationThreshold=snpImputationThreshold);
		
		cat("Files for chr:", j, "written.\n\n");
	}
}