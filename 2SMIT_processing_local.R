
# Load libraries ----------------------------------------------------------

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(seqinr)
#library(evobiR)


#     !/usr/bin/env Rscript
# Various input files -----------------------------------------------------

# TO DO ITEMS

# 1. change paths
# 2. list intronless control genes used in the list below
intronlessgenes <- c('YAL012W', 'YBR152W', "YOR153W")
# 3. enter sample names to analyze below
analyzeme <- c('hs1A_50', 'hs1B_51', 'hs1C_52', 'hs1D_53')
#analyzeme <- c('hs1A_50')
# 4. set output directory (must end with '/)
outputdir <- c("/disk_aug18_a/jan2019/processed_data100/")
# analyzeme <- c('rrp6')
# set output directory






# creates two subdirectories in output directory specified above called smit_files and smit_curves. Individual samples
# will each create their own subfolders inside of these parent dirs. 
#smit_file_path <- paste(outputdir, 'smit_files', '/', sep='')
#if (dir.exists( smit_file_path ) == FALSE){
#  dir.create( smit_file_path )
#}
#smit_plot_path <- paste(outputdir, 'smit_curves', '/', sep='')
#if (dir.exists( smit_plot_path ) == FALSE){
#  dir.create( smit_plot_path )
#}

terminalExonLengths <- read.table('/home/common/SMITs/R/SMIT_accessory_files/Scer3_snyder_TElength.txt', sep='\t')
# source('/home/common/SMITs/R/SMIT_accessory_files/expCdfFunction.R')
source('/home/common/SMITs/R/SMIT_accessory_files/smitData.R')
source('/home/common/SMITs/R/SMIT_accessory_files/commonFunctions.R')



# input data
# replicate 1 data
hs1A_50 <- fread("/disk_aug18_a/jan2019/neug_scripts/hs1A_50_sorted.bed" , sep = '\t')
#hs1B_51 <- fread("/disk_aug18_a/jan2019/neug_scripts/hs1B_51_sorted.bed" , sep = '\t')
#hs1C_52 <- fread("/disk_aug18_a/jan2019/neug_scripts/hs1C_52_sorted.bed" , sep = '\t')
#hs1D_53 <- fread("/disk_aug18_a/jan2019/neug_scripts/hs1D_53_sorted.bed" , sep = '\t')

names(hs1A_50) <- c("chr", "start", "end", "id", "score", "strand", "thickStart", "thickEnd", "color", "block_number", "block_size", "block_start")
#names(hs1B_51) <- c("chr", "start", "end", "id", "score", "strand", "thickStart", "thickEnd", "color", "block_number", "block_size", "block_start")
#names(hs1C_52) <- c("chr", "start", "end", "id", "score", "strand", "thickStart", "thickEnd", "color", "block_number", "block_size", "block_start")
#names(hs1D_53) <- c("chr", "start", "end", "id", "score", "strand", "thickStart", "thickEnd", "color", "block_number", "block_size", "block_start")
# names(rrp6) <- c("chr", "start", "end", "id", "score", "strand", "thickStart", "thickEnd", "color", "block_number", "block_size", "block_start")


# a vector of the genes that you used primers for. Important to include because you will likely get
# data mapping to genes that were not in your primer pool.
genes_primed <- read.csv("/disk_aug18_a/jan2019/62.list" , sep = '\t' , header = T, na.strings = "")
genes_primed <- as.character(genes_primed$Systematic_name)
genes_primed <- genes_primed[!is.na(genes_primed)]


# yeast mRNA annotation in bed format
scer3_bed_new <- read.table("/home/common/SMITs/R/SMIT_accessory_files/snyder_all.txt")
# scer3_bed_new <- read.table("C:/Users/taraa/Documents/_lab/annotations/snyder/snyder_all.txt")
scer3_bed_new$score <- rep(0, length(scer3_bed_new$V1))
colnames(scer3_bed_new) <- c('chr', 'start', 'end', 'strand', 'name', 'score')
scer3_bed_new <- scer3_bed_new[c('chr', 'start', 'end', 'name', 'strand', 'score')]
scer3_bed <- scer3_bed_new



# scer3_bed above lists the start coordinate as the leftmost coordinate regardless of strand. 
# this data frame below flips the start and stop coordinates for those genes on the Crick strand.
scer3adjustedforstrand <- scer3_bed
scer3adjustedforstrand$strand <- as.character(scer3adjustedforstrand$strand)
for (i in 1:length(scer3adjustedforstrand$start)){
  if (scer3adjustedforstrand$strand[i] == '-'){
    start <- scer3adjustedforstrand$start[i]
    scer3adjustedforstrand$start[i] <- scer3adjustedforstrand$end[i]
    scer3adjustedforstrand$end[i] <- start
  } 
}


# filter the annotation to just include the SMIT genes without intronless controls
genes_nointronless <- genes_primed[ !( genes_primed %in% intronlessgenes ) ]
genes <- scer3adjustedforstrand[ scer3adjustedforstrand$name %in% genes_nointronless , ]
genes$chr <- as.character( genes$chr )


# filter the annotation to just include the SMIT genes with intronless controls
genes_withintronless <- scer3adjustedforstrand[ scer3adjustedforstrand$name %in% genes_primed , ]
genes_withintronless$chr <- as.character( genes_withintronless$chr )


# intron annotation file
introns <- read.table("/home/common/SMITs/R/SMIT_accessory_files/Gene_FeatureType_Introns.tsv", header = F)
#introns <- read.table("C:/Users/taraa/Documents/_lab/_machinelearning/data/2018 First Half of ML data 1345/Gene_FeatureType_Introns.tsv", header=F)
names(introns) <- c("name", "fiveSS", "threeSS")
introns$intron_size <- introns$threeSS - introns$fiveSS
introns$name <- as.character(introns$name)
introns <- introns[-grep("^Q", introns$name),]
for (i in 1:nrow(introns)){
  strand <- strsplit(introns$name[i], split = "")[[1]][7]
  # print(paste(i,': ', strand))
  #if (strand == 'C'){
  if (strand == '-'){
    fivess <- introns$fiveSS[i]
    threess <- introns$threeSS[i]
    introns$fiveSS[i] <- threess
    introns$threeSS[i] <- fivess
  } 
}


outputdir <- paste(outputdir, 'smit_files/', sep='')
outputdir_plots <- gsub( 'smit_files/' , 'smit_curves/', outputdir)



# Functions ---------------------------------------------------------------

mean2 <- function(x){
  mean(x,na.rm=TRUE)
}
sd2 <- function(x){
  sd(x, na.rm = TRUE)
}


# finds R1 and R2 reads and merges them into one line in BED file
find.mate <- function(data){
  temp <- sub("/[12]", "", data$id)
  temp2 <- sub(".*/", "", data$id)
  data$id <- temp
  data$pair <- temp2
  data_one <- data[data$pair == 1,]
  data_two <- data[data$pair == 2,]
  test <- data.frame()
  test <- dplyr::inner_join(data_one, data_two, by="id")
  # return(test)
  assign( paste(analyzeme[i], '_paired', sep ='') , test, envir = .GlobalEnv )
}



# calculates the insert size for each paired read. Insert size is the length of fragment between the R1 and R2 reads 
calc.insert.size <- function(data){
  
  # Pair 2 (.y) spans the intron and Pair 1 (.x) is the 3' end
  # For watson strand genes, end.x - start.y should be greater
  # For crick strand genes, end.y - start.x should be greater
  
  data$strand_gene <- substring(data$name, 7, 7)
  data$block_number.y <- as.integer(data$block_number.y)
  data$block_start.y <- as.character(data$block_start.y)
  data$block_size.y <- as.character(data$block_size.y)
  
  spliced <- data[ data$block_number.y == 2, ]
  unspliced <- data[ data$block_number.y == 1, ]
  spliced$intron_size <- as.integer(unlist(strsplit(spliced$block_start.y, ','))[2*(1:length(spliced$block_start.y))]) -
    as.integer(unlist(strsplit(spliced$block_size.y, ','))[2*(1:length(spliced$block_size.y))-1])
  unspliced$intron_size <- rep(NA)
  
  #spliced$insert_size <- ifelse( spliced$strand_gene == 'W',
  spliced$insert_size <- ifelse( spliced$strand_gene == '+',
                                 (spliced$end.x - spliced$start.y) - spliced$intron_size ,
                                 (spliced$end.y - spliced$start.x) - spliced$intron_size )
  unspliced$insert_size <- ifelse( unspliced$strand_gene == '+' ,
                                   unspliced$end.x - unspliced$start.y,
                                   unspliced$end.y - unspliced$start.x)
  
  final <- bind_rows(spliced, unspliced)
  
  assign( paste(analyzeme[i], '_insertsize', sep='') , final, envir = .GlobalEnv )
  # return(final)
}


# add name of gene to the bed file
add.name <- function(data){
  for ( j in 1:length(genes_withintronless$name )){
    tttt <- length(genes_withintronless$name)
	print(tttt)
    print(j)
    name <- genes_withintronless$name[j]
    strand <- substr(name, 7, 7)
	print( name)
	print( strand)
    if (strand == '+'){
      reads <- data[ (data$chr.y == genes_withintronless$chr[j] & data$start.y > genes_withintronless$start[j] - 100 &
                        data$start.y < genes_withintronless$end[j]) , ]
      reads$name <- rep( name, length(reads$chr.x) )
    } else if (strand == '-'){
      reads <- data[ data$chr.y == genes_withintronless$chr[j] & data$end.y < genes_withintronless$start[j] + 100 &
                       data$end.y > genes_withintronless$end[j] , ]
      reads$name <- rep( name, length(reads$chr.x) )
    }
    
    if ( j == 1 ){
      final <- reads
    } else{
      final <- bind_rows( final, reads )
    }
  }
  # return(final)
  assign( paste(analyzeme[i], '_named', sep ='') , final, envir = .GlobalEnv )
}




# reformat data structure to make downstream analysis easier
reformat.2 <- function(data){
  for (h in 1:length(genes_nointronless)){
    print(paste("Processing", genes_nointronless[h]))
    name <- genes_nointronless[h]
    threeSS <- introns[introns$name == name,3]
    fiveSS <- introns[introns$name == name,2]
    # start <- scer3_bed[scer3_bed$name == name,2]
    start <- scer3adjustedforstrand[scer3adjustedforstrand$name == name, 2]
    # end <- scer3_bed[scer3_bed$name == name,3]
    end <- scer3adjustedforstrand[scer3adjustedforstrand$name == name, 3]
    intron_size <- introns[introns$name == name,4]
    
    a <- data[data$name == name,]
    a <- a[a$chr.x == a$chr.y,]
    
    if (length(a$end.x) > 2){
      if (a$strand.x[3] == '+'){   # Crick Strand
        df.smit <- data.frame(position = seq(from=(threeSS + 100), to=(end - 500)))
        df.smit$relPosition <- (threeSS - df.smit$position)
        df.smit$splicedLength <- abs(df.smit$position - start) - intron_size
        df.smit$unsplicedLength <- abs(df.smit$position - start)
        a_s <- a[a$block_number.y == 2,]
        a_un <- a[a$block_number.y == 1,]
        splicedReads <- data.frame(table(a_s$start.x))
        unsplicedReads <- data.frame(table(a_un$start.x))
        
        
      } else if (a$strand.x[3] == '-') {   # Watson Strand
        df.smit <- data.frame(position = seq(from=(threeSS - 100), to=(end + 500)))
        df.smit$relPosition <- (df.smit$position - threeSS)
        df.smit$splicedLength <- (df.smit$position - start) - intron_size
        df.smit$unsplicedLength <- df.smit$position - start
        a_s <- a[a$block_number.y == 2,]
        a_un <- a[a$block_number.y == 1,]
        splicedReads <- data.frame(table(a_s$end.x))
        unsplicedReads <- data.frame(table(a_un$end.x))
        
      }
      names(splicedReads) <- c('position', 'freq')
      names(unsplicedReads) <- c('position', 'freq')
      splicedReads$position <- as.numeric(as.character(splicedReads$position))
      unsplicedReads$position <- as.numeric(as.character(unsplicedReads$position))
      
      t <- left_join(df.smit, splicedReads, by = 'position')
      df.smit <- left_join(t, unsplicedReads, by = 'position')
      names(df.smit) <- c('position', 'relPosition', 'splicedLength', 'unsplicedLength', 
                          'splicedReads', 'unsplicedReads')
      df.smit[is.na(df.smit)] <- 0
      write.table(df.smit, paste(name, ".smit", sep=""), row.names = F, quote=F, sep='\t')
      rm("df.smit")
    }
  }
}



smitAnalysis.process  <- function( file ) {
  # Read raw data 
  print(paste("Processing",file))
  smitData  <- smit.read( file )
  insertLength.para  <- read.table("insertLengthBias.txt",sep='\t',head=T)
  
  # Remove outliers 
  smitData  <- smit.defineXRange(smitData, minX=-2000, maxX=5000 )
  if( length(smit.getRelPosition(smitData)) != 0 ) {
    # Bin data 
    #smitData  <- smit.bin(smitData,binSize=5)
    
    # Calculate splicing values 
    smitData  <- smit.setSplicingValues(smitData, insertLengthCutOff=insertLength.para$cutoff, insertLengthDecayRate=-insertLength.para$rate)
    
    # Write splicing values to file 
    smit.writeSplicingValues( smitData )
    #smit.plot(smitData, weight=1)
    # Shift relative position to the exit channel of Pol II ~36nt downstream of the 3'SS
    #smitData  <- smit.shiftRelPosition( smitData, -36 )
    
    # Fit exponential to the data
    #smitData  <- smit.fitExpCDF( smitData, 1 )  
  } else {
    print(paste("No valid positions for gene", file) )
  }
}


# Define exponential cdf function 
expCDF  <- function(x,k,delay=0) {
  a  <- x + delay
  return( pexp(a,rate=k) )
}

# Fit exp cdf to data
expCDF.fit  <- function( x, y, weights=rep(1,length( y ) ), kStart=0.005 ) {
  out <- tryCatch(
    out <- nls( y ~ ( expCDF( x, k) ), start=list(k=kStart),weights=weights,control = nls.control(maxiter=5000),algorithm="port",lower=list(k=0),trace=F),
    error=function( cond ) {
      message( cond )
      return( NA )
    },
    warning=function( cond ) {
      message( cond )
      return( NA )
    },
    finally={
    }
  )    
  return(out)
}

# Predict the exp cdf values for a rate
expCDF.predict <- function( x, k ) {
  return( expCDF( x, k ) )
}



# This the binned plotting method that can be seen in the 2016 Carrillo, Herzel et al Cell publication
# I prefer the sliding window based plotting method included in the for loop immedately below the
# "Execute" section below

# plot.splicingValues_new <- function( statistics, id, model=NA, write=T, pathOut="", xrange=NA ) {
#   
#   exonEnd  <- getTerminalExonEnd( id )
#   saturation  <- getSaturation( statistics, id )
#   
#   splicingValuesReached  <- data.frame( getCrossingPoint( statistics, id, 0.5 ) ) 
#   
#   df <- statistics
#   df$label  <- rep("data",nrow(df) ) 
#   
#   if( !is.na( model ) ) {
#     xPredict <- seq(from=0,to=max( df$meanPos, na.rm=T ), by=1 )
#     yPredict <- predict( object=model, newdata=data.frame( x=xPredict ) )
#     dfPredict <- data.frame( meanPos=xPredict, meanY=yPredict, sdY=rep( 0, length( xPredict ) ), pValue=rep(NA,length( xPredict ) ), reads=rep(NA,length( xPredict ) ),  memberCount=rep(NA,length( xPredict ) ), label=rep("predict",length( xPredict ) ) )
#     
#     df <- rbind( df, dfPredict )
#   }
#   
#   limits  <- aes( ymax=meanY+sdY,ymin=meanY-sdY)
#   title <- paste( id, "splicingValue" )
#   p  <- ggplot( df, aes(x=meanPos,y=meanY ) ) # Read data 
#   p  <- p + ylab(label="fraction co-txn splicing") + ylim(0,1) # Label the y axis
#   p  <- p + xlab(label="pos (rel to 3'SS)")  + expand_limits(x=0) # Label the y axis
#   p  <- p + ggtitle( title )
#   
#   # Plot data layer
#   p  <- p + geom_point(data=df[ df$label == "data", ], color="black", size = 4) 
#   p  <- p + geom_errorbar(data=df[ df$label == "data", ], limits, width=10, size = 1.5)
#   
#   # Plot model layer 
#   p  <- p + geom_line(data=df[ df$label == "predict", ], color="black" )
#   
#   p <- p + geom_vline( xintercept = saturation[ 1 ], colour="black", linetype = "dashed", size =2 )
#   p <- p + geom_hline( yintercept = saturation[ 2 ], colour="#399088", linetype = "dashed", size = 2 )
#   p <- p + geom_vline( xintercept = splicingValuesReached$x , colour="#aa859c", linetype = "dashed", size = 2 )
#   if( exonEnd <= max( statistics$meanPos, na.rm=T ) ) {
#     p <- p + geom_vline( xintercept = exonEnd, colour="#868887", linetype = "dashed", size =2 )
#   }
#   
#   if( !is.na(xrange) ) {
#     p <- p + expand_limits(x = xrange)
#   } 
#   
#   if( write ) {
#     fileOut  <- paste( c( pathOut, id, "_fit.pdf" ), collapse="")
#     pdf(file = fileOut, width = 7, height = 7 , useDingbats = F )
#     p
#     dev.off()
#     # ggsave( filename=fileOut,plot=p, width=7, height=7, useDingbats=FALSE )
#   }
#   
#   return( p )
# }




# Execute -----------------------------------------------------------------

for (i in 1:length( analyzeme )){
  
  # find mated pair
  find.mate( get(analyzeme[i]) )
  print( 'Mates found' )
  
  # add name of gene to each line of Bed file
  add.name( get(paste(analyzeme[i], '_paired', sep='')) )
  print( 'Completed gene identification.')
  
  # calculate insert size
  calc.insert.size( get(paste(analyzeme[i], '_named', sep='')) )
  print( 'Completed insert size calculation.' )
  
  
  # check if a sample-named directory already exists and if not, create it and move into it
  path <- paste(outputdir, analyzeme[i], '/', sep = '')
  if (dir.exists( path ) == FALSE){
    dir.create( path )
  }
  setwd( path )

  # reformat the files to cooperate with fernando's code
  reformat.2( get(paste(analyzeme[i], '_insertsize', sep='')) )
  print( 'Completed reformating.')
  
  # rm(list = c( paste(analyzeme[i])))
  # rm(list = c( paste(analyzeme[i], '_insertsize', sep='')))
  
  df <- get( paste(analyzeme[i], '_insertsize', sep ='') )
  # assign( paste('intronless_', analyzeme[i] , sep=''), df[ df$name %in% intronlessgenes , ] )
  # assign( paste('intron_containing_', analyzeme[i], sep=''), df[ !(df$name %in% intronlessgenes ) , ] )
  intronless <- df[ df$name %in% intronlessgenes , ]
  intron_containing <- df[ !( df$name %in% intronlessgenes ) , ]
  
  
  ## this is Fernando's code to normalize to intronless genes
  a <- data.frame( insertLengths = intronless$insert_size)
  a$bin <- cut(x=a[,1], breaks=seq(from=0, to=max(a[,1]), by=10), labels=seq(from=5, to=max(a[,1])-5, by=10 ) )
  data.frame( table( a$bin ) ) -> b #Get counts over each position 
  xLow <- 100
  xHigh <- 500
  b[,1] <- as.numeric( as.character( b[,1] ) ) #Make levels numeric
  plot( b[,1], log(b[,2]), cex=0.1, xlim=c(0,1000), ylab="counts (Log e)",xlab="insert size (nt)") # Log-linear plot
  arrows( xLow,0,xLow,40, col='red' )
  arrows( xHigh,0,xHigh,40, col='red' )
  # Extract data within boundaries
  b[b[,1]>xLow & b[,1]<=xHigh,1] -> x
  b[b[,1]>xLow & b[,1]<=xHigh,2] -> y 
  fit <- lm( log(y+1) ~x ) # I added the +1 after y that was not there in fernando's code because I had infinite values when calculating log(y)
  abline( fit, col='blue')
  
  write.table(data.frame(cutoff=xLow, rate=fit$coefficients[2]), 
              file = "insertLengthBias.txt",
              quote = F, sep = "\t", row.names = F) 
  insertLength.para <- read.table("insertLengthBias.txt", header = T)
  
  rm(list = c( paste(analyzeme[i])))
  rm(list = c( paste(analyzeme[i], '_insertsize', sep='')))
  rm(list = c( paste(analyzeme[i], '_named', sep='')))
  rm(list = c( paste(analyzeme[i], '_paired', sep='')))
  
  # begin plotting
  files <- list.files('.' , pattern = '*.smit' )
  files <- as.character( files )
  
  # this calls the functions from the source files loaded in the beginning that normalize the data
  # to the intronless controls (code by Fernando Carrillo Oesterreich)
  sapply( files, smitAnalysis.process )
  
  pathOut <- paste( outputdir_plots, analyzeme[i], '/', sep='' )
  if (dir.exists( pathOut ) == FALSE){
    dir.create( pathOut )
  }
  
  files  <- list.files(".", pattern="splicingValues.txt$")
  summaryStatistics <- ldply(lapply(files, summary_statistics))
  
  
  # Plot SMIT curves using a sliding window
  for( k in 1:length( files ) ){
    file  <- files[ k ]
    print( paste( as.character(k), file  ) )
    rcs  <- read_count_summary(file) 
    #print( rcs  )
    id  <- strsplit( file, split="_")[[1]][ 1 ]
    # Consider only transcripts with more or equal valid positions than bins
    binCount <- 30
    if( rcs$position_count >= binCount ) {
      data <- read.table(file, header = T)
      data <- data[data$relPos >= -30,]
      values <- SlidingWindow(mean2, data[,4],30,1)
      stds <- SlidingWindow(sd2, data[,4], 30, 1)
      limits <- aes(ymax = (values + stds)[1:range], ymin=(values - stds)[1:range])
      min <- min(data[,1])+30
      #max <- max(data[,1])
      end <- scer3adjustedforstrand[scer3adjustedforstrand$name == id,3]
      relend <- abs(end - introns[introns$name == id,3])
      range <- relend+100 - min + 1
      df <- data.frame(x = seq(min, relend+100), y = values[1:range], stds = stds[1:range])
      xend <- max(df[!is.na(df$y),1])
      
      counter <- 0
      for (k in 1:length(df$x)){
        if (is.na(df$y[k])){
          counter <- counter + 1
          xend <- df$x[k] - 10
        } else { counter <- 0 }
        if (counter >= 10){ break } 
      }
      
      if (xend > relend){  satend <- relend  } else { satend <- xend }
      y <- df$y[df$x %in% seq((satend-119), satend)]
      sat <- mean(y, na.rm = T)
      
      p <- ggplot(df)+
        geom_line(aes(x,y), size = 2)+
        geom_ribbon(aes(x =x , ymin = (y-stds), ymax = (y+stds)), alpha = 0.2)+
        geom_vline(aes(xintercept = relend), colour ='blue', linetype='longdash')+
        geom_hline(aes(yintercept = sat), colour = 'red', linetype = 'longdash')+
        ggtitle(paste('Splicing Kinetics for', id))+xlab("Position relative to 3' SS")+
        ylab("Fraction spliced co-transcriptionally")+
        coord_cartesian(xlim=c(0,xend), ylim=c(0,1))
      
      pdf(paste(pathOut, id, '.pdf', sep=''), useDingbats = F, width = 10)
      show(p)  
      dev.off()
    }
  
  }
} 



