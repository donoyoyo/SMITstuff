##########################
# Common functions for methods in non-parameter analysis. 
########################## 

getPValue  <- function( v1, v2 ) {
  out <- tryCatch(
{
  t.test( v1, v2, alternative="less" )$p.value # Try to calculate the t.test and retrieve the pValue
},
error=function(cond) { # catch an error and return NA 
  return(NA)
},
warning=function(cond) { # catch a warning and return NA   
  return(NA)
},
finally={ } )    
  return(out)
}

binAndSummarize <- function( file, id, remove_positions_after_exonEnd=F, binCount=25 ) {
  
  # Read input
  df  <- read.table(file=file,sep="\t",head=T)
  
  # Remove positions after exon if specified 
  if(remove_positions_after_exonEnd) {
    end <- getTerminalExonEnd(id)
    if( length(end ) > 0 ) { df <- df[ df$relPos <= (getTerminalExonEnd(id) - 5), ]  }
  }
  
  # Filter read counts. Remove all values upstream of the 3'SS 
  minReadCount  <- 10
  df <- df[ df$readSum >= minReadCount, ]
  df  <- df[ df$relPos >= 0, ]
  
  # Bin with constant number of bins. This means that the number of positions in bins varies in response to length. 
  #binCount  <- 25
  df$bin  <- cut( df$relPos, breaks=binCount )
  
  # Get summary statistic of specified bin: meanPos, meanY, sdY and pValue to next bin 
  summaryStatistic <- function( bin_n, bin_nPlus1, df ) {
    a <- df[ df$bin == bin_n, ] 
    b <- df[ df$bin == bin_nPlus1, ] 
    
    pValue  <- getPValue( a$splicingValueNorm, b$splicingValueNorm )
    
    return( c( mean( a$relPos), mean( a$splicingValueNorm ), var( a$splicingValueNorm ), pValue, readSum=sum( a$readSum ), binMembers=nrow( a ) ) )
  }
  
  # Calculate the summary statistics for each bin. 
  bins  <- as.character( unique( df$bin ) )
  statistics  <- data.frame( meanPos=rep(NA,length( bins ) ),
                             meanY=rep(NA,length( bins ) ), 
                             sdY=rep(NA,length( bins ) ), 
                             pValue=rep(NA,length( bins ) ), 
                             reads=rep(NA,length( bins ) ),
                             memberCount=rep(NA,length( bins ) ) )
  for( i in 2:length(  bins ) ) { statistics[ i-1, ] <- summaryStatistic( bins[ i - 1 ], bins[ i ], df  ) }
  
  # Correct pValue
  statistics$pValue <- p.adjust( statistics$pValue,method="hochberg" )
  
  return( statistics )
}

binAndSummarize_constant_bin_width <- function( file, id, remove_positions_after_exonEnd=F, bin_width=20 ) {
  
  # Read input
  df  <- read.table(file=file,sep="\t",head=T)
  
  # Remove positions after exon if specified 
  if(remove_positions_after_exonEnd) { df <- df[ df$relPos <= (getTerminalExonEnd(id) - 5), ] }
  
  # Filter read counts. Remove all values upstream of the 3'SS 
  minReadCount  <- 10
  df <- df[ df$readSum >= minReadCount, ]
  df  <- df[ df$relPos >= 0, ]
  
  # Bin with constant bin width.
  df$bin <- cut(df$relPos, breaks=seq(min(df$relPos)-1,max(df$relPos)+bin_width,bin_width))
  
  # Get summary statistic of specified bin: meanPos, meanY, sdY and pValue to next bin 
  summaryStatistic <- function( bin_n, bin_nPlus1, df ) {
    a <- df[ df$bin == bin_n, ] 
    b <- df[ df$bin == bin_nPlus1, ] 
    
    pValue  <- getPValue( a$splicingValueNorm, b$splicingValueNorm )
    
    return( c( mean( a$relPos), mean( a$splicingValueNorm ), var( a$splicingValueNorm ), pValue, readSum=sum( a$readSum ), binMembers=nrow( a ) ) )
  }
  
  # Calculate the summary statistics for each bin. 
  bins  <- as.character( unique( df$bin ) )
  statistics  <- data.frame( meanPos=rep(NA,length( bins ) ),
                             meanY=rep(NA,length( bins ) ), 
                             sdY=rep(NA,length( bins ) ), 
                             pValue=rep(NA,length( bins ) ), 
                             reads=rep(NA,length( bins ) ),
                             memberCount=rep(NA,length( bins ) ) )
  for( i in 2:length(  bins ) ) { statistics[ i-1, ] <- summaryStatistic( bins[ i - 1 ], bins[ i ], df  ) }
  
  # Correct pValue
  statistics$pValue <- p.adjust( statistics$pValue,method="hochberg" )
  
  return( statistics )
}



# Get position of terminal exon. 
getTerminalExonEnd <- function( id ) {
  xEnd <- terminalExonLengths[ terminalExonLengths[,1] == id, 2]
  return( xEnd )
}

# Get Saturation point defined as the last 4 bins 
getSaturation <- function( statistics, id ) {
  
  # Define the last position. 
  xEnd <- max( statistics$meanPos, na.rm=T )   
  exonEnd  <- getTerminalExonEnd( id )
  if( exonEnd < xEnd ) { xEnd  <- exonEnd }
  
  # Remove bins with lowest 10% members
  minMember  <- ceiling( quantile( statistics$memberCount, na.rm=T, probs=c( 0.1 ) )[[ 1 ]] )
  df  <- statistics[ which( statistics$memberCount >= minMember ), ]
  
  lastBin <- max( which( !df$meanPos >= xEnd ) )
  yEnd <- mean( statistics[ (lastBin-1):lastBin, ]$meanY )
  
  return( c( xEnd, yEnd ) ) 
}

getCrossingPoint <- function( statistics, id, fractionSaturation=0.5, infer_crossing_point=T ) {
  saturation  <- getSaturation( statistics, id )
  
  # Remove bins with lowest 10% members, don't touch the first 2 bins
  minMember  <- ceiling( quantile( statistics$memberCount, na.rm=T, probs=c( 0.1 ) )[[ 1 ]] )
  bins_to_remove <- which( statistics$memberCount < minMember | is.na(statistics$memberCount))
  bins_to_remove  <- bins_to_remove[ bins_to_remove >= 3 ]
  df  <- statistics[ -bins_to_remove, ]
  
  y_target  <-  saturation[ 2 ]*fractionSaturation
  
  if( y_target > 0 ) {
    bin_crossed <- min( which( df$meanY > y_target ) )
    
    bin_pre_cross <- (bin_crossed - 1) 
    if(infer_crossing_point) {
      if( bin_pre_cross > 0 ) {
        x <- df$meanPos[ (bin_crossed - 1 ): bin_crossed  ]
        y <- df$meanY[ (bin_crossed - 1 ): bin_crossed  ]  
      } else {
        x_pre <- max(0,df$meanPos[bin_crossed] - (df$meanPos[(bin_crossed+1)] - df$meanPos[bin_crossed]))
        x <- c(x_pre, df$meanPos[bin_crossed] )
        y <- c(0, df$meanY[bin_crossed])
      }       
    } else {
      x <- df$meanPos[ (bin_crossed - 1 ): bin_crossed  ]
      y <- df$meanY[ (bin_crossed - 1 ): bin_crossed  ]  
    }
    
    model <- lm( y ~ x )
    x_target <- ( y_target - model$coefficients[[ 1 ]] ) / model$coefficients[[ 2 ]]
    return( data.frame( x=x_target, y=y_target ) )  
  } else {
    return( data.frame( x=NA, y=0 ) )
  }
  
}

# Plot splicing bins vs. positions. 
plot.splicingValues <- function( statistics, id, model=NA, write=F, pathOut="", xrange=NA ) {
  
  exonEnd  <- getTerminalExonEnd( id )
  saturation  <- getSaturation( statistics, id )
  
  splicingValuesReached  <- data.frame( getCrossingPoint( statistics, id, 0.5 ) ) 
  
  df <- statistics
  df$label  <- rep("data",nrow(df) ) 
  
  if( !is.na( model ) ) {
    xPredict <- seq(from=0,to=max( df$meanPos, na.rm=T ), by=1 )
    yPredict <- predict( object=model, newdata=data.frame( x=xPredict ) )
    dfPredict <- data.frame( meanPos=xPredict, meanY=yPredict, sdY=rep( 0, length( xPredict ) ), pValue=rep(NA,length( xPredict ) ), reads=rep(NA,length( xPredict ) ),  memberCount=rep(NA,length( xPredict ) ), label=rep("predict",length( xPredict ) ) )
    
    df <- rbind( df, dfPredict )
  }
  
  limits  <- aes( ymax=meanY+sdY,ymin=meanY-sdY)
  title <- paste( id, "splicingValue" )
  p  <- ggplot( df, aes(x=meanPos,y=meanY ) ) # Read data 
  p  <- p + ylab(label="fraction co-txn splicing") + ylim(0,1) # Label the y axis
  p  <- p + xlab(label="pos (rel to 3'SS)")  + expand_limits(x=0) # Label the y axis
  p  <- p + ggtitle( title )
  
  # Plot data layer
  p  <- p + geom_point(data=df[ df$label == "data", ], color="red") 
  p  <- p + geom_errorbar(data=df[ df$label == "data", ], limits, width=10)
  
  # Plot model layer 
  p  <- p + geom_line(data=df[ df$label == "predict", ], color="black" )
  
  p <- p + geom_vline( xintercept = saturation[ 1 ], colour="red", linetype = "longdash" )
  p <- p + geom_hline( yintercept = saturation[ 2 ], colour="red", linetype = "longdash" )
  p <- p + geom_vline( xintercept = splicingValuesReached$x , colour="green", linetype = "longdash" )
  if( exonEnd <= max( statistics$meanPos, na.rm=T ) ) {
    p <- p + geom_vline( xintercept = exonEnd, colour="blue", linetype = "longdash" )
  }
  
  if( !is.na(xrange) ) {
    p <- p + expand_limits(x = xrange)
  } 
  
  if( write ) {
    fileOut  <- paste( c( pathOut, id, "_fit.pdf" ), collapse="")
    ggsave( filename=fileOut,plot=p, width=7, height=7, useDingbats=FALSE )
  }
  
  return( p )
}

# Plot pValues vs. positions. 
plot.pValues <- function( statistics, id, write=F ) {
  
  p <- NA 
  
  if( !all( is.na( statistics$pValue ) ) ) {
    
    # Assign an id for plotting 
    title <- paste( id, "pValue" )
    p  <- ggplot( statistics, aes(x=meanPos,y=pValue) )
    p <- p + ylab(label="pValue") + scale_y_log10() + annotation_logticks(sides = "l")
    p  <- p + xlab(label="pos (rel to 3'SS)") + expand_limits(x=0) # Label the y axis
    p  <- p + ggtitle( title )
    p  <- p + geom_point(color="red")
    
    if( write ) {
      fileOut  <- paste( c( id, "_pValue.png" ), collapse="")
      png( fileOut )
      show( p )
      dev.off()   
    }
  }
  
  return( p )
}


####################
# Methods dealing with the first observed splicing values
#####################

# Get entries of first occured splicing 
getFirstNSplicedEntries <- function( file, n ) {
  
  minReadCount  <- 10
  df  <- read.table(file=file,sep="\t",head=T)
  df <- df[ df$readSum >= minReadCount, ]
  df  <- df[ df$relPos >= 0, ]
  
  splicedX  <- df[ df$splicingValueRaw > 0, ]
  
  return( head( splicedX, n ) )
}

# Get first n spliced positions
getFirstNSplicedPositions <- function( file, n ) {
  df  <- getFirstNSplicedEntries( file, n )
  return( df[,1])
}

# Get first n splicing values 
getFirstNSplicingValue <- function( file, n ) {
  df  <- getFirstNSplicedEntries( file, n )
  return( df[,4])
}

######################
# Fitting models
######################

#
# 1. Assume splicing to occur with a constant rate k 
# 2. Assume a maximal saturation level of s_max
# 3. Assume a delay in signal x_delay
# co-txn splicing(x) = s_max - exp{ -k*( x - x_delay ) }
#
exponential  <- function( x, k, s_max, x_delay  ) {
  y <- ( 1 - exp( -k*( x - x_delay ) ) ) * s_max 
  y[ which( y < 0 ) ] <- 0
  return( y )
}

fitSingleExponential <- function( file, statistics ) {
  
  df <- statistics[ !is.na( statistics$meanPos ), ]
  
  x <- df$meanPos
  y <- df$meanY
  weights  <- df$sdY^-1
  weights[ is.infinite(weights) ]  <-  max( weights[ !is.infinite( weights ) ] )
  
  k1Start=0.01
  k2Start=0.0001
  s_max <-  getSaturation( statistics, id )[ 2 ]
  x_delay  <- getFirstNSplicedPositions( file, 1 )
  
  
  out <- nls( y ~ ( exponential( x, k, s_max, x_delay ) ), start=list( k=k1Start ), weights=weights, control = nls.control(maxiter=5000),algorithm="port",lower=list(k=0,s_max=0,delay=0),trace=F)
  
  return( out )
}

######################
# Summary statistics
######################

# Get basic read-count information on the splicingvalues file. 
read_count_summary <- function(file) {
  
  id  <- strsplit( file, split="_")[[1]][ 1 ]
  
  # Read input values, filter read counts. Remove all values upstream of the 3'SS 
  minReadCount  <- 10
  df  <- read.table(file=file,sep="\t",head=T)
  df <- df[ df$readSum >= minReadCount, ]
  df  <- df[ df$relPos >= 0, ]
  
  out  <- data.frame(gene=id, position_count=nrow(df), min_pos=min(df[["relPos"]]), max_pos=max(df[["relPos"]]), 
                     pos_mean_readcount=mean(df[["readSum"]]), pos_sd_readcount=sd(df[["readSum"]]))
  
  if( out$position_count > 3 ) {
    out <- cbind(out, data.frame(pos_to_pos_correlation=as.numeric(cor.test(df$splicingValueNorm[1:(nrow(df)-1)],df$splicingValueNorm[2:nrow(df)])$estimate) ))  
  } else {
    out <- cbind(out, data.frame(pos_to_pos_correlation=NA ))  
  }
  
  return(out)
}

summary_statistics <- function(file) {
  rcs  <- read_count_summary(file)
  id  <- strsplit( file, split="_")[[1]][ 1 ]
  binCount <- 30
  if( rcs$position_count >= binCount ) {
    statistics  <- binAndSummarize( file, id, remove_positions_after_exonEnd=T, binCount=binCount )    
    
    y <- statistics$meanY[!is.na(statistics$meanY)]
    df <- cbind( rcs, data.frame(bin_mean_readcount=mean(statistics$reads, na.rm=T), 
                                 bin_sd_readcount=sd(statistics$reads, na.rm=T), 
                                 bin_to_bin_sd=sd(statistics$meanY,na.rm=T), 
                                 bin_to_bin_correlation=as.numeric(cor.test(y[2:length(y)],y[1:length(y)-1])$estimate)))
  } else {
    df <- cbind( rcs, data.frame(bin_mean_readcount=NA, 
                                 bin_sd_readcount=NA, 
                                 bin_to_bin_sd=NA, 
                                 bin_to_bin_correlation=NA))
  }
  return(df)
}

