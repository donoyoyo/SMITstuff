####################
# Define class for SMIT data
####################

# Dependencies
#source("~/Desktop/Slow_mutant/expCdfFunction.R")
# source("C:/Users/taraa/Documents/_lab/SMIT/scripts/misc/expCdfFunction.R")

library(methods)
library(binom)

####################
# Exponential CDF function 
####################
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


setClass(
  Class="smitData", 
  representation(
    id="character", 
    relPosition="integer", 
    splicedLength="integer",
    unsplicedLength="integer",
    splicedReads="integer", 
    unsplicedReads="integer",
    readSum="integer",
    splicedValueRaw="numeric",
    splicedValueNorm="numeric",
    ciLwr="numeric",
    ciUpr="numeric", 
    splicing.rate="numeric"
  )
)

####################
# smitData specific functions
####################

require( binom )
# Read smit data from file
smit.read  <- function( file ) {
  a  <- read.table( file, sep='\t',head=T)
  smitData  <- smit.dataFrameToSmitData( a ) 
  smitData@id  <- strsplit(file, "\\.")[[ 1 ]][ 1 ]
  return( smitData )
}

# Generate a smitData from data frame with correct columns
smit.dataFrameToSmitData <- function( a ) {
  return( 
    new("smitData",
        relPosition=as.integer( a$relPosition ), 
        splicedLength=as.integer( a$splicedLength ),
        unsplicedLength=as.integer( a$unsplicedLength ),
        splicedReads=as.integer( a$splicedReads ), 
        unsplicedReads=as.integer( a$unsplicedReads ),
        readSum=as.integer( a$splicedReads + a$unsplicedReads ),
        splicedValueRaw=rep(as.integer( NA ),length(a$relPosition)),
        splicedValueNorm=rep(as.integer( NA ),length(a$relPosition)),
        ciLwr=rep(as.integer( NA ),length(a$relPosition)),
        ciUpr=rep(as.integer( NA ),length(a$relPosition)), 
        splicing.rate=as.numeric( NA )
    )
  )
}

# Set the splicing values for the smit data 
# To this end we calculate the raw values which is the binomila fraction of the spliced reads
# and the insert length normalized binomial fraction 
smit.setSplicingValues <- function( smitData, insertLengthCutOff, insertLengthDecayRate ) {
  
  ##### Set the unnormalized splicing value #####
  # Calculate the binomial fraction
  splicingValue.raw  <- function( splicedCount, unsplicedCount ) {
    return( splicedCount / (splicedCount + unsplicedCount ) )
  }
  smitData@splicedValueRaw  <- splicingValue.raw( smitData@splicedReads, smitData@unsplicedReads )
  
  ##### Set the insert length normalized splicing value #####
  # Calculate the normalized binomial fraction 
  splicingValue.norm  <- function( splicedCount, unsplicedCount, splicedLength, unsplicedLength, insertLengthCutOff, insertLengthDecayRate ) {
    # Define insert length probability function, the parameters are learned from the data
    probability.insert  <- function( x, insertLengthCutOff, insertLengthDecayRate ) {
      if( x <= insertLengthCutOff ) {
        return( 0 )
      } else {
        # Normalization factor is the integral from cutoff to infity 
        norm.factor <- ( exp( - insertLengthCutOff * insertLengthDecayRate ) ) / insertLengthDecayRate
        return( exp( -insertLengthDecayRate * x ) / norm.factor )
      }
    }
    
    prob.spliced  <- sapply( splicedLength, probability.insert, insertLengthCutOff, insertLengthDecayRate )
    prob.unspliced <- sapply( unsplicedLength, probability.insert, insertLengthCutOff, insertLengthDecayRate )
    
    
    splicingValueNorm  <- as.numeric( rep(NA,length( splicedCount ) ) )
    for( i in 1:length( splicedCount ) ) {
      if( prob.spliced[ i ] == 0 ) {
          #splicingValueNorm[ i ]  <- as.numeric( NA ) # Set to zero, for logistic regression. 
          splicingValueNorm[ i ]  <- 0 # Set to NA, for exponential fit. 
      } else if( splicedCount[ i ] == 0 & unsplicedCount[ i ] == 0 ) {
        splicingValueNorm[ i ]  <- as.numeric( NA )
      } else {
          if( splicedCount[ i ] == 0 ) {
            splicingValueNorm[ i ]  <- 0
          } else {
            count.norm.spliced  <- splicedCount[ i ] / prob.spliced[ i ]
            count.norm.unspliced  <- unsplicedCount[ i ] / prob.unspliced[ i ]
            splicingValueNorm[ i ]  <- ( count.norm.spliced / ( count.norm.spliced + count.norm.unspliced ) )
          }
      }
    }
    return( splicingValueNorm )
  }
  smitData@splicedValueNorm  <- splicingValue.norm( smitData@splicedReads, smitData@unsplicedReads, smitData@splicedLength, smitData@unsplicedLength, insertLengthCutOff, insertLengthDecayRate )
  
  
  
  # Copy the internal paramter and return
  s = smitData
  return( s )
}

# # Sets the confidence intervals of the normalized splicing value
# smit.setSplicedValueNormConfInt  <- function( smitData ) {
#   # calculates the ci from the binomial fraction (f) and the number of observations (n)
#   binomFraction.ci <- function(f,n) {
#     if( is.na( f ) | is.na( n ) ) {
#       out  <- NA
#     } else {
#       f  <- as.numeric( f )
#       n  <- as.numeric( n )
#       x  <-  round( f * n )
#       out  <-  binom.wilson(x,n)
#     }
#     return( out )
#   }
#   
#   #  Calculate confidence intervals 
#   for( j in 1:length(smitData@relPosition ) ) {
#     # fake add splicing value norm = 0 for all ocurrances with relPos < 0   
#     if( smitData@relPosition[ j ] < 0 ) {
#       smitData@splicedValueNorm[ j ]  <- 0
#     }
#     
#     ci  <- binomFraction.ci( smitData@splicedValueNorm[ j ], smitData@readSum[ j ] )
#     message( ci )
#     if( length( ci ) == 1  ) {
#       smitData@ciLwr[ j ]  <- NA
#       smitData@ciUpr[ j ]  <- NA
#     } else {
#       smitData@ciLwr[ j ]  <- ci$lower
#       smitData@ciUpr[ j ]  <- ci$upper
#     }
#   }
#   
#   s = smitData 
#   return( s )
# }

# Keep values within defined relative positions
smit.defineXRange  <- function( smitData, minX, maxX ) {
  subsetVector  <- smitData@relPosition <= maxX & smitData@relPosition >= minX 
  smitData  <- 
  return( 
    new("smitData",
        id=smitData@id,
        relPosition=smitData@relPosition[subsetVector], 
        splicedLength=smitData@splicedLength[subsetVector],
        unsplicedLength=smitData@unsplicedLength[subsetVector],
        splicedReads=smitData@splicedReads[subsetVector], 
        unsplicedReads=smitData@unsplicedReads[subsetVector],
        readSum=smitData@readSum[subsetVector],
        splicedValueRaw=smitData@splicedValueRaw[subsetVector],
        splicedValueNorm=smitData@splicedValueNorm[subsetVector],
        ciLwr=smitData@ciLwr[subsetVector],
        ciUpr=smitData@ciUpr[subsetVector],
        splicing.rate=smitData@splicing.rate
    )
  )
}

# Bin smit data by relative position 
smit.bin  <- function( smitData, binSize ) {
  # Fill position information 
  relPosition <- as.numeric( seq(min( smitData@relPosition ),max( smitData@relPosition )+binSize, binSize ) ) # create the relPosRow
  splicedLength  <- relPosition - (smitData@relPosition[1]-smitData@splicedLength[1] ) # create the splicedLength 
  unsplicedLength  <- relPosition - (smitData@relPosition[1]-smitData@unsplicedLength[1] ) 
  
  # Create data frame to hold values 
  b <- data.frame( relPosition, splicedLength, unsplicedLength, splicedReads=rep(0,length(relPosition)), unsplicedReads=rep(0,length(relPosition)) )

  
  # Fill spliced and unspliced read values 
  currentBinRow = 1
  currentMaxSize = b$relPosition[ currentBinRow + 1 ]
  for( j in 1:length( smitData@relPosition ) ) { #Add the read count info 
    while( smitData@relPosition[ j ] >= currentMaxSize ) {
      currentBinRow  <- currentBinRow + 1 
      currentMaxSize = b$relPosition[ currentBinRow + 1 ]
    } 
    b$splicedReads[ currentBinRow ]  <- b$splicedReads[ currentBinRow ] + smitData@splicedReads[ j ]
    b$unsplicedReads[ currentBinRow ]  <- b$unsplicedReads[ currentBinRow ] + smitData@unsplicedReads[ j ]
  }
  # Add sum 
  b$readSum  <- b$splicedReads + b$unsplicedReads

  # Create a smitData object from data 
  a   <- smit.dataFrameToSmitData( b )
  a@id  <- smitData@id

  return( a )
}

# Shift relative position by a defined integer
smit.shiftRelPosition <- function( smitData, shift ) {
  relPosition.new  <- as.integer( smitData@relPosition + shift )
  smitData@relPosition  <- relPosition.new 
  a = smitData 
  return( a )
}

# Fit exponential cdf to the data
smit.fitExpCDF <- function( smitData, weightPow=1 ) {
  x  <- as.numeric( smitData@relPosition )
  y  <- as.numeric( smitData@splicedValueNorm )
  w  <- smitData@readSum^weightPow
  
  mod  <- expCDF.fit(x,y,w)
  if( length( mod ) > 1 ) {
    smitData@splicing.rate  <- coef( mod )[[ "k" ]]
  }
  
  a = smitData
  
  return( a )
}

# Plot the smit data
smit.plotOld <- function( smitData, weight=1 ) 
  {
  x <- smitData@relPosition
  yObs  <- smitData@splicedValueNorm
  w  <- smitData@readSum

  symbols(x=x, y=yObs, circles=w^weight, inches=1/10, ann=F, bg="firebrick2", fg=NULL, ylim=c(0,1) )
  if( !is.na( smitData@splicing.rate ) ) {
    yPredict  <- expCDF.predict(x,smitData@splicing.rate) 
    points( x, yPredict , type='l' )  
  }
}

smit.plot  <- function( smitData, weight = 1 ) {
  require( ggplot2 )
  
  w  <- smitData@readSum^weight
  df  <- data.frame( x=smitData@relPosition, theta=smitData@splicedValueNorm,reads=w )
  p  <- ggplot( df, aes(x=x,y=theta,size=reads ) ) # Read data 
  p  <- p + ylab(label=expression(theta)) + ylim(0,1) # Label the y axis
  p  <- p + ggtitle( smitData@id )
  #p  <- p + geom_point(color="white",fill="red",shape=21) + scale_area(range=c(0,5)) + theme(legend.position="none")
  p  <- p + geom_point(color="red") + scale_size_area(max_size=5) + theme(legend.position="none")
  
  if( !is.na( smitData@splicing.rate ) ) {
    test  <- function(x) { return(expCDF(x,k=smitData@splicing.rate) ) }
    p  <- p + stat_function( fun=test ) 
  }
  
  show( p )

}

smit.writeSplicingValues <- function( smitData, fileOut=NA ) {
  df <- data.frame( relPos = smit.getRelPosition( smitData ),
                    readSum = smit.getReadSum( smitData ),
                    splicingValueRaw = smit.getSplicedValueRaw( smitData ), 
                    splicingValueNorm = smit.getSplicedValueNorm( smitData ) )
  
  if( is.na( fileOut ) )
  {
    fileOut  <- paste( smit.getId( smitData ), "_splicingValues.txt", collapse=".", sep="" )
  }
  
  write.table( df, fileOut, sep="\t", quote=F, row.names=F )
}

##############################
# Getter 
##############################
smit.getRelPosition  <- function( smitData ) { return( smitData@relPosition ) }
smit.getSplicingRate <- function( smitData ) { return( smitData@splicing.rate) }
smit.getSplicedValueRaw <- function( smitData ) { return( smitData@splicedValueRaw ) }
smit.getSplicedValueNorm <- function( smitData ) { return( smitData@splicedValueNorm ) }
smit.getReadSum  <- function( smitData ) { return( smitData@readSum ) }
smit.getId <- function( smitData ) { return( smitData@id ) }
