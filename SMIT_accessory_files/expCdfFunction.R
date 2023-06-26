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