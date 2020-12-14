#' @title Linear interpolation function
#' @description A \code{nimbleFunction} emulating BUGS/JAGS's interp.lin. 
#' @param z value where the interpolation take place 
#' @param x numeric vector giving the coordinates of the points to be interpolated. 
#' @param y numeric vector giving the coordinates of the points to be interpolated.
#' @return interpolated value
#' @examples 
#' data(intcal20)
#' interpLin(4500,intcal20$CalBP,intcal20$C14Age)
#' # equivalent to:
#' approx(x=intcal20$CalBP,y=intcal20$C14Age,xout=4500)$y
#' @rdname interpLin
#' @import nimble 
#' @export

interpLin <- nimbleFunction(
  run = function(z = double(0),x=double(1),y=double(1)) {
    index.vector= which(x<z)
    index = index.vector[1]
    xp = x[index]
    xp1 = x[index-1]
    yp = y[index]
    yp1 = y[index-1]
    res = yp+(((yp1-yp)*(z-xp))/(xp1-xp))
    returnType(double(0))
    return(res)
  })    