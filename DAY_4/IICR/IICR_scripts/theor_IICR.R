####  Theoretical IICR two genes ####

IICR2nislandsS <- function(nbOfIslands = 10, geneFlow = 1, timeVector){
  # Compute the IICR of two genes under an n-islands model
  #
  # Args:
  #   nbOfIslands: number of islands (n)
  #   geneFlow: Migration rate (M)
  #   timeVector: a vector of times to evaluate the IICR
  #
  # Returns
  #   a vector IICR(t) of the IICR values at t
  n <- nbOfIslands
  M <- geneFlow
  t <- timeVector
  cGamma <- M/(n-1)
  cAlpha <- 0.5 * (1 + n*cGamma + sqrt((1 + n*cGamma)^2 - 4*cGamma))
  cBeta <- 0.5 * (1 + n*cGamma - sqrt((1 + n*cGamma)^2 - 4*cGamma))
  return( ((1 - cBeta)*exp(-cAlpha*t) + (cAlpha - 1)*exp(-cBeta*t)) / 
    ((cAlpha - cGamma)*exp(-cAlpha*t) + (cGamma - cBeta)*exp(-cBeta*t)) )
}
