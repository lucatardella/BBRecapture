#' Partial Capture Histories from a Binary Data Matrix
#'
#' \code{pch} is used to obtain all the observed partial capture histories corresponding to an observed binary data matrix.
#'
#' @usage pch(data.matrix)
#'
#' @param data.matrix a binary data matrix
#'
#' @return \code{pch} returns a matrix of mode \code{"character"} where each element represents the partial capture history associated to the respective element of the input binary data matrix.
#'
#' @seealso \code{\link{BBRecap.custom.part}} and \code{\link{LBRecap.custom.part}}
#'
#' @author Danilo Alunni Fegatelli and Luca Tardella
#'
#' @examples
#'
#' data(greatcopper) # load greatcopper data
#' pch(greatcopper)
#'
#' @export

pch=function(data.matrix){
partial.capture.histories=matrix(NA,nrow=nrow(data.matrix),ncol=ncol(data.matrix))

for(i in 1:nrow(data.matrix)){
for(j in 1:(ncol(data.matrix)-1)){


partial.capture.histories[i,1]=""
partial.capture.histories[i,j+1]=paste(data.matrix[i,1:j],collapse="")


	}}

out=partial.capture.histories
return(out)

}
