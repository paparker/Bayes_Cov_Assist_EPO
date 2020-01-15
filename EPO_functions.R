epoTransform <- function(M, Ref, nfact, AddM=NULL){
  # M is the matrix containing external variation
  # Ref is the reference matrix with no external variation
  # nfact is the number of components to orthogonalize against
  # AddM is an optional matrix that can be EPO transformed
  
  ## Standardize matrices
  mRef <- apply(Ref,1,mean)
  sdRef <- apply(Ref,1,sd)
  Ref <- (Ref-mRef)/sdRef
  
  mM <- apply(M,1,mean)
  sdM <- apply(M,1,sd)
  M <- (M-mM)/sdM
  
  ## Find projection matrix 
  D <- M - Ref
  DTD <- t(D) %*% D
  SVD <- svd(DTD)
  Q <- SVD$v[,1:nfact] %*% t(SVD$v[,1:nfact])
  P <- diag(nrow(Q)) - Q ## EPO transform matrix
  
  ## Apply transformation
  newM <- M %*% P
  
  output <- list(EPO_M=newM, Proj=P)
  if(!is.null(AddM)){
    newAdd <- AddM %*% P
    output[["EPO_Add"]] <- newAdd
  }
  return(output)
}


