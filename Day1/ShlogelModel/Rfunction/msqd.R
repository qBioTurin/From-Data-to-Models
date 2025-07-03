
msqd<-function(reference, output)
{
  reference[1,] -> times_ref
  reference[2,] -> X1_ref
  
  # We will consider the same time points
  X1 <- output[which(output$Time %in% times_ref),"X1"]
  X1_ref <- X1_ref[which( times_ref %in% output$Time)]
  
  diff.X1 <- 1/length(times_ref)*sum(( X1 - X1_ref )^2 )
  
  return(diff.X1)
}
