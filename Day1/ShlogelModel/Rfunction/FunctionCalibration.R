#
# init_generation<-function(optim_vector)
# {
#     return(c(optim_vector[1],1e-4,0))
# }


X1Calibration<-function(x,n)
{
  return(x[1]*n)
}

k1Calibration<-function(x,n)
{
  return(x[2]*n)
}

k2Calibration<-function(x,n)
{
  return(x[3]*n)
}
