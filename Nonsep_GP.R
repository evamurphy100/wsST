#Simulate a Gaussian process with mean zero and a non-separable Gneiting correlation function

Gneiting<-function(h,u, al = 0.75, s2 = 0.968, a=0.5, c=1,beta)
{
  #al = temporal smoothness parameter
  #s2 = variance of the spatio-temporal process
  # a,c = temporal and spatialscaleing parameters
  #beta = space-time interaction
  part1=a*abs(u)^(2*al)+1;
  denom=part1^(beta/2);
  part2=exp(-c*h/denom);
  return(s2*part2/part1);
}

### Simulate the process with the given correlation function
# Time points
time <- 1:6
nt = 6
#spatial locations
loc <- cbind(expand.grid(seq(-83, -74, len = 16), seq(30, 36, len = 16)))
locdist <- rdist.earth(loc)
ns = dim(loc)[1]
#Compute the Covariance
Sigma=matrix(0,ns*nt,ns*nt)
for(space1 in 1:ns){
    for(space2 in 1:ns)
    {
      hlag=locdist[space1,space2]
      for(time1 in 1:nt)
      {
        for(time2 in 1:nt)
        {
          ulag=time2 - time1
          index1=(time1-1)*ns+space1;
          index2=(time2-1)*ns+space2;      
          Sigma[index1,index2]=Gneiting(h = hlag, u=ulag, beta = .6)
        }
      }
    }
  }

#Simulate the covariance
set.seed(123)
xx.gneiting <- crossprod(chol(Sigma), rep(rnorm(nrow(loc) * 6)))
## Create the time spatial data frame
simdf.gneiting <- data.frame(lon = rep(loc[,1], each = 6), 
                             lat = rep(loc[,2], each = 6), 
                            dat= xx.gneiting, 
                            time = rep(1:6,  nrow(loc)) )

# plot
lattice::levelplot(dat~ lon + lat|time, data = simdf.gneiting)
