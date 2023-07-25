#Simulate a Gaussian process with mean zero and a non-separable Gneiting correlation function

Gneiting<-function(h,u,a=.5,c=1,beta)
{
  part1=a*abs(u)+1;
  denom=part1^(beta/2);
  part2=exp(-c*h/denom);
  return(part2/part1);
}

## Set up the locations and time points
N = 16
x <- seq(0,3, len = N)
y <- seq(0,3, len = N)
loc = expand.grid(x,y)
plot(loc)
t = 1:30 # number of days
h = rdist(loc)
u = rdist(t)
bet = 1 
n = dim(loc)[1] # number of locations
p = length(t) # number of time points

## Generate the nonseparable data
data=integer(n*p)

# Generate the covariance set-up
Sigma11=matrix(0,n*p,n*p)
for(space1 in 1:n)
{
  for(space2 in 1:n)
  {
    hlag=h[space1,space2]
    for(time1 in 1:p)
    {
      for(time2 in 1:p)
      {
        ulag=u[time1, time2];
        index1=(time1-1)*n+space1;
        index2=(time2-1)*n+space2;      
        Sigma11[index1,index2]=Gneiting(h=hlag,u=ulag,beta=bet)
      }
    }
  }
}
## Simulate a Gaussian process at loc locations for 
library(parallel)
library(foreach)

cl = makeCluster(detectCores())
#Activate cluster for each library
registerDoParallel(cl)
clusterEvalQ(cl, {
  library(foreach)
  library(doParallel)
})
## Simulate 100 replicates
sim.ws.data = foreach(i = 1:100) %dopar% {
## simulate 30 years of a month
Z = replicate(30, rnorm(n*p, mean = rep(0, n*p), sd = ))
data.ws.sim = array(dim = c(7680, 30))
for( j in 1:30){
  #data=t(L11)%*%Z[,i]
  data.ws.sim[,j] = t(L11)%*%Z[,j]
}
Z.matrix = array(dim = c(256, 30, 30)) #30 days over 30 years
for(j in 1:30){
  Z.matrix[,,j] = matrix(data.ws.sim[,j], nrow = 256, ncol = 30)
}
Z.matrix
}
stopCluster(cl)
#Plot a few replicates
library(fields)
quilt.plot(loc, sim.ws.data[[1]][,1,1])
for(i in 1:10) {
  image.plot(x,y, matrix(sim.ws.data[[1]][,1,i], nrow = 16, ncol = 16))
}
