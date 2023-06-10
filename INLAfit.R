library(INLAspacetime)
library(INLA)
library(inlabru)
library(fields)
library(maps)
library(viridisLite)
library(gtools)
library(mapdata)
load("Wind_SC.RData")
#load("Wind_smaller.RData")
# Save the locations
loc = cbind(c(lonSE_SC), c(latSE_SC))
#loc = cbind(c(lonSE_smaller), c(latSE_smaller))

## Number of time points used
## Subsets the data
t.max = 50
# Set-up the data frame
#ws =c(r22.smaller[,,1:t.max,1])
ws = c(r22.SC[,,1,1:t.max])
ws.df <- data.frame(lon = loc[,1], 
                    lat = loc[,2], 
                    ws= ws, 
                    time = 1:t.max )

## Max edge in spatial mesh
## Small number makes algorithm very slow
## ex: 2, 1.5
max.e = 2

## Prior model information
M = list()
for (i in 1:2) M[[i]] = list()
## Posteriors
fits = list()
## stacks (see later)
stack = list()

# Create temporal mesh
mesh.t = inla.mesh.1d(1:t.max)

#Create spatial mesh
mesh.s = inla.mesh.2d(loc = cbind(ws.df$lon, ws.df$lat), 
                      max.edge=c(1, 2),
                      cutoff=max.e/5)

plot(mesh.s);axis(1);axis(2)
points(ws.df$lon, ws.df$lat, pch = 16,cex = 0.25, col=adjustcolor("black", alpha.f = 0.5))
maps::map("state", add=TRUE, col = grey(.5))
mesh.s$n*mesh.t$n
# Call the functions for the nonseprable method from Haakon's website
download.file("https://raw.githubusercontent.com/haakonbakkagit/haakonbakkagit.github.io/master/functions-rgeneric-121-march2020.R", destfile="functions-rgeneric-121-march2020.R")
source("functions-rgeneric-121-march2020.R")
## Separable method
## Model component in space
mco.space = inla.spde2.pcmatern(mesh = mesh.s,
                                prior.range = c(600, .6), prior.sigma = c(18, 0.01))


#stack for separable model
#A.s1 = inla.spde.make.A(mesh = mesh.s, loc = cbind(ws.df$lon, ws.df$lat))
iset = inla.spde.make.index('i', n.spde = mesh.s$n, n.group = mesh.t$m)
A.st = inla.spde.make.A(mesh = mesh.s, loc = cbind(ws.df$lon, ws.df$lat), group = ws.df$time, group.mesh =mesh.t) 

stack[[1]] = inla.stack(
  data = list(y = ws.df$ws), 
  effects = list(i = iset), 
  A = list(A.st), 
  tag = 'stdata') 


# Separable model formula
M[[1]]$shortname = "Separable"
## The separable model
#hyper22 = list(prec = list(initial = -2*log(5), fixed=T))
hyper.ar1.rho = list(rho = list(prior = "pc.cor1", param=c(0,0.9)))

## The separable formula
form.sep = y ~ -1 + f(i, model = mco.space, group = i.group, control.group = list(model="ar1", hyper = hyper.ar1.rho))

M[[1]]$formula = form.sep

fits[[1]] = inla(M[[1]]$formula,
                 family = "gaussian",
                 data =inla.stack.data(stack[[1]]),
                 control.predictor=list(compute=F, A=inla.stack.A(stack[[1]])))

## Non-seprable method
## Mesh in space and time
## Lambdas for exponential prior on transformed hyper-param (1/rt, 1/rs and sig)
rgen.obj = list(mesh.space = mesh.s,
                mesh.time = mesh.t,
                lambdas = c(1,1,5))
## Nonsep model definition
nm = mesh.s$n*mesh.t$n

## The non-separable random effect / random field
mco.nonsep = inla.rgeneric.define(
  model = stmodel121.interpret, debug = FALSE, n = nm, obj = rgen.obj)

i.nonsep = 1:(mesh.s$n * mesh.t$n)
A.st = inla.spde.make.A(mesh = mesh.s, loc = cbind(ws.df$lon, ws.df$lat), group = ws.df$time, group.mesh =mesh.t) 

stack[[2]] = inla.stack(
  data = list(y = ws.df$ws), 
  effects = list(i.nonsep = i.nonsep), 
  A = list(A.st), 
  tag = 'stdata') 

mco.nonsep.fix = inla.rgeneric.define(model = stmodel121.interpret, 
                                      debug = FALSE, n = nm, obj = rgen.obj)
## The formula
## It is possible to add a spatial random effect (uncomment the f(s))
form.nonsep = y ~ -1 + f(i.nonsep, model = mco.nonsep.fix, n = nm)

M[[2]]$formula = form.nonsep

fits[[2]] = inla(M[[2]]$formula,
                 family = "gaussian",
                 data =inla.stack.data(stack[[2]]),
                 control.predictor=list(compute=TRUE, A=inla.stack.A(stack[[2]])))
## Plot
plot_CanRCM4_SC <- function(dat, main, col = tim.colors(256), ...){
  par(las = 1)
  maps::map("state", xlim = range(lonSE_SC), ylim = range(latSE_SC), ...)
  image.plot(lonSE_SC, latSE_SC, dat, main = main,
           col = col, xlab = "Lon", ylab = "Lat", add = T,
           horizontal = T, ...)
  maps::map("state", xlim = range(lonSE_SC), ylim = range(latSE_SC),
      add = T, ...)
}

#Extract the space-time variables. We add on the intercept for easier comparison.
time.k = 1
rf.st.nonsep = fits[[2]]$summary.random$i$mean 
field.nonsep =rf.st.nonsep[1:mesh.s$n + (time.k-1)*mesh.s$n]
zlim.temp = range(field.nonsep)
xlim = c(min(lonSE_SC), max(lonSE_SC)); ylim = c(min(latSE_SC), max(latSE_SC))
proj = inla.mesh.projector(mesh.s, xlim = xlim, 
                             ylim = ylim, dims=c(25, 25))
field.nonsep.proj = inla.mesh.project(proj, field.nonsep)

pdf("Non_sep_Estimate_SC.pdf")
image.plot(list(x = proj$x, y=proj$y, z = field.nonsep.proj), col = tim.colors(256),
             xlim = xlim, ylim = ylim, zlim = c(0, 9))  
title("Estimated process using a non-separable covariance")
maps::map("state", add=TRUE, col = grey(.01))
dev.off()

rf.st.sep = fits[[1]]$summary.random$i$mean
time.k = 1

field =rf.st.sep[1:mesh.s$n + (time.k-1)*mesh.s$n]
zlim.temp.sep = range(field)
xlim = c(min(lonSE_SC), max(lonSE_SC)); ylim = c(min(latSE_SC), max(latSE_SC))
proj = inla.mesh.projector(mesh.s, xlim = xlim, 
                             ylim = ylim, dims=c(25, 25))
field.proj = inla.mesh.project(proj, field)

pdf("Sep_vs_True_SC_50yrs.pdf")
image.plot(list(x = proj$x, y=proj$y, z = field.proj), col = tim.colors(256),
             xlim = xlim, ylim = ylim,  zlim = c(0, 9))  
title("Estimated process using a separable covariance")
#points(ws.df$lon, ws.df$lat, col=adjustcolor("black", alpha.f = 0.1))
maps::map("state", add=TRUE, col = grey(.01))
dev.off()

pdf("True_Process_SC.pdf")
plot_CanRCM4_SC(r22.SC[,,1,51], main="", zlim = c(0, 9))
title("True process")
dev.off()
