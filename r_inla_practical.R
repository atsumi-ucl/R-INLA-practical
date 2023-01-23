library(INLA)
data(SPDEtoy)
names(SPDEtoy)
dim(SPDEtoy)
coords<- as.matrix(SPDEtoy[,1:2]); p5<-coords[1:5,]
pl.dom<-cbind(c(0,1,1,0.7, 0),c(0,0,0.7, 1,1))
args(inla.mesh.2d)
mesh1<- inla.mesh.2d(p5, max.edge = c(0.5, 0.5))
dev.new(noRStudioGD = T)
plot(mesh1)
mesh2 <- inla.mesh.2d(p5, max.edge = c(0.5, 0.5), cutoff=0.2)
plot(mesh2)
mesh3 <- inla.mesh.2d(p5, max.edge = c(0.1, 0.5), cutoff=0.1)
mesh4 <- inla.mesh.2d(p5, max.edge = c(0.1, 0.5), offset=c(0, -0.65))
mesh5 <- inla.mesh.2d(, pl.dom, max.edge = c(0.3, 0.5), offset=c(0.03, 0.5))
mesh6 <- inla.mesh.2d(,loc.domain = pl.dom, max.edge = c(0.3, 0.5), offset=c(0.03, 0.5), cutoff=0.1)
mesh7 <- inla.mesh.2d(, loc.domain=pl.dom, max.edge=c(0.3, 0.5), n=5, offset=c(.05,.1))  
mesh8 <- inla.mesh.2d(, loc.domain=pl.dom, max.edge=c(.3, 0.5), n=7, offset=c(.01,.3))  
mesh9 <- inla.mesh.2d(, loc.domain=pl.dom, max.edge=c(.3, 0.5), n=4, offset=c(.05,.3))  


non_conves_bdry <- inla.nonconvex.hull(p5, -0.3)
mesh10<- inla.mesh.2d(boundary = non_conves_bdry, max.edge = c(0.5 ,1), offset=c(0.2, 0.2), cutoff=0.2)

par(mfrow = c(3, 4), mar=c(0,0,1,0))

for (i in 1:10) {
  plot(pl.dom, type='l', col=3, lwd=2*(i>4), xlim=c(-0.57, 1.57),
       main = paste('mesh', i, sep=''), asp=1, axes=FALSE)
  plot(get(paste('mesh', i, sep='')), add=TRUE)
  points(p5, pch=19, col=2)
}

args(inla.spde2.matern)

mesh11<- inla.mesh.2d(, loc.domain=pl.dom, max.edge = c(0.092, 0.2))
spde11 <- inla.spde2.pcmatern(
  mesh=mesh11, alpha=2, 
  prior.range=c(0.3, 0.5), 
  prior.sigma = c(1, 0.01))

A11 <- inla.spde.make.A(mesh11, loc=coords)
dim(A11)

table(rowSums(A11>0))
table(rowSums(A11))

stkll <- inla.stack(data=list(resp=SPDEtoy$y), A=list(A11, 1),
                    effects=list(i=1:spde11$n.spde,
                                 b0=rep(1, nrow(SPDEtoy))))
res11 <- inla(resp ~ 0 + b0 + f(i, model=spde11),
              data=inla.stack.data(stkll),
              control.predictor=list(A=inla.stack.A(stkll),compute=T))
names(res11)
res11$summary.fixed
