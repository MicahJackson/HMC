## split HMC 2
split.HMC <- function(x, y, q.hat, hessian, S,R, L, eps, n.iter){


     nx <- ncol(x)
     S.inv <- solve(S)
     q.mat <- matrix(c(q.hat),nrow=1)

     curr.q <- q.hat
     curr.U <- get.U(x, y, curr.q, S.inv)
     curr.dU1 <- dU1(x, y, curr.q, q.hat, S.inv, hessian)
     count <- 0

     for(n in 2:n.iter){

          move <- leap(x, y, curr.q, q.hat, S.inv, curr.dU1, hessian, R, nx,L,eps,count)

          curr.q <- move$q
          curr.U <- move$yoo
          curr.dU1 <- move$du1
          count <- move$ct

          q.mat <- rbind(q.mat,c(curr.q))
     }
     print(count/n.iter)
     return(q.mat)
}


leap <- function(x, y, current.q, q.hat,  S.inv, current.dU1, hessian, R, n.covs, L, eps, count){

     pos <- current.q
     mom <- t(rmvnorm(1,sigma=diag(length(pos))))
     du <- current.dU1

     curr.H <- get.U(x, y, pos, S.inv) + get.K(mom,diag(length(pos)))

     for(l in 1:L){
          mom <- mom - (eps/2)*du
          q.star <- current.q - q.hat
          xnot <- matrix(c(q.star,mom), ncol=1)
          go <- R%*%xnot
          pos <- go[1:n.covs] + q.hat
          mom <- go[(n.covs+1):(2*n.covs)]
          du <- dU1(x, y, pos, q.hat, S.inv, hessian)
          mom <- mom - (eps/2)*du
     }

     new.H <- get.U(x, y, pos, S.inv) + get.K(mom, diag(length(pos)))
     ratio <- new.H - curr.H
     accept <- min(1, exp(-ratio))
     u <- runif(1)

     if(u < accept){
          q <- pos; yoo <- get.U(x, y, pos, S.inv); du1 <- dU1(x, y, pos, q.hat, S.inv, hessian); count <- count+1
     }
     else{
          q <- current.q; yoo <- get.U(x, y, current.q, S.inv); du1 <- dU1(x, y, current.q, q.hat,S.inv, hessian)
     }

     return(list(q=q, yoo=yoo, du1=du1, ct=count))
}


get.U <- function(x,y,pos,S.inv){
     ## Potential Energy
     # Using normal (0,S) prior on parameters q (position). pre proc S.inv; serve in for S
     eta <- x%*%pos
     log.prior <- -(1/2)*(crossprod(pos, S.inv)%*%pos)
     log.like <- sum(y*eta-log(1+exp(eta)))

     return(-(log.prior+log.like))
}

get.K <- function(mom,M.inv){
     ## Kinetic Energy from momentum dbn normal(0,M)
     # define M.inv outside of this fcn if stationary
     return(-(1/2)*(t(mom)%*%M.inv%*%mom))
}

dU1 <- function(x,y,q,hat,S.inv, J){
     mu <- exp(x%*%q)/(1+exp(x%*%q))
     S.inv%*%q - crossprod(x,(y-mu)) - t(J)%*%(q-hat)
}

#####
# with varying step size
split.HMCe <- function(x, y, q.hat, hessian, S,R, L, eps, n.iter, e.scale=FALSE){


     nx <- ncol(x)
     S.inv <- solve(S)
     q.mat <- matrix(c(q.hat),nrow=1)

     curr.q <- q.hat
     curr.U <- get.U(x, y, curr.q, S.inv)
     curr.dU1 <- dU1(x, y, curr.q, q.hat, S.inv, hessian)
     count <- 0

     for(n in 2:n.iter){
          move <- leape(x, y, curr.q, q.hat, S.inv, curr.dU1, hessian, R, nx,L,eps, e.scale,count)

          curr.q <- move$q
          curr.U <- move$yoo
          curr.dU1 <- move$du1
          count <- move$ct

          q.mat <- rbind(q.mat,c(curr.q))
     }
     print(count/n.iter)
     return(q.mat)
}


leape <- function(x, y, current.q, q.hat,  S.inv, current.dU1, hessian, R, n.covs, L, eps, e.scale,ct){

     pos <- current.q
     mom <- t(rmvnorm(1,sigma=diag(length(pos))))
     du <- current.dU1

     curr.H <- get.U(x, y, pos, S.inv) + get.K(mom,diag(length(pos)))
     e <- ifelse(e.scale, runif(1,e.scale*eps,eps), eps)
     for(l in 1:L){
          mom <- mom - (e/2)*du
          q.star <- current.q - q.hat
          xnot <- matrix(c(q.star,mom), ncol=1)
          go <- R%*%xnot
          pos <- go[1:n.covs] + q.hat
          mom <- go[(n.covs+1):(2*n.covs)]
          du <- dU1(x, y, pos, q.hat, S.inv, hessian)
          mom <- mom - (e/2)*du
     }

     new.H <- get.U(x, y, pos, S.inv) + get.K(mom, diag(length(pos)))
     ratio <- new.H - curr.H
     accept <- min(1, exp(-ratio))
     u <- runif(1)

     if(u < accept){
          q <- pos; yoo <- get.U(x, y, pos, S.inv); du1 <- dU1(x, y, pos, q.hat, S.inv, hessian); ct=ct+1
     }
     else{
          q <- current.q; yoo <- get.U(x, y, current.q, S.inv); du1 <- dU1(x, y, current.q, q.hat,S.inv, hessian)
     }

     return(list(q=q, yoo=yoo, du1=du1,ct=ct))
}

#####
# tryouts

####
# simulated data
I <- 10000
J <- 10
Sigma <- diag(10)# sd for simmed observations



data.mat <- matrix(NA,ncol=J,nrow=I)
for(j in 1:J){
     for(i in 1:I){
          data.mat[i,] <- rmvnorm(1,sigma=Sigma)
     }}
data.mat <- cbind(1,data.mat)

set.seed(2)
betas <- rnorm(11)

y.sim <- numeric(I)
for(i in 1:I){
     y.sim[i] <- rbinom(1,1,(exp(data.mat[i,]%*%betas)/(1+exp(data.mat[i,]%*%betas))))
}


# prior info
S <- 100*diag(11)
M <- diag(11)
S.inv <- solve(S)
M.inv <- solve(M)

q.hat <- post.logistic.MAP(data.mat,y.sim,rep(.01,11))
hess <- logistic.hessian(data.mat,q.hat) + diag(c(S.inv%*%q.hat))

nc <- ncol(data.mat)
A <- matrix(0,2*nc,2*nc)
A[1:nc,(nc+1):(2*nc)] <- diag(nc)
A[(nc+1):(2*nc),1:nc] <- -hess

e <- .004

d <- e*eigen(A)$values
ed <- exp(d)
big.D <- diag(ed)

gamma <- eigen(A)$vectors
g.inv <- solve(gamma)
R <- gamma%*%big.D%*%g.inv ### changes with epsilon!!!!

std.sim <- logistic.HMC(data.mat, y.sim, q.hat, S, M, 20, .00013,2000)
system.time(std.sim <- logistic.HMC(data.mat, y.sim, q.hat, S, M, 20, .00013,2000))
getPerformance(data.mat, y.sim, std.sim, 48.782, 2000)

split.sim <- split.HMC(data.mat,y.sim,q.hat, hess, S, R, L=2, eps=.004,2000)
system.time(split.HMC(data.mat,y.sim,q.hat, hess, S, R, L=2, eps=.004,2000))
getPerformance(data.mat, y.sim, split.sim, 8.254, 2000)

splite.sim <- split.HMCe(data.mat,y.sim,q.hat, hess, S, R, L=2, eps=.004,2000,.9)
system.time(split.HMCe(data.mat,y.sim,q.hat, hess, S, R, L=2, eps=.004,2000,.9))
getPerformance(data.mat, y.sim, splite.sim, 8.066, 2000)



#####
## setup for example run
haberman <- read.csv('~/Documents/school/Winter\ 17/Computing/haberman.csv', header=FALSE)
names(haberman) <- c('age','year','nodes','survive5')
haberman$survive5[haberman$survive5==2] <- 0
summary(glm(survive5 ~ age + nodes, data=haberman, family='binomial'))
scaled <- scale(haberman[,c(1,3)])
scaled <- cbind(1,scaled)
X <- cbind(1,haberman[,c(1,3)])
Y <- haberman$survive5
S <- diag(rep(100,3))
S.inv <- solve(S)
M <- diag(rep(1,3))

q.hat <- post.logistic.MAP(scaled,Y,rep(.01,3))
hess <- logistic.hessian(scaled,q.hat) + diag(c(S.inv%*%q.hat))

nc <- ncol(scaled)
A <- matrix(0,2*nc,2*nc)
A[1:nc,(nc+1):(2*nc)] <- diag(nc)
A[(nc+1):(2*nc),1:nc] <- -hess

e <- .008

d <- e*eigen(A)$values
ed <- exp(d)
big.D <- diag(ed)

gamma <- eigen(A)$vectors
g.inv <- solve(gamma)
R <- gamma%*%big.D%*%g.inv ### changes with epsilon!!!!

####

std.hab <- logistic.HMC(scaled, Y, q.hat,S,M=diag(3),20,.003, 2000)
system.time(logistic.HMC(scaled, Y, q.hat,S,M=diag(3),20,.003, 2000))
getPerformance(scaled, Y, std.hab, 2.268,2000)

split.hab <- split.HMC(scaled, Y, q.hat, hess, S, R,2, .2,2000)
system.time(split.HMC(scaled, Y, q.hat, hess, S, R,2, .2,2000))
## with trajectory length .2, (2,.1) is best combo with 35 acceptance.
## using (2,.2) is much better. 90 acceptance

splite.hab <- split.HMCe(scaled, Y, q.hat, hess, S, R,8, .008,2000,.9)
system.time(split.HMCe(scaled, Y, q.hat, hess, S, R,8, .008,2000,.9))
## performs similarly with e.scale .9. If it's .8, performs worse
getPerformance(scaled, Y, splite.hab, 1.306,2000)





