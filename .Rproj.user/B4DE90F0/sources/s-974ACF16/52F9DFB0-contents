# Test Latent Differential Graphical Model (LDGM) on simulation data

# Created by: Elysia Saputra
# Date: May 10, 2020

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.

#source("functions_preprocessing.R")
library(ggplot2)
library(gplots)
library(huge)
library(igraph)
source("functions_LDGM.R")

# Load datasets
X = read.table("simulated_graphs/sim20_1_1/save/data/data.1.txt",header=T)
Y = read.table("simulated_graphs/sim20_1_2/save/data/data.1.txt",header=T)

# Convert to non-paranormal distribution
X.npn = huge.npn(X)
Y.npn = huge.npn(Y)


LDGM=function(X, Y, lambda=1, alpha=0.001, maxIter=5000, seed=1, loss='L1', weight=NULL){
  E_X = estimateLatentCorMatrix(X)
  E_Y = estimateLatentCorMatrix(Y)

  Q = E_X %x% E_Y
  b = c(E_X-E_Y)
  
  out = accProximalGradientDescent(Q,b,alpha=alpha,seed=seed,loss=loss, weight=weight, lambda=lambda)
  diff_g = symmetrize(out$w, X)
  out$diff_g = diff_g
  out
}

out_test = LDGM(X.npn, Y.npn)
delta_out = symmetrize(out_test$w, X.npn)


### Step 2: Estimating the latent covariance matrix of X and Y
E_X = estimateLatentCorMatrix(X)
E_Y = estimateLatentCorMatrix(Y)

### Step 3: Optimization
Q = E_X %x% E_Y
b = c(E_X-E_Y)


line_search=function(x, alpha, beta=0.5){
  alpha_old = alpha
  termination_condition = FALSE
  while (!termination_condition){
    
  }
}
  
gradientDescent=function(Q,b,tol=10^(-8), alpha=0.01, maxIter=5000, lambda=0.9, seed=1){
  set.seed=seed
  w_old = rnorm(length(b), 0, 1)
  obj_all = NULL
  w_diff_all = NULL
  for (i in 1:maxIter){
    print(i)
    g = as.vector(calculateGradient(w_old,Q,b))
    w_new = w_old - alpha*g
    obj_all = c(obj_all, calculateObjective(w_new, Q, b))
    w_diff_all = c(w_diff_all, sum(abs(w_new-w_old)))
    if (i>1 && sum(abs(w_new-w_old)) < tol){
      break
    }
    w_old = w_new
    alpha= lambda*alpha
  }
  out = list('w'=w_new, 'obj'=obj_all, 'convergence'=w_diff_all)
}

out = gradientDescent(Q,b,seed=13)
obj = out$obj
convergence = out$convergence
plot(seq(1,length(obj),1), obj)
plot(seq(1,length(convergence),1), convergence)

w = out$w

delta = matrix(w, nrow = nrow(E_X), byrow = TRUE)
### Step 4: Symmetrizing the matrix
for(i in 1:(ncol(delta)-1)){
  for (j in (i+1):ncol(delta)){
    delta.ij = delta[i,j]
    delta.ji = delta[j,i]
    if (abs(delta.ij) >= abs(delta.ji)){
      delta[j,i] = delta.ij
    } else {
      delta[i,j] = delta.ji
    }
  }
}

symmetrize=function(w, X){
  delta = matrix(w, nrow = ncol(X), byrow = TRUE)
  ### Step 4: Symmetrizing the matrix
  for(i in 1:(ncol(delta)-1)){
    for (j in (i+1):ncol(delta)){
      delta.ij = delta[i,j]
      delta.ji = delta[j,i]
      if (abs(delta.ij) >= abs(delta.ji)){
        delta[j,i] = delta.ij
      } else {
        delta[i,j] = delta.ji
      }
    }
  }
  delta
}

delta.round = round(delta)
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "RdYlBu"))(75)
heatmap(out_test$diff_g, Rowv=NA,Colv="Rowv",symm=T, col = redgreen(75))



#### Proximal Gradient Descent ####

proxL1=function(x, alpha, lambda){
  out = sign(x)*pmax(0, abs(x)-alpha*lambda)
  out
}

proxElastic=function(x, alpha, lambda, weight){
  t = x/((1-weight)*lambda*alpha+1) # need to figure out where alpha sets in
  out = sign(t)*pmax(0, abs(t)-(weight*lambda*alpha)/((1-weight)*lambda*alpha+1))
  out
}

w_elastic = proxElastic(w_old - alpha*g,alpha,lambda, 0.4)
w_l1 = proxL1(w_old - alpha*g,alpha,lambda)
plot(w_old-alpha*g, w_l1)
plot(w_old-alpha*g, w_elastic)


tol=0.00001 
alpha=0.5
maxIter=5000
beta=0.5
lambda=0.9

proximalGradientDescent=function(Q,b,tol=10^(-8), alpha=0.01, maxIter=5000, lambda=0.9, beta=0.5, seed=1, loss='L1',weight=NULL){
  w_old = rnorm(length(b), 0, 1)
  obj_all = NULL
  w_diff_all = NULL
  for (i in 1:maxIter){
    obj_old = calculateObjective(w_old, Q, b)
    g = as.vector(calculateGradient(w_old,Q,b))
    if (loss=='L1'){
      w_new = proxL1(w_old - alpha*g,alpha,lambda)
    } else if (loss=='elastic'){
      w_new = proxElastic(w_old - alpha*g,alpha,lambda,weight)
    }
    obj_new = calculateObjective(w_new, Q, b)
    
    # Line search with majorization-minimization
    termination_condition = FALSE
    while (!termination_condition){
      obj_upper_bound = obj_old + g%*%(w_new-w_old) + (w_new-w_old)%*%(w_new-w_old)/(2*alpha)
      if (obj_new > obj_upper_bound){
        alpha = alpha*beta
        w_new = proxL1(w_old - alpha*g,alpha,lambda)
        obj_new = calculateObjective(w_new, Q, b)
      } else {
        termination_condition = TRUE
      }
    }
    
    obj_all = c(obj_all, obj_new)
    w_diff_all = c(w_diff_all, sum(abs(w_new-w_old)))
    if (i>1 && sum(abs(w_new-w_old)) < tol){
      break
    }
    w_old = w_new
  }
  out = list('w'=w_new, 'obj'=obj_all, 'convergence'=w_diff_all)
  out
}

accProximalGradientDescent=function(Q,b,tol=10^(-8), alpha=0.01, maxIter=5000, lambda=0.9, beta=0.5, seed=1, loss='L1',weight=NULL){
  w_old = rnorm(length(b), 0, 1)
  obj_all = NULL
  w_diff_all = NULL
  for (i in 1:maxIter){
    omega = i/(i+3)
    if (i == 1){
      w_old2 = w_old
    }
    y_new = w_old + omega*(w_old - w_old2)
    obj_old = calculateObjective(w_old, Q, b)
    g = as.vector(calculateGradient(y_new,Q,b))
    if (loss=='L1'){
      w_new = proxL1(y_new - alpha*g,alpha,lambda)
    } else if (loss=='elastic'){
      w_new = proxElastic(w_old - alpha*g,alpha,lambda,weight)
    }
    obj_new = calculateObjective(w_new, Q, b)
    
    # Line search with majorization-minimization
    termination_condition = FALSE
    while (!termination_condition){
      obj_upper_bound = obj_old + g%*%(w_new-y_new) + (w_new-y_new)%*%(w_new-y_new)/(2*alpha)
      if (obj_new > obj_upper_bound){
        alpha = alpha*beta
        w_new = proxL1(w_old - alpha*g,alpha,lambda)
        obj_new = calculateObjective(w_new, Q, b)
      } else {
        termination_condition = TRUE
      }
    }
    
    obj_all = c(obj_all, obj_new)
    w_diff_all = c(w_diff_all, sum(abs(w_new-w_old)))
    if (i>1 && sum(abs(w_new-w_old)) < tol){
      break
    }
    w_old = w_new
    w_old2 = w_old
  }
  out = list('w'=w_new, 'obj'=obj_all, 'convergence'=w_diff_all)
  out
}


out = proximalGradientDescent(Q,b,alpha=0.001,seed=3,loss='L1', weight=0.3)
out2 = accProximalGradientDescent(Q,b,alpha=0.001,seed=3,loss='elastic', weight=0.8)
obj = out$obj
obj2 =out2$obj
convergence = out$convergence
convergence2 = out2$convergence
plot(seq(1,length(obj),1), obj)

plot(seq(1,length(convergence),1), convergence)
plot(seq(1,length(convergence2),1), convergence2, add=T)

w_out = out2$w
delta = matrix(w_out, nrow = nrow(E_X), byrow = TRUE)
### Step 4: Symmetrizing the matrix
for(i in 1:(ncol(delta)-1)){
  for (j in (i+1):ncol(delta)){
    delta.ij = delta[i,j]
    delta.ji = delta[j,i]
    if (abs(delta.ij) >= abs(delta.ji)){
      delta[j,i] = delta.ij
    } else {
      delta[i,j] = delta.ji
    }
  }
}

delta.round = round(delta)
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "RdYlBu"))(75)
heatmap(delta.round, Rowv=NA,Colv="Rowv",symm=T, col = redgreen(75))






true.edges.X = read.delim("simulated_graphs/sim20_1_1/save/graph/graph.1.txt",header=F,sep=" ")
true.edges.X = true.edges.X[,c(2,4)]
true.edges.X = true.edges.X[-c(1,2,3),]

true.edges.Y = read.delim("simulated_graphs/sim20_1_2/save/graph/graph.1.txt",header=F,sep=" ")
true.edges.Y = true.edges.Y[,c(2,4)]
true.edges.Y = true.edges.Y[-c(1,2,3),]

true.diff.edge.mat = matrix(0,nrow=ncol(X),ncol=ncol(X))

intersect.edge.table = matrix(NA, nrow=0,ncol=2)
colnames(intersect.edge.table) = c("1","2")
for (i in 1:nrow(true.edges.X)){
  row.i = true.edges.X[i,]
  edges.i = c(as.character(row.i[1][[1]]), as.character(row.i[2][[1]]))
  edges.i.edited = sub("X","",edges.i)
  edges.i.edited = as.numeric(edges.i.edited)
  diff.edge = TRUE
  for (j in 1:nrow(true.edges.Y)){
    row.j = true.edges.Y[j,]
    edges.j = c(as.character(row.j[1][[1]]), as.character(row.j[2][[1]]))
    edges.j.edited = sub("X","",edges.j)
    edges.j.edited = as.numeric(edges.j.edited)
    
    intersect.edge = intersect(edges.i.edited, edges.j.edited)
    if (length(intersect.edge) == 2){
      diff.edge = FALSE
    }
  }
  if (diff.edge == TRUE){
    edges.i.ind = as.integer(sub("X","",edges.i))
    true.diff.edge.mat[edges.i.ind[1],edges.i.ind[2]] = 1
    true.diff.edge.mat[edges.i.ind[2],edges.i.ind[1]] = 1
  }
}

delta.Z = (delta-mean(delta))/sd(delta)

heatmap(delta.round,Rowv=NA,Colv="Rowv",symm=T)
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "RdYlBu"))(75)
heatmap(true.diff.edge.mat, Rowv=NA,Colv="Rowv",symm=T, col = redgreen(75))








