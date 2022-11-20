
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


calculateObjective=function(w,Q,b){
  obj = 1/2*t(w)%*%Q%*%w - t(b)%*%w
  obj[1,1]
}

calculateGradient=function(w,Q,b){
  g = w%*%Q-b
  g
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

proxL1=function(x, alpha, lambda){
  out = sign(x)*pmax(0, abs(x)-alpha*lambda)
  out
}

proxElastic=function(x, alpha, lambda, weight){
  t = x/((1-weight)*lambda*alpha+1) # need to figure out where alpha sets in
  out = sign(t)*pmax(0, abs(t)-(weight*lambda*alpha)/((1-weight)*lambda*alpha+1))
  out
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


estimateLatentCorMatrix=function(data){
  E = matrix(NA,nrow=ncol(data),ncol=ncol(data))
  for (i in 1:ncol(data)){
    for (j in i:ncol(data)){
      if (i == j){
        E[i,j] = 1
      } else {
        m = cbind(data[,i], data[,j])
        cor.res = cor(m,method="kendall",use="pairwise")
        E[i,j] = sin(pi/2*cor.res[1,2])
        E[j,i] = sin(pi/2*cor.res[1,2])
      }
    }
  }
  E
}
