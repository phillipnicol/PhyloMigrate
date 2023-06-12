


#Distance matrix M
# Position

#' @export
estimateMigration <- function(cnr,
                              pos,
                              A,
                              nrep=100,
                              perturb=FALSE,
                              dist.method="euclidean",
                              boot=FALSE) {
  cnr <- cbind(rep(0,nrow(cnr)),cnr)
  colnames(cnr)[1] <- "PD"
  M <- as.matrix(dist(t(cnr),method=dist.method))

  Q <- A2Q(A)

  trw <- fastme.bal(M)
  trw <- root(trw,outgroup="PD",resolve.root=TRUE)
  trw <- drop.tip(trw,tip="PD")

  #trw$edge.length <- rep(1,length(trw$edge.length))
  cnr <- cnr[,-1]
  M <- M[-1,]; M <- M[,-1]
  if(is.null(colnames(cnr))) {
    labels <- 1:ncol(cnr)
  } else {
    labels <- colnames(cnr)
    true_pos <- pos
    names(true_pos) <- labels
  }
  #Prune
  wmp <- Sankoff(trw,pos,labels,A)
  min.nz <- min(wmp$meta[wmp$meta[,3] > 0,3])
  zero.rows <- which(wmp$meta[,3] < 10^{-10})
  print(length(zero.rows))
  if(perturb) {
    wmp$meta[,3] <- wmp$meta[,3] + 0.001
  } else if(length(zero.rows) > 0) {
    wmp$meta <- wmp$meta[-zero.rows,]
  }

  out <- optim(par=wmp$ntran/sum(wmp$meta[,3]),fn=fn,gr=gr,Q=Q,
        meta=wmp$meta,method="L-BFGS-B",lower=10^{-10},hessian=TRUE)
  alpha.mt <- out$par
  alpha_sd <- 1/out$hessian


  Ac <- A
  Q <- A2Q(Ac)
  wmp <- Sankoff(trw,pos,labels,Ac)
  alpha.pars <- wmp$ntran/sum(trw$edge.length)
  alpha.m <- alpha.mt

  if(boot == FALSE) {
    out <- list()
    out$alpha <- alpha.mt
    out$ntrans <- wmp$ntran
    return(out)
  }
  nreg <- length(unique(pos))
  true_gd <- 0
  sim_gd.mt <- rep(0,nrep)
  sim_gd.pars <- rep(0,nrep)
  sim_gd.m <- rep(0,nrep)
  markov <- rep(0,nrep)
  mp <- rep(0,nrep)
  for(j in 1:nrep) {
    pos.sim <- simulateData(alpha.mt,trw,root.pos=wmp$root.pos,nreg,A)
    names(pos.sim) <- trw$tip.label
    pos.sim <- pos.sim[names(true_pos)]
    Comp <- Compare(true_pos,pos.sim,M)
    sim_gd.mt[j] <- Comp$sim
    true_gd <- Comp$true

    pos.sim <- simulateData(alpha.m,trw,root.pos=wmp$root.pos,nreg,Ac)
    names(pos.sim) <- trw$tip.label
    pos.sim <- pos.sim[names(true_pos)]
    Comp <- Compare(true_pos,pos.sim,M)
    sim_gd.m[j] <- Comp$sim
    true_gd <- Comp$true

    pos.sim <- simulateData(alpha.pars,trw,root.pos=wmp$root.pos,nreg,Ac)
    names(pos.sim) <- trw$tip.label
    pos.sim <- pos.sim[names(true_pos)]
    Comp <- Compare(true_pos,pos.sim,M)
    sim_gd.pars[j] <- Comp$sim
    true_gd <- Comp$true

  }

  out <- list()
  out$alpha <- alpha.mt
  out$ntrans <- trans
  out$sim_gd.mt <- sim_gd.mt
  out$sim_gd.m <- sim_gd.m
  out$sim_gd.pars <- sim_gd.pars
  out$true_gd <- true_gd
  out$tran.matrix <- wmp$tran.matrix
  out$alpha_sd <- alpha_sd
  return(out)
}

Sankoff <- function(trw,region,labels,A) {
  # Get optimal score
  node.order <- order(node.depth(trw),decreasing=FALSE)


  ncell <- length(region)
  pseudoregion <- region
  for(i in 1:ncell) {
    true.i <- which(labels == trw$tip.label[i])
    pseudoregion[i] <- which(region[true.i]==colnames(A))
  }

  pseudoregion <- as.numeric(pseudoregion)

  #Initialization
  S <- matrix(Inf,nrow=length(node.order),ncol=nrow(A))
  for(i in node.order) {
    if(node.depth(trw)[i] == 1) {
      # It is a leaf
      S[i,pseudoregion[i]] <- 0
    }
  }

  # Now main loop
  for(i in node.order) {
    if(node.depth(trw)[i] == 1) {
      next
    }
    children <- trw$edge[which(trw$edge[,1] == i),2]

    for(j in 1:nrow(A)) {
      S[i,j] <- min(S[children[1],]+A[j,])+min(S[children[2],]+A[j,])
    }
  }


  out <- list()
  out$S <- S
  root <- node.order[length(node.order)]
  out$ntran <- min(S[root,])
  out$root.pos <- which.min(S[root,])

  ## Top down refinement
  state <- rep(0, length(node.order))
  state[root] <- which.min(S[root,])
  reverse.order <- order(node.depth(trw),decreasing=TRUE)
  for(i in reverse.order) {
    if(node.depth(trw)[i] == 1) {
      next
    }
    children <- trw$edge[which(trw$edge[,1] == i),2]
    min.s <- min(S[children[1],] + A[state[i],])
    state[children[1]] <- which(A[state[i],]+S[children[1],]==min.s)[1]
    min.s <- min(S[children[2],] + A[state[i],])
    state[children[2]] <- which(A[state[i],]+S[children[2],]==min.s)[1]
  }


  meta <- matrix(0,0,ncol=3)
  for(i in node.order) {
    v <- c(0,0,0)
    v[1] <- state[i]
    children <- trw$edge[which(trw$edge[,1] == i),2]
    ixs <- which(trw$edge[,1] == i)
    if(length(children)==0) {next}
    v[2] <- state[children[1]]
    #Edge length maybe
    v[3] <- trw$edge.length[ixs[1]]
    meta <- rbind(meta,v)
    v[2] <- state[children[2]]
    v[3] <- trw$edge.length[ixs[2]]
    meta <- rbind(meta,v)
  }
  out$meta <- meta
  return(out)
}

Newton <- function(meta, Q, alphastart) {
  alpha <- alphastart

  var <- 0

  while(TRUE) {
    alpha.new <- alpha-deriv(meta,alpha,Q)/second.deriv(meta,alpha,Q)
    if(abs(alpha.new-alpha) < 10^{-5}) {
      alpha <- alpha.new
      var <- -1/second.deriv(meta,alpha,Q)
      break
    }

    if(alpha.new < 0) {
      alpha <- alpha/2
    } else {
      alpha <- alpha.new
    }
    print(alpha)
  }
  out <- list()
  out$var <- var
  out$alpha <- alpha
  return(out)
}

log.lik <- function(meta, alpha, Q) {
  sum <- 0
  for(i in 1:nrow(meta)) {
    u <- meta[i,1]; v <- meta[i,2]
    sum <- sum+log(expm(alpha*meta[i,3]*Q)[u,v])
  }
  sum
}

fn <- function(alpha,meta,Q) {
  scores <- apply(meta,1,function(x) {
    u <- x[1]; v <- x[2]
    log(expm(alpha*x[3]*Q)[u,v])
  })
  -sum(scores)
}

gr <- function(alpha,meta, Q) {
  scores <- apply(meta,1,function(x) {
    u <- x[1]; v <- x[2]
    P <- expm(alpha*x[3]*Q)
    num <- (Q %*% P)[u,v]
    denom <- P[u,v]
    x[3]*num/denom
  })
  -sum(scores)
}

deriv <- function(meta, alpha, Q) {
  sum <- 0
  for(i in 1:nrow(meta)) {
    u <- meta[i,1]; v <- meta[i,2]
    P <- expm(alpha*meta[i,3]*Q)
    num <- (Q %*% P)[u,v]
    denom <- P[u,v]
    sum <- sum + meta[i,3]*num/denom
  }
  sum
}

second.deriv <- function(meta, alpha, Q) {
  sum <- 0
  Q2 <- Q %*% Q
  for(i in 1:nrow(meta)) {
    u <- meta[i,1]; v <- meta[i,2]
    w <- meta[i,3]
    P <- expm(alpha*meta[i,3]*Q)
    denom <- P[u,v]^2

    num1 <- w^2*P[u,v]*(Q2 %*% P)[u,v]
    num2 <- w*((Q %*% P)[u,v])^2
    sum <- sum + w*((num1-num2)/denom)
  }
  sum
}





countTransitions <- function(trw,region,labels,A) {
  node.order <- order(node.depth(trw),decreasing=FALSE)
  ncell <- length(region)
  pseudoregion <- rep("0",ncell+trw$Nnode)
  for(i in 1:ncell) {
    true.i <- which(labels == trw$tip.label[i])
    pseudoregion[i] <- region[true.i]
  }
  trans <- 0
  trans.matrix <- matrix(0, nrow=8, ncol=8)
  for(i in node.order) {
    if(node.depth(trw)[i] == 1) {next}
    #Get child
    children <- trw$edge[which(trw$edge[,1] == i),2]
    region.child <- pseudoregion[children]
    i1 <- which(unique(region)==region.child[1])
    i2 <- which(unique(region)==region.child[2])
    trans <- trans+A[i1,i2]
    print(A[i1,i2])
    pseudoregion[i] <- region.child[1]
  }
  print(trans.matrix)
  return(trans)
}

setParent <- function(node,region,trw) {
  children <- trw$edge[which(trw$edge[,1] == node),2]
  trans <- 0
  trans.vector <- c()
  for(j in children) {
    if(region[j] == "0") {
      out <- setParent(j,region,trw)
      region <- out$region
      trans <- trans + out$trans
      trans.vector <- c(trans.vector,out$tv)
      region[node] <- region[j]
    } else {
      if(region[node] != "0" & region[node]!=region[j]) {
        trans <- trans+1
        hgt <- node.depth.edgelength(trw)[node]
        trans.vector <- c(trans.vector,hgt)
      }
      region[node] <- region[j]
    }
  }
  out <- list()
  out$region <- region
  out$trans <- trans
  out$tv <- trans.vector
  return(out)
}


simulateData <- function(alpha, trw, root.pos,nreg,A) {
  root <- which.max(node.depth(trw))
  cells <- max(node.depth(trw))
  pos <- rep(0, length(node.depth(trw)))
  pos[root] <- root.pos

  pos <- simDataNode(root,alpha,trw,pos,nreg,A)
  return(pos[1:cells])
}

simDataNode <- function(node,alpha,trw,pos,nreg,A) {
  child <- which(trw$edge[,1]==node)
  if(length(child) == 0) {
    return(pos)
  }
  for(j in child) {
    EL <- trw$edge.length[j]
    nmove <- rpois(1,lambda=EL*alpha)
    #cat(node, " ", EL, " ", move, "\n")
    pos.cur <- pos[node]
    if(nmove == 0) {
      pos[trw$edge[j,2]] <- pos.cur
    } else {
      for(k in 1:nmove) {
        nbhd <- which(A[pos.cur,] == 1)
        pos.cur <- sample(nbhd,1)
      }
    }
    pos[trw$edge[j,2]] <- pos.cur
    pos <- simDataNode(trw$edge[j,2],alpha,trw,pos,nreg,A)
  }
  return(pos)
}


Compare <- function(true_pos, pos, D) {
  k <- 10

  #True GD
  ncell <- nrow(D)
  KNN <- t(apply(D, 1, order)[2:(k+1),])
  matches <- vapply(1:ncell,FUN.VALUE=integer(1), function(x) {
    sum(true_pos[x] == true_pos[KNN[x,]])
  })
  true.metric <- sum(matches)/(k*ncell)

  #Simulated GD
  matches <- vapply(1:ncell,FUN.VALUE=integer(1), function(x) {
    sum(pos[x] == pos[KNN[x,]])
  })
  sim.metric <- sum(matches)/(k*ncell)

  out <- list()
  out$sim <- sim.metric; out$true <- true.metric
  return(out)
}


accuracy <- function(pos.true, pos.sim) {
  correct <- 0
  for(i in 1:1000) {
    cells <- sample(1:length(pos.true),size=2,replace=FALSE)
    region1 <- pos.true[cells[1]]
    region2 <- pos.true[cells[2]]
    if(region1==region2) {
      if(pos.sim[cells[1]] == pos.sim[cells[2]]) {
        correct <- correct + 1
      }
    } else {
      if(pos.sim[cells[1]] != pos.sim[cells[2]]) {
        correct <- correct + 1
      }
    }
  }
  return(correct/1000)
}

A2Q <- function(A) {
  Q <- A
  for(i in 1:nrow(A)) {
    ixs <- which(A[i,] == 1)
    Q[i,ixs]=1/length(ixs)
    Q[i,-ixs] <- 0
    Q[i,i] <- -1
  }
  return(Q)
}



#Distance matrix M
# Position

#' @export
estimateMigrationCov <- function(cnr,
                                pos,
                                A.list,
                                patients,
                                dist.method="euclidean",
                                pat.x) {

  meta <- matrix(0,nrow=0,ncol=5)
  pats <- unique(patients)
  cntr <- 1
  Q.list <- list()
  for(pa in pats) {
    ixs <- which(patients==pa)
    cnr.sub <- cnr[,ixs]
    cnr.sub <- cbind(rep(0,nrow(cnr.sub)),cnr.sub)
    colnames(cnr.sub)[1] <- "PD"
    M <- as.matrix(dist(t(cnr.sub),method=dist.method))

    Q <- A2Q(A.list[[cntr]])

    trw <- fastme.bal(M)
    trw <- root(trw,outgroup="PD",resolve.root=TRUE)
    trw <- drop.tip(trw,tip="PD")

    cnr.sub <- cnr.sub[,-1]
    M <- M[-1,]; M <- M[,-1]
    labels <- colnames(cnr.sub)
    true_pos <- pos[ixs]
    names(true_pos) <- labels

    wmp <- Sankoff(trw,pos[ixs],labels,A.list[[cntr]])
    min.nz <- min(wmp$meta[wmp$meta[,3] > 0,3])
    zero.rows <- which(wmp$meta[,3] == 0)
    print(length(zero.rows))
    if(length(zero.rows > 0)) {
      wmp$meta <- wmp$meta[-zero.rows,]
    }

    wmp$meta <- cbind(wmp$meta,rep(pat.x[cntr],nrow(wmp$meta)),
                      rep(cntr,nrow(wmp$meta)))
    Q.list[[cntr]] <- Q
    cntr <- cntr+1

    meta <- rbind(meta,wmp$meta)
  }

  mle <- optim(par=c(0,0),fn=log.lik.cov,gr=gr.covariate,
               method="BFGS",meta,Q.list,hessian=TRUE,
               control=list(trace=TRUE))

  #if(perturb) {
  #  wmp$meta[,3] <- wmp$meta[,3] + 0.01
  #}


  out <- list()
  out$cov <- mle$par
  out$se <- sqrt(diag(solve(mle$hessian)))
  return(out)
}

log.lik.cov <- function(x,meta,Q.list) {
  sum <- 0
  for(i in 1:nrow(meta)) {
    mu <- exp(x[1]+x[2]*meta[i,4])
    u <- meta[i,1]; v <- meta[i,2]
    sum <- sum+log(expm(mu*meta[i,3]*Q.list[[meta[i,5]]])[u,v])
  }
  -sum
}

gr.covariate <- function(x,meta,Q.list) {
  gr <- c(0,0)
  sum1 <- 0
  sum2 <- 0
  for(i in 1:nrow(meta)) {
    u <- meta[i,1]; v <- meta[i,2]
    mu <- exp(x[1]+x[2]*meta[i,4])

    P <- expm(meta[i,3]*mu*Q.list[[meta[i,5]]])
    num <- (Q.list[[meta[i,5]]] %*% P)[u,v]
    denom <- P[u,v]
    sum1 <- sum1 + mu*meta[i,3]*num/denom
    sum2 <- sum2+meta[i,4]*meta[i,3]*mu*num/denom
  }
  return(c(-sum1,-sum2))
}


getVu <- function(t, d, Q) {
  Mat <- matrix(0, nrow=nrow(Q), ncol=ncol(Q))
  for(i in 1:nrow(Mat)) {
    for(j in 1:ncol(Mat)) {
      if(i != j) {
        Mat[i,j] <- (exp(d[i]*t)-exp(d[j]*t))/(d[i]-d[j])
      } else{
        Mat[i,j] <- t*exp(d[i]*t)
      }
    }
  }
  return(Mat)
}

#Here alpha is a vector
gr2 <- function(alpha,meta, Q) {
  A <- Q
  ixs <- which((row(A) != col(A)) & (Q > 10^{-10}))
  A[ixs] <- alpha
  diag(A) <- 0
  diag(A) <- -rowSums(A)

  decomp <- eigen(A)
  D <- diag(decomp$values)
  X <- decomp$vectors

  r <- length(alpha)
  nabla <- rep(0, r)

  for(u in 1:r) {
    E <- matrix(0, nrow=nrow(Q), ncol=ncol(Q))
    E[which(row(A) != col(A))[u]] <- 1
    Gu <- solve(X) %*% E %*% X

    scores <- apply(meta,1,function(x) {
      Vu <- Gu*getVu(x[3], decomp$values, Q)
      Deriv <- X %*% Vu %*% solve(X)

      num <- Deriv[x[1],x[2]]
      denom <- expm(x[3]*A)[x[1],x[2]]

      x[3]*Re(num)/Re(denom)
    })
    nabla[u] <- -sum(scores)
  }
  nabla
}

fn2 <- function(alpha, meta, Q) {
  A <- Q
  ixs <- which((row(A) != col(A)) & (Q > 10^{-10}))
  A[ixs] <- alpha
  diag(A) <- 0
  diag(A) <- -rowSums(A)

  scores <- apply(meta,1,function(x) {
    u <- x[1]; v <- x[2]
    log(expm(x[3]*A)[u,v])
  })
  -sum(scores)
}

estimateMigration2 <- function(cnr,
                              pos,
                              A,
                              nrep=100,
                              perturb=FALSE,
                              dist.method="euclidean") {
  cnr <- cbind(rep(0,nrow(cnr)),cnr)
  colnames(cnr)[1] <- "PD"
  M <- as.matrix(dist(t(cnr),method=dist.method))

  Q <- A2Q(A)

  trw <- fastme.bal(M)
  trw <- root(trw,outgroup="PD",resolve.root=TRUE)
  trw <- drop.tip(trw,tip="PD")

  #trw$edge.length <- rep(1,length(trw$edge.length))
  cnr <- cnr[,-1]
  M <- M[-1,]; M <- M[,-1]
  if(is.null(colnames(cnr))) {
    labels <- 1:ncol(cnr)
  } else {
    labels <- colnames(cnr)
    true_pos <- pos
    names(true_pos) <- labels
  }
  #Prune
  wmp <- Sankoff(trw,pos,labels,A)
  min.nz <- min(wmp$meta[wmp$meta[,3] > 0,3])
  zero.rows <- which(wmp$meta[,3] < 10^{-10})
  print(length(zero.rows))
  if(perturb) {
    wmp$meta[,3] <- wmp$meta[,3] + 0.001
  } else if(length(zero.rows) > 0) {
    wmp$meta <- wmp$meta[-zero.rows,]
  }

  alpha.dim <- sum(Q > 10^{-10})
  alpha.start <- rep(wmp$ntran/sum(wmp$meta[,3]),alpha.dim) + rnorm(n=alpha.dim,sd=0.05)
  alpha.start[alpha.start < 0] <- 10^{-3}

  out <- optim(par=alpha.start,fn=fn2,gr=gr2,Q=Q,
               meta=wmp$meta,method="L-BFGS-B",lower=10^{-10},hessian=TRUE)
  cat("Convergence code ", out$convergence, "\n")
  alpha.mt <- out$par
  alpha_sd <- sqrt(diag(solve(out$hessian)))

  A <- Q
  ixs <- which((row(A) != col(A)) & (Q > 10^{-10}))
  A[ixs] <- alpha.mt
  diag(A) <- 0

  ig <- graph_from_adjacency_matrix(adjmatrix=A,weighted=TRUE)
  E(ig)$curved <- 0.2
  plot(ig, edge.width=10*E(ig)$weight, edge.label = round(E(ig)$weight,2))

  out <- list()
  out$alpha <- alpha.mt
  out$ntrans <- wmp$ntran
  out$alpha.sd <- alpha_sd
  return(out)
}

