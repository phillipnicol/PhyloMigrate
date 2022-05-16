


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
  labels <- colnames(cnr)
  true_pos <- pos
  names(true_pos) <- labels
  #Prune
  wmp <- Sankoff(trw,pos,labels,A)
  zero.rows <- which(wmp$meta[,3] < 10^{-2})
  print(length(zero.rows))
  if(length(zero.rows > 0)) {
    wmp$meta <- wmp$meta[-zero.rows,]
  }
  out <- Newton(wmp$meta,Q,wmp$ntran/sum(trw$edge.length))
  alpha.mt <- out$alpha
  alpha_sd <- out$var


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
    if(abs(alpha.new-alpha) < 0.001) {
      alpha <- alpha.new
      var <- -1/second.deriv(meta,alpha,Q)
      break
    }
    alpha <- alpha.new
    if(alpha < 0) {
      alpha <- 0.01
    }
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

    num1 <- w*P[u,v]*(Q2 %*% P)[u,v]
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



