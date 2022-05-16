

##Functions


getS <- function(pos,c1, c2,c3,c4) {
  if(pos[1] > 0) {
    if(pos[2] < c1) {
      1
    } else if(pos[2] < c2) {
      2
    } else if(pos[2] < c3) {
      3
    } else {
      4
    }
  } else {
    if(pos[2] < c1) {
      5
    } else if(pos[2] < c2) {
      6
    } else if(pos[2] < c3) {
      7
    } else {
      8
    }
  }
}


computeGD <- function(pos,D) {
  k <- 10
  #True GD
  ncell <- nrow(D)
  KNN <- t(apply(D, 1, order)[2:(k+1),])
  matches <- vapply(1:ncell,FUN.VALUE=integer(1), function(x) {
    sum(pos[x] == pos[KNN[x,]])
  })
  true.metric <- sum(matches)/(k*ncell)
  return(true.metric)
}





## Generate data
set.seed(348590)
library(SITH)
library(ggplot2)
out <- simulateTumor(max_pop = 10000,
                     death_rate=0.23,
                     driver_prob = 1,
                     selective_adv = 1.05,
                     mut_rate = 10^{-4})


#Slice
out$cell_ids <- out$cell_ids[out$cell_ids$z == 0,]


reps <- 25
val <- rep(0, reps)
for(i in 1:reps) {
  sc.sample <- randomSingleCells(tumor=out,ncells=300)
  c <- seq(min(out$cell_ids$y),max(out$cell_ids$y),length.out=4)
  c1 <- c[1]; c2 <- c[2]; c3 <- c[3]; c4 <- c[4]
  sc.S <- apply(sc.sample$positions,1,function(x) getS(x,c1,c2,c3,c4))

  D <- as.matrix(dist(sc.sample$sequencing))
  val[i] <- computeGD(sc.S,D)
}

saveRDS(val, file="sith_gdA_raw.RDS")

dfgg <- data.frame(vals=val)
library(ggplot2)
p <- ggplot(dfgg,aes(x=vals))
p <- p + geom_density(alpha=.2, fill="#FF6666")
p <- p + geom_vline(xintercept=0.175,linetype="dashed",color="blue")
p <- p + geom_vline(xintercept=0.29,linetype="dashed",color="blue")
p <- p + geom_vline(xintercept=0.31,linetype="dashed",color="blue")
p <- p + geom_vline(xintercept=0.35,linetype="dashed",color="blue")
p <- p + geom_vline(xintercept=0.36,linetype="dashed",color="blue")
p <- p + geom_vline(xintercept=0.45,linetype="dashed",color="blue")
p <- p + geom_vline(xintercept=0.625,linetype="dashed",color="blue")
p <- p + xlab("GD")+ylab("Density")
save(p, file="sith_gdA_fig.RData")
