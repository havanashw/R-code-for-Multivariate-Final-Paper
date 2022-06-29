rm(list = ls())
library(MASS)

main <- function(N=100, p=5, rho=0, r.type=1, B=1000) {
  Sig <- toeplitz(rho^seq(0, p-1))
  # Sig <- outer(rep(rho, p), rep(1, p)); diag(Sig) <- 1
  dens_func_r <- function(r, r.new, N) {
    n <- N-1
    if(r == 0) {
      dens <- gamma(n/2)*(1-r.new^2)^((n-3)/2)/(gamma((n-1)/2)*sqrt(pi))
    } else {
      temp <- 2^(n-2)*(1-r^2)^(n/2)*(1-r.new^2)^((n-3)/2)/(pi*factorial(n-2))
      dens.i <- c()
      for (i in 1:N) {
        dens.i[i] <- temp*(2*r*r.new)^(i-1)*(gamma((n+i-1)/2))^2/factorial(i-1)
      }
      dens <- sum(na.omit(dens.i))
    }
    return(dens)
  }
  
  dens_func_R <- function(R.bar, R.new, N) {
    n <- N-1
    if (R.bar == 0) {
      dens <- 2*gamma(n/2)/(gamma((p-1)/2)*gamma((n+1-p)/2))*R.new^(p-2)*(1-R.new^2)^((n+1-p)/2-1)
    } else {
      temp <- (1-R.new^2)^((n-1-p)/2)*(1-R.bar^2)^(n/2)/(gamma((n-p+1)/2)*gamma(n/2))
      dens.i <- c()
      for (i in 1:50) {
        dens.i[i] <- 2*temp*((R.bar^2)^(i-1)*R.new^(p+2*i-4)*(gamma(n/2+i-1))^2)/((factorial(i-1))*gamma((p-1)/2+i-1))
      }
      dens <- sum(na.omit(dens.i))
    }
    return(dens)
  }
  
  if(r.type == 1) {
    r12 <- c(); dens <- c(); r.vec <- seq(-1, 1, length.out=B)
    for (b in 1:B) {
      x <- mvrnorm(N, mu=rep(0, p), Sigma=Sig)
      r12[b] <- cor(x[,1], x[,2])
      r <- Sig[1,2]
      dens[b] <- dens_func_r(r, r.vec[b], N)
    }
    return(data.frame(sample_cor=r12, sample_cor_x=r.vec, d=dens, cor=r, N=N))
  }
  
  if(r.type == 2) {
    q <- 2
    Sig11 <- Sig[1:q,1:q]; Sig12 <- Sig[1:q,((q+1):p)]; Sig21 <- t(Sig12); Sig22 <- Sig[((q+1):p),((q+1):p)]
    Sig11.2 <- Sig11 - Sig12%*%solve(Sig22)%*%Sig21
    sig12.par <- Sig11.2[1,2]; sig11.par <- Sig11.2[1,1]; sig22.par <- Sig11.2[2,2];
    rho12.par <- sig12.par/sqrt(sig11.par*sig22.par)
    
    r12 <- c(); dens <- c(); r.vec <- seq(-1, 1, length.out=B)
    for (b in 1:B) {
      x <- mvrnorm(N, mu=rep(0, p), Sigma=Sig)
      A <- var(x)*(N-1)
      A11 <- A[1:q,1:q]; A12 <- A[1:q,((q+1):p)]; A21 <- t(A12); A22 <- A[((q+1):p),((q+1):p)]
      A11.2 <- A11 - A12%*%solve(A22)%*%A21
      a12.par <- A11.2[1,2]; a11.par <- A11.2[1,1]; a22.par <- A11.2[2,2]
      r12.par <- a12.par/sqrt(a11.par*a22.par)
      r12[b] <- r12.par
      
      dens[b] <- dens_func_r(rho12.par, r.vec[b], N-(p-q))
    }
    return(data.frame(sample_cor=r12, sample_cor_x=r.vec, d=dens, cor=rho12.par, N=N))
  }
  
  if(r.type == 3) {
    q <- 1
    Sig11 <- Sig[1:q,1:q]; Sig12 <- Sig[1:q,((q+1):p),drop=F]; Sig21 <- t(Sig12); Sig22 <- Sig[((q+1):p),((q+1):p)]
    R.bar <- as.numeric(sqrt((Sig12%*%solve(Sig22)%*%Sig21)/Sig11))
    
    r12 <- c(); dens <- c(); r.vec <- seq(0, 1, length.out=B)
    for (b in 1:B) {
      x <- mvrnorm(N, mu=rep(0, p), Sigma=Sig)
      A <- var(x)*(N-1)
      A11 <- A[1:q,1:q]; A12 <- A[1:q,((q+1):p),drop=F]; A21 <- t(A12); A22 <- A[((q+1):p),((q+1):p)]
      R <- as.numeric(sqrt((A12%*%solve(A22)%*%A21)/A11))
      r12[b] <- R
      dens[b] <- dens_func_R(R.bar, r.vec[b], N)
    }
    return(data.frame(sample_cor=r12, sample_cor_x=r.vec, d=dens, cor=R.bar, N=N))
  }
}

res <- list()
k <- 1
for (rho in c(0, 0.2, 0.5, 0.7)) {
  for (N in c(20, 50, 100)) {
    res[[k]] <- main(N=N, p=10, rho=rho, r.type=3, B=5000) 
    k <- k+1
  }
}

plot.res <- Reduce("rbind", res)
N.labs <- c("N=20", "N=50", "N=100")
names(N.labs) <- c("20", "50", "100")
cor.labs <- paste0("corrcoef=", unique(round(plot.res$cor, 2)))
names(cor.labs) <- paste0(unique(plot.res$cor))
library(ggplot2)
ggplot(plot.res) + 
  geom_density(aes(x=sample_cor, color="blue"), alpha=0.2, fill="blue") +
  geom_line(aes(x=sample_cor_x, y=d, color="red")) +
  geom_vline(aes(xintercept=cor), color="black", alpha=0.5,
             linetype="dashed")+
  scale_color_manual(name="", 
                     values=c("blue", "red"),
                     labels=c("Estimated", "Real")) +
  facet_grid(N ~ cor,
             labeller=labeller(N=N.labs, cor=cor.labs)) +
  theme_bw() +
  theme(axis.title=element_text(face="bold", size=14),
        axis.text.x=element_text(colour="black",
                                 size=12),
        axis.text.y=element_text(colour="black",
                                 size=12),
        legend.text=element_text(size=14),
        legend.title=element_text(face="bold", size=14)) +
  theme(panel.grid.minor=element_blank(),
        # panel.grid.major=element_blank(),
        panel.border=element_rect(size=1),
        legend.position="top") +
  theme(strip.text.x=element_text(size=14, color="black", face="bold"),
        strip.text.y=element_text(size=14, color="black", face="bold"),
        strip.background=element_rect(size=0)) +
  labs(x="Sample Correlation Coefficient", y="Density")

