library(matrixcalc)

hessian_mmvn <- function(Y,out){
  M <- out$M
  U <- out$U
  V <- out$V
  w <- out$w
  f.z.y <- out$f.z.y
  r <- dim(Y)[1]
  p <- dim(Y)[2]
  n <- dim(Y)[3]
  k <- dim(M)[3]
  if(r > 1) D_r <- duplication.matrix(r)[,-1] else D_r <- as.matrix(1)
  if(p > 1) D_p <- duplication.matrix(p) else D_p <- as.matrix(1)
  ppc <- p*r + ncol(D_r) + ncol(D_p)
  c_k <- matrix(0, ppc, k)
  C_k <- matrix(0, ppc, ppc)
  if(k > 2) A <- diag(1/w[-k]) else A <- 1/w[-k]
  if(k > 1) {
    A <- rbind(A, -rep(1/w[k], k-1))
    a_bar <- f.z.y %*% A
  }
  q <- q_t <- rep(0, k*(1+ppc)-1)
  I1 <- I2 <- Q_t <- matrix(0, k*(1+ppc)-1, k*(1+ppc)-1)
  
  for(u in 1:n){
    if(k > 1){
      q_t[1:(k-1)] <- a_bar[u,]
      Q_t[1:(k-1),1:(k-1)] <- - q_t[1:(k-1)] %*% t(q_t[1:(k-1)]) 
    }
    
    for(i in 1:k){
      V_k.inv <- ginv(V[,,i] * U[1,1,i])
      U_k.inv <- ginv(U[,,i] / U[1,1,i])
      VxU_k.inv <- ginv(V[,,i] %x% U[,,i])
      B_k.r <- U_k.inv %*% (Y[,,u] - M[,,i])
      B_k.p <- V_k.inv %*% t(Y[,,u] - M[,,i])
      R_k <- p*U_k.inv - B_k.r %*% V_k.inv %*% t(B_k.r) 
      P_k <- r*V_k.inv - B_k.p %*% U_k.inv %*% t(B_k.p)
      
      c_k[1:(p*r), i] <- VxU_k.inv %*% c(Y[,,u] - M[,,i])
      
      c_k[(p*r+1):(p*r + ncol(D_r)), i] <- - t(D_r) %*% c(R_k)/2
      
      c_k[(p*r + ncol(D_r) + 1) : (p*r + ncol(D_r) + ncol(D_p)), i] <- 
        - t(D_p) %*% c(P_k)/2 
      
      C_k[1:(p*r), 1:(p*r)] <- VxU_k.inv 
      
      C_k[1:(p*r), (p*r+1) : (p*r + ncol(D_r))] <- 
        VxU_k.inv %*% (t(B_k.r) %x% diag(r)) %*% D_r
      
      C_k[1:(p*r), (p*r + ncol(D_r) + 1) : (p*r + ncol(D_r) + ncol(D_p))] <- 
        VxU_k.inv %*% (diag(p) %x% t(B_k.p)) %*% D_p
      
      C_k[(p*r+1) : (p*r + ncol(D_r)), (p*r+1) : (p*r + ncol(D_r))] <- 
        t(D_r) %*% ((p*U_k.inv - 2 * R_k) %x% U_k.inv) %*% D_r / 2
      
      C_k[(p*r+1) : (p*r + ncol(D_r)),
          (p*r + ncol(D_r) + 1) : (p*r + ncol(D_r) + ncol(D_p))] <- 
        t(D_r) %*% ((B_k.r %*% V_k.inv) %x% (U_k.inv %*% t(B_k.p))) %*% D_p / 2
      
      C_k[(p*r + ncol(D_r) + 1) : (p*r + ncol(D_r) + ncol(D_p)),
          (p*r + ncol(D_r) + 1) : (p*r + ncol(D_r) + ncol(D_p))] <- 
        t(D_p) %*% ((r*V_k.inv - 2 * P_k) %x% V_k.inv) %*% D_p / 2
      
      C_k[lower.tri(C_k)] <- t(C_k)[lower.tri(C_k)] 
      
      if(k == 1){
        q_t <- c_k
        Q_t <- -C_k
      }else{
        
        q_t[(k+ppc*(i-1)):(k+ppc*i-1)] <- f.z.y[u,i]*c_k[,i]
        
        Q_t[1:(k-1), (k+ppc*(i-1)):(k+ppc*i-1)] <-
          f.z.y[u,i]*(A[i,] - a_bar[u,]) %*% t(c_k[,i])
        
        Q_t[(k+ppc*(i-1)):(k+ppc*i-1), (k+ppc*(i-1)):(k+ppc*i-1)] <-
          -(f.z.y[u,i]*C_k - f.z.y[u,i] * (1 - f.z.y[u,i]) * c_k[,i] %*% t(c_k[,i]))
        
        if(i > 1) {
          for(j in 1:(i-1)){
            Q_t[(k+ppc*(j-1)):(k+ppc*j-1), (k+ppc*(i-1)):(k+ppc*i-1)] <-
              -f.z.y[u,i]*f.z.y[u,j]*c_k[,j] %*% t(c_k[,i])
          }
        }
      }
    }
    q <- q + q_t
    I1 <- I1 + q_t%*%t(q_t)
    I2 <- I2 - Q_t
  }
  I2[lower.tri(I2)] <- t(I2)[lower.tri(I2)] 
  return(list(q = q,
              I1 = I1,
              I2 = I2,
              I3 = I2%*%ginv(I1)%*%(I2))
  )
}
