####################################
####                             ###
####   Multivariate functions    ###
####                             ###
####################################

# Var optimization of x variables
spls_fx <- function(X, Y, n, kx = NULL, tol = 1e-6) {
  
  # Validate kx
  if ( is.null(kx)) {
    kx <- rep((ncol(X)-1), n)
  } else {
    for (j in 1:length(kx)) {
      if (kx[j] >= ncol(X)) {
        kx[j] <- (ncol(X) -1)
      }
    }
  }
  
  # Validate n
  if( n >= nrow(X)) {
    n <- qr(X)$rank-1
  }
  
  # Center
  # meanX <- apply(X, 2, mean)
  # X <- t(apply(X, 1, function(x) x-meanX))
  #
  # meanY <- apply(Y, 2, mean)
  # Y <- t(apply(Y, 1, function(y) y-meanY))
  
  T <- NULL
  W <- NULL
  Q <- NULL
  U <- NULL
  P <- NULL
  C <- NULL
  W <- NULL
  D <- NULL
  B <- array(0, c(ncol(X), ncol(Y),n))
  
  for ( i in 1:n) {
    #pls <- pls2(X, Y, 1)
    #uo <- pls$U
    uo <- Y[,max(which(apply(Y,2,var) == max(apply(Y,2,var))))] #select only 1 column with max var in Y
    ende <- FALSE
    while (!ende) {
      wo <- t(X) %*% uo
      lambda1 <- sort(abs(wo), decreasing = TRUE)[kx[i]+1]
      zeros <- (abs(wo) <= lambda1)
      wn <- apply(wo, 2, function(g) sign(g)*(abs(g)-lambda1))
      wn[zeros] <- 0
      wn <- wn/as.vector(sqrt(t(wn)%*%wn))
      th <- as.matrix(X) %*% wn
      ch <- crossprod(Y, th)/as.vector(t(th) %*% th)
      cn <- ch
      #lambda2 <- sort(abs(ch), decreasing = TRUE)[ky[i]+1]
      #zeros <- (abs(ch) <= lambda2)
      #cn <- sapply(ch, function(g) sign(g)*(abs(g)-lambda2))
      #cn[zeros] <- 0
      #cn <- cn/as.vector(sqrt(t(cn) %*% cn))
      uhnew <- Y %*% cn
      deltau <- uhnew - uo
      unorm <- as.numeric(sqrt(t(deltau) %*% deltau))
      if (unorm < tol){
        ende <- TRUE
      }
      wo <- wn
      uo <- uhnew
    }
    ph <- t(X) %*% th /as.vector(t(th) %*% th)
    qh <- t(Y) %*% uhnew/as.vector(t(uhnew) %*% uhnew)
    dh <- t(Y) %*% th/as.vector(t(th) %*% uhnew)
    X <- X - th %*% t(ph)
    Y <- Y - th%*%t(dh)
    T <- cbind(T, th)
    Q <- cbind(Q, qh)
    U <- cbind(U, uhnew)
    P <- cbind(P, ph)
    D <- cbind(D, dh)
    C <- cbind(C, ch)
    W <- cbind(W, wn)
    B[,,i] <- W %*% solve(t(P) %*% W) %*% t(C)
    
  }
  list(P = P, T = T, Q = Q, U = U, W = W, C = C, B = B, Xf = X, Yf = Y, D = D)
}

# var optimization of x and y variables
spls_fy <- function(X, Y, n, kx = NULL, ky = NULL, tol = 1e-6, mode = "regression") {
  
  # Validate kx
  if ( is.null(kx)) {
    kx <- rep((ncol(X)-1), n)
  } else {
    for (j in 1:length(kx)) {
      if (kx[j] >= ncol(X)) {
        kx[j] <- (ncol(X) -1)
      }
    }
  }
  
  # Validate ky
  if (is.null(ky)) {
    ky <- rep((ncol(Y)-1), n)
  } else {
    for (h in 1:length(ky)) {
      if (ky[h] >= ncol(Y)) {
        ky[h] <- (ncol(Y) -1)
      }
    }
  }
  
  # Validate n
  if( n >= nrow(X)) {
    n <- qr(X)$rank-1
  }
  
  # Center
  # meanX <- apply(X, 2, mean)
  # X <- t(apply(X, 1, function(x) x-meanX))
  #
  # meanY <- apply(Y, 2, mean)
  # Y <- t(apply(Y, 1, function(y) y-meanY))
  
  T <- NULL
  W <- NULL
  Q <- NULL
  U <- NULL
  P <- NULL
  C <- NULL
  W <- NULL
  D <- NULL
  B <- array(0, c(ncol(X), ncol(Y),n))
  
  for ( i in 1:n) {
    #pls <- pls2(X, Y, 1)
    #uo <- pls$U
    uo = Y[,max(which(apply(Y,2,var) == max(apply(Y,2,var))))] #select only 1 column with max var in Y
    #uo <- Y[,1]
    ende <- FALSE
    while (!ende) {
      wo <- t(X) %*% uo
      lambda1 <- sort(abs(wo), decreasing = TRUE)[kx[i]+1]
      zeros <- (abs(wo) <= lambda1)
      wn <- sapply(wo, function(g) sign(g)*(abs(g)-lambda1))
      wn[zeros] <- 0
      wn <- wn/as.vector(sqrt(t(wn)%*%wn))
      th <- as.matrix(X) %*% wn
      ch <- crossprod(Y, th)/as.vector(t(th) %*% th)
      #ch <- ch
      lambda2 <- sort(abs(ch), decreasing = TRUE)[ky[i]+1]
      zeros <- (abs(ch) <= lambda2)
      cn <- sapply(ch, function(g) sign(g)*(abs(g)-lambda2))
      cn[zeros] <- 0
      cn <- cn/as.vector(sqrt(t(cn) %*% cn))
      uhnew <- Y %*% cn
      deltau <- uhnew - uo
      unorm <- as.numeric(sqrt(t(deltau) %*% deltau))
      if (unorm < tol){
        ende <- TRUE
      }
      wo <- wn
      uo <- uhnew
    }
    ph <- t(X) %*% th /as.vector(t(th) %*% th)
    qh <- t(Y) %*% uhnew/as.vector(t(uhnew) %*% uhnew)
    dh <- t(Y) %*% th/as.vector(t(th) %*% uhnew)
    X <- X - th %*% t(ph)
    if ( mode == "regression") {
      Y <- Y - th%*%t(dh)
    } else {
      Y <- Y - uhnew%*%t(qh)
    }
    T <- cbind(T, th)
    Q <- cbind(Q, qh)
    U <- cbind(U, uhnew)
    P <- cbind(P, ph)
    D <- cbind(D, dh)
    C <- cbind(C, ch)
    W <- cbind(W, wn)
    B[,,i] <- W %*% solve(t(P) %*% W) %*% t(C)
    
  }
  B <- W %*% solve(t(P) %*% W) %*% t(C)
  list(P = P, T = T, Q = Q, U = U, W = W, C = C, B = B, Xf = X, Yf = Y, D = D)
}

# PLS nipals algorithm
pls2 <- function (X, Y, a, it = 1000, tol = 1e-06, scale = FALSE)
{
  drop=F
  if( a >= nrow(X)) {
    a <- qr(X)$rank-1
  }
  Xh <- X
  Yh <- Y
  T <- NULL
  W <- NULL
  Q <- NULL
  U <- NULL
  P <- NULL
  D <- NULL
  C <- NULL
  W <- NULL
  for (h in 1:a) {
    perm <- 0
    nr <- 0
    uh = Yh[,max(which(apply(Yh,2,var) == max(apply(Yh,2,var))))] #select only 1 column with max var in Y
    #uh <- Yh[, 1]
    ende <- FALSE
    while (!ende) {
      nr <- nr + 1
      wh <- t(Xh) %*% uh
      wh <- wh/as.vector(sqrt(t(wh) %*% wh))
      th <- Xh %*% wh
      ch <- crossprod(Yh, th)/drop(crossprod(th))
      uhnew <- Yh %*% ch
      deltau <- uhnew - uh
      unorm <- as.numeric(sqrt(t(deltau) %*% deltau))
      if (unorm < tol) {
        ende <- TRUE
      }
      uh <- uhnew
    }
    ph <- t(Xh) %*% th/as.vector(t(th) %*% th)
    qh <- t(Yh) %*% uh/as.vector(t(uh) %*% uh)
    dh <- t(uh) %*% th/as.vector(t(th) %*% th)
    Xh <- Xh - th %*% t(ph)
    Yh <- Yh - (th %*% t(ch)) * as.vector(dh)
    T <- cbind(T, th)
    Q <- cbind(Q, qh)
    U <- cbind(U, uh)
    P <- cbind(P, ph)
    D <- c(D, dh)
    C <- cbind(C, ch)
    W <- cbind(W, wh)
    B <- W %*% solve(t(P) %*% W) %*% t(C)
  }
  list(P = P, T = T, Q = Q, U = U, D = D, W = W, C = C, B = B)
}


# MB PLS algorithm
mbpls1 <- function(X, Y, a, tol = 1e-6) {
  
  design <- matrix(c(1,0,0,1),2,2)
  
  # Fold X
  Tb <- list("1" = matrix(0, nrow(X[[1]]), a))
  Wb <- list("1" = matrix(0, ncol(X[[1]]), a))
  Pb <- list("1" = matrix(0, ncol(X[[1]]), a))
  Bb <- list("1" = array(0, c(ncol(X[[1]]), ncol(Y), a)))
  for ( i in 2:length(X)) {
    Tb[[toString(i)]] <- matrix(0, nrow(X[[i]]), a)
    Wb[[toString(i)]] <- matrix(0, ncol(X[[i]]), a)
    Pb[[toString(i)]] <- matrix(0, ncol(X[[i]]), a)
    Bb[[toString(i)]] <- array(0, c(ncol(X[[i]]), ncol(Y), a))
  }
  W <- NULL
  T <- matrix(0, nrow(Y), (length(X)))
  Tsup <- NULL
  C <- NULL
  U <- NULL
  
  for (h in 1:a) {
    #tsup <- as.matrix(eigen(Xsup%*%t(Xsup))$vector[,1])
    #tsup <- tsup * as.vector(1/sqrt(t(tsup)%*%tsup))
    #print(tsup)
    u <- Y[,1]
    ende <- FALSE
    while (!ende)  {
      for (s in 1:length(X)) {
        Wb[[s]][,h] <- t(X[[s]]) %*% u / as.vector(t(u)%*%u)
        Wb[[s]][,h] <- Wb[[s]][,h] / as.vector(sqrt(t(Wb[[s]][,h])%*%Wb[[s]][,h]))
        Tb[[s]][,h] <- (X[[s]] %*% Wb[[s]][,h]) / sqrt(ncol(X[[s]]))
        T[,s] <- Tb[[s]][,h]
      } 
      #Z <- matrix(0, nrow(X[[s]]), length(X))
      #for (s in length(X)) {
      #  Z[,s] <- rowSums(mapply("*", design[s,], as.data.frame(T)))
      #  Wb[[s]][,h] <- t(X[[s]]) %*% Z[,s]
      #  Wb[[s]][,h] <- Wb[[s]][,h] / as.vector(sqrt(t(Wb[[s]][,h])%*%Wb[[s]][,h]))
      #  Tb[[s]][,h] <- (X[[s]] %*% Wb[[s]][,h])# / sqrt(ncol(X[[s]]))
      #  T[,s] <- Tb[[s]][,h]
      #}
      wh <- t(T) %*% u / as.vector(t(u)%*%u)
      wh <- wh / as.vector(sqrt(t(wh)%*%wh))
      tsup <- T%*%wh / as.vector(t(wh)%*%wh)
      ch <- t(Y) %*% tsup / as.vector(t(tsup)%*%tsup)
      unew <- Y%*%ch/as.vector(t(ch)%*%ch)
      deltau <- unew - u
      unorm <- as.numeric(sqrt(t(deltau) %*% deltau))
      if (unorm <= (tol)) {
        ende <- TRUE
      }
      u <- unew
    }
    dh <- t(u) %*% tsup/as.vector(t(tsup) %*% tsup)
    #print(dh)
    W <- cbind(W, wh)
    Tsup <- cbind(Tsup, tsup)
    C <- cbind(C, ch)
    U <- cbind(U, u)
    for (s in 1:length(X)) {
      Pb[[s]][,h] <- t(X[[s]]) %*% as.matrix(Tb[[s]][,h]) / as.vector(t(Tb[[s]][,h])%*%Tb[[s]][,h])
      X[[s]] <- X[[s]] - Tb[[s]][,h]%*%t(Pb[[s]][,h])
      #print(dim(Wb[[s]][,1:h] %*% solve(t(Pb[[s]][,1:h])%*%Wb[[s]][,1:h]) %*% t(C)))
      Bb[[s]][,,h] <- Wb[[s]][,1:h] %*% solve(t(Pb[[s]][,1:h])%*%Wb[[s]][,1:h]) %*% t(C)
    }
    Y <- Y - (tsup%*%t(ch) * as.vector(dh))
  }
  D <- t(U) %*% Tsup / as.vector(t(Tsup)%*%Tsup)
  list(W = W, Tsup = Tsup, C = C, U = U, Tb = Tb, Wb = Wb, Pb = Pb, B = Bb, D = D)
}




####################################
####                             ###
####    Tune sparse functions    ###
####                             ###
####################################

## S-PLS
tune.comp <- function(X, Y, factor, fold, ncomp = 10, rep = 50, low.limit = 0.5, mode = "regression",
                      optimal.threshold = 0.0975,  option = "Q") {
 
  
  # Test mode
  if (!is.element(mode, c("regression", "canonical"))){
    stop("Not valid mode: regression or canonical")
  }
  
  if (nrow(X) != nrow(Y)) {
    stop("nrow of X and Y must be equal")
  }
  # Centrar
  meanX <- apply(X, 2, mean)
  X <- t(apply(X, 1, function(x) x-meanX))

  meanY <- apply(Y, 2, mean)
  Y <- t(apply(Y, 1, function(y) y-meanY))
  
  # Resultado total
  q2t <- matrix(0, ncomp, rep)
  r2t <- matrix(0, ncomp, rep)
  
  # See fold suitability
  if (nrow(X) <= 5) {
    message("* Low number of replicates. Leave-1-out cv will be performed.")
  } 
  # if ((nlevels(factor)/length(factor))==1 & option == "Q") {
  #   message("* Only 1 sample for each condition is not suitable for prediction (Q^2 optimization).")
  # }
  
  message("\n Running optimization:")
  pb <- txtProgressBar(style = 3)
  
  if (option == "Q") {
    for (r in 1:rep){
      #message(paste("\t- Permutation:", r))
      if (nrow(X) <= 5) {
        gr <- create.fold(factor, length(factor))
      }
      if ((nlevels(factor)/length(factor))==1) {
        gr <- create.fold(factor, length(factor))
      } else {
        gr <- create.fold(factor, fold)
      }
      # Opt. nº var(X)
      var.x <- c()
      q2_rep <- c()
      for ( n in 1:ncomp) {
        Y.hat <- matrix(0, nrow(Y), ncol(Y))
        Y.model <- matrix(0, nrow(Y), ncol(Y))
        for ( s in 1:length(gr)) {
          pls <- pls2(X[-gr[[s]],,drop=F], as.matrix(Y[-gr[[s]],,drop=F]), a=n)
          Bt <- pls$B
          Y.hat[gr[[s]],] <- X[gr[[s]],,drop=F] %*% Bt
        }
        # Q2 optimization
        E <- Y - Y.hat
        press <- sum(apply(E,2,function(x) sum(x**2)))
        sct <- sum(apply(Y,2,function(x) sum(x**2)))
        q2_rep <- c(q2_rep, (1-(press/sct)))
        setTxtProgressBar(pb, ((n/ncomp)*(r/rep)+((r-1)/rep)))
      }
      q2t [,r] <- q2_rep
    }
    # output final optimal 
    setTxtProgressBar(pb, 1)
    # par(mfrow=c(1,2))
    opt.comp <- optimal.comp(q2t, "opt", p = "Q^2", low.limit = low.limit,
                             optimal.threshold = optimal.threshold, plot = TRUE)
    Q2 = q2t

  } else if (option == "R") {
    for (r in 1:rep){
      #message(paste("\t- Permutation:", r))
      if (nrow(X) <= 5) {
        gr <- create.fold(factor, length(factor))
      }
      if ((nlevels(factor)/length(factor))==1) {
        gr <- create.fold(factor, length(factor))
      } else {
        gr <- create.fold(factor, fold)
      }
      # Opt. nº var(X)
      var.x <- c()
      r2_rep <- matrix(0, ncomp, length(gr))
      for ( n in 1:ncomp) {
        Y.model <- matrix(0, nrow(Y), ncol(Y))
        for ( s in 1:length(gr)) {
          pls <- pls2(X[-gr[[s]],,drop=F], as.matrix(Y[-gr[[s]],,drop=F]), a=n)
          Y.model <- pls$T%*%t(pls$C)
          E <- as.matrix(Y[-gr[[s]],]) - Y.model
          scr <- sum(apply(E,2,function(x) sum(x**2)))
          sct <- sum(apply(as.matrix(Y[-gr[[s]],]),2,function(x) sum(x**2)))
          r2 <- 1 - (scr/sct)
          r2_rep[n,s] <- r2
        }
        setTxtProgressBar(pb, max(((n/ncomp)*(r/rep)+((r-1)/rep)),0.9))
      }
      r2t [,r] <- apply(r2_rep, 1, mean)
    }
    # output final optimal 
    setTxtProgressBar(pb, 1)
    # par(mfrow=c(1,2))
    opt.comp <- optimal.comp(r2t, "opt", p = "R^2", low.limit = low.limit,
                             optimal.threshold = optimal.threshold, plot = TRUE)
    names(opt.comp)[3] <- "R2"
    R2 = r2t
    #dev.off()
  }
  else {
    print("option must be Q or R")
  }
  #list(Q2 = q2t, R2 = r2t)
  return(opt.comp)
}


tune.var <- function(X, Y, varX, varY=NULL, factor, fold, ncomp = 10, option = "Q", rep = 50, 
                     low.limit = 0.5, sel.var = "opt", mode = "regression",
                     optimal.threshold = 0.0975) {
  # Test mode
  if (!is.element(mode, c("regression", "canonical"))){
    stop("Not valid mode: regression or canonical")
  }
  
  # Test rows
  if (nrow(X) != nrow(Y)) {
    stop("nrow of X and Y must be equal")
  }
  
  # Centrar
  # meanX <- apply(X, 2, mean)
  # X <- t(apply(X, 1, function(x) x-meanX))
  # 
  # meanY <- apply(Y, 2, mean)
  # Y <- t(apply(Y, 1, function(y) y-meanY))
  
  # Resultado total
  q2t <- list("1" = array(0, c(ncomp, length(varX), max(2,length(varY)))))
  r2t <- list("1" = array(0, c(ncomp, length(varX), max(2,length(varY)))))
  
  # Ordenar numero de Var
  org.var <- varX
  varX <- sort(varX, decreasing = F)
  
  # See fold suitability
  if (nrow(X) <= 5) {
    message("* Low number of replicates. Leave-1-out cv will be performed.")
  } 
  if ((nlevels(factor)/length(factor))==1 & option == "Q") {
    message("* Only 1 sample for each condition is not suitable for prediction (Q^2 optimization).")
  }
  message("\n Running optimization: ")
  pb <- txtProgressBar(style = 3)
  
  if (option == "Q") {
    for (r in 1:rep){
      # Opt. nº var(X)
      
      if (nrow(X) <= 5) {
        gr <- create.fold(factor, length(factor))
      } 
      if ((nlevels(factor)/length(factor))==1) {
        gr <- create.fold(factor, length(factor))
      } else {
        gr <- create.fold(factor, fold)
      }
      var.x <- c()
      var.y <- c()
      q2_rep <- matrix(0, ncomp, length(varX))
      q2y_rep <- array(0, c(ncomp, length(varX), max(2,length(varY))))
      for ( n in 1:ncomp) {
        nvar <- 1
        var.yforx <- c()
        for ( v in 1:length(varX)) {
          if (!is.null(varY)) {
            for ( q in 1:length(varY)) {
              Y.hat <- matrix(0, nrow(Y), ncol(Y))
              for ( s1 in 1:length(gr)) {
                s.pls1 <- spls_fy(X[-gr[[s1]],], as.matrix(Y[-gr[[s1]],]), n=n, kx = c(var.x, varX[v]), ky = c(var.y, varY[q]), mode = mode)
                Bt <- s.pls1$B
                Y.hat[gr[[s1]],] <- X[gr[[s1]],] %*% Bt
              }
              # Q2 optimization
              E <- Y - Y.hat
              #rmsep <- c(rmsep, sqrt(sum(apply(E,2,function(x) sum(x**2)))))
              press <- sum(apply(E,1,function(x) sum(x**2)))
              sct <- sum(apply(Y,2,function(x) sum(x**2)))
              q2y_rep[n,v,q] <- (1 - (press/sct))
            }
            # Var opt.
            var.yforx <- c(var.yforx, optimal.var.n(q2y_rep[n,v,],option=sel.var)$var)
          } else  {
            var.y <- rep(ncol(Y), n)
            Y.hat <- matrix(0, nrow(Y), ncol(Y))
            for ( s in 1:length(gr)) {
              s.pls <- spls_fx(X[-gr[[s]],], as.matrix(Y[-gr[[s]],]), n=n, kx = c(var.x, varX[v]))
              Bt <- s.pls$B[,,n]
              Y.hat[gr[[s]],] <- X[gr[[s]],] %*% Bt
            }
            # Q2 optimization
            E <- Y - Y.hat
            #rmsep <- c(rmsep, sqrt(sum(apply(E,2,function(x) sum(x**2)))))
            press <- sum(apply(E,1,function(x) sum(x**2)))
            sct <- sum(apply(Y,2,function(x) sum(x**2)))

            q2y_rep[n,v,] <- rep((1 - (press/sct)), 2)
          }
          var.yforx <- c(var.yforx, optimal.var.n(q2y_rep[n,v,], option=sel.var,
                                                  optimal.threshold)$var)
        }
        # Var opt.
        var.x <- c(var.x, optimal.var.nx(q2y_rep[n,,], var.yforx, option=sel.var,
                                         optimal.threshold)$var)
        var.y <- c(var.y, var.yforx[var.x])
        setTxtProgressBar(pb, min(1,((n/ncomp)*(r/rep)+((r-1)/rep))))
      }
      q2t[[toString(r)]] <- q2y_rep
    }
    # output final optimal 
    if ( is.null(varY)) {
      opt.var <- optimal.var.null(q2t, sel.var, p = "Q^2", varX, varY, low.limit = low.limit,
                                  optimal.threshold)
    } else {
      opt.var <- optimal.var(q2t, sel.var, p = "Q^2", varX, varY, low.limit = low.limit,
                             optimal.threshold)
    }
    
  } else {
    for (r in 1:rep){
      # Opt. nº var(X)
      
      if (nrow(X) <= 5) {
        gr <- create.fold(factor, length(factor))
      } 
      if ((nlevels(factor)/length(factor))==1) {
        gr <- create.fold(factor, length(factor))
      } else {
        gr <- create.fold(factor, fold)
      }
      var.x <- c()
      var.y <- c()
      #q2_rep <- matrix(0, ncomp, length(varX))
      #q2y_rep <- array(0, c(ncomp, length(varX), max(2,length(varY))))
      r2y_rep <- array(0, c(ncomp, length(varX), max(2,length(varY))))
      for ( n in 1:ncomp) {
        nvar <- 1
        var.yforx <- c()
        for ( v in 1:length(varX)) {
          if (!is.null(varY)) {
            for ( q in 1:length(varY)) { 
              #Y.hat <- matrix(0, nrow(Y), ncol(Y))
              r2_parcial <- c()
              for ( s1 in 1:length(gr)) {
                s.pls1 <- spls_fy(X[-gr[[s1]],], as.matrix(Y[-gr[[s1]],]), n=n, kx = c(var.x, varX[v]), ky = c(var.y, varY[q]), mode = mode)
                #Bt <- s.pls1$B[,,n]
                #Y.hat[gr[[s1]],] <- X[gr[[s1]],] %*% Bt
                Y.model <- s.pls1$T%*%t(s.pls1$C)
                E <- as.matrix(Y[-gr[[s1]],]) - Y.model
                scr <- sum(apply(E,2,function(x) sum(x**2)))
                sct <- sum(apply(as.matrix(Y[-gr[[s1]],]),2,function(x) sum(x**2)))
                r2 <- 1 - (scr/sct)
                r2_parcial <- c(r2_parcial, r2)
              }
              r2y_rep[n,v,q] <- mean(r2_parcial)
            }
            # Var opt.
            var.yforx <- c(var.yforx, optimal.var.n(r2y_rep[n,v,], sel.var,
                                                    optimal.threshold)$var[1])
          } else  {
            var.y <- rep(ncol(Y), n)
            r2_parcial <- c()
            #Y.hat <- matrix(0, nrow(Y), ncol(Y))
            for ( s in 1:length(gr)) {
              s.pls <- spls_fx(X[-gr[[s]],], as.matrix(Y[-gr[[s]],]), n=n, kx = c(var.x, varX[v]))
              Y.model <- s.pls$T%*%t(s.pls$C)
              E <- as.matrix(Y[-gr[[s]],]) - Y.model
              scr <- sum(apply(E,2,function(x) sum(x**2)))
              sct <- sum(apply(as.matrix(Y[-gr[[s]],]),2,function(x) sum(x**2)))
              r2 <- 1 - (scr/sct)
              r2_parcial <- c(r2_parcial, r2)
            }
            r2y_rep[n,v,] <- rep(mean(r2_parcial), 2)
          }
          var.yforx <- c(var.yforx, optimal.var.n(r2y_rep[n,v,], sel.var,
                                                  optimal.threshold)$var)
        }
        # Var opt.
        var.x <- c(var.x, optimal.var.nx(r2y_rep[n,,], var.yforx, sel.var,
                                         optimal.threshold)$var)
        var.y <- c(var.y, var.yforx[var.x])
        setTxtProgressBar(pb, min(1,((n/ncomp)*(r/rep)+((r-1)/rep))))
      }
      r2t[[toString(r)]] <- r2y_rep
    }
    # output final optimal 
    if ( is.null(varY)) {
      opt.var <- optimal.var.null(r2t, sel.var, p = "R^2", varX, varY, low.limit = low.limit,
                                  optimal.threshold)
    } else {
      opt.var <- optimal.var(r2t, sel.var, p = "R^2", varX, varY, low.limit = low.limit,
                             optimal.threshold)
    }
    
  }
  setTxtProgressBar(pb, 1)
  if ( is.null(varY)) {
    resopt <- varX[opt.var$optvar]
    resmax <- varX[opt.var$maxvar]
    names <- c()
    for ( i in 1:ncomp) {
      names <- c(names, paste("Comp.", i))
    }
    names(resopt) <- names
    names(resmax) <- names
    list(var.opt.x = resopt, var.max.x = resmax, validation.data = opt.var$rq)
  } else {
    resoptx <- varX[opt.var$optvarx]
    resmaxx <- varX[opt.var$maxvarx]
    resopty <- varY[opt.var$optvary]
    resmaxy <- varY[opt.var$maxvary]
    names <- c()
    for ( i in 1:ncomp) {
      names <- c(names, paste("Comp.", i))
    }
    names(resoptx) <- names
    names(resmaxx) <- names
    names(resopty) <- names
    names(resmaxy) <- names
    list(var.opt.x = resoptx, var.opt.y = resopty, 
         var.max.x = resmaxx, var.max.y = resmaxy, 
         validation.data = opt.var$rq)
  }
}


## MB-S-PLS
tune.comp.block <- function(X, Y, design, factor, fold=2, ncomp = 5, rep = 5, low.limit = 0.5, mode = "regression",
                            optimal.threshold = 0.0975) {
  option = "Q"
  
  # Test mode
  if (!is.element(mode, c("regression", "canonical"))){
    stop("Not valid mode: regression or canonical")
  }
  
  # Resultado total
  q2t <- matrix(0, ncomp, rep)

  # See fold suitability
  if (nrow(X[[1]]) <= 5) {
    message("* Low number of replicates. Leave-1-out cv will be performed.")
  } 
  if ((nlevels(factor)/length(factor))==1 & option == "Q") {
    message("* Only 1 sample for each condition is not suitable for prediction (Q^2 optimization).")
  }
  
  message("\n Running optimization:")
  pb <- txtProgressBar(style = 3)
  
  for (r in 1:rep){
    #message(paste("\t- Permutation:", r))
    if (nrow(X[[1]]) <= 5) {
      gr <- create.fold(factor, length(factor))
    } 
    if (fold==1) {
      gr <- sapply(1:nrow(X[[1]]),list)
    }
    #if ((nlevels(factor)/length(factor))==1) {
    #  gr <- create.fold(factor, length(factor))
    #} else {
    #  gr <- create.fold(factor, fold)
    #}
    # Opt. nº var(X)
    var.x <- c()
    q2_rep <- c()
    Y.hat <- array(0, c(nrow(Y), ncol(Y), ncomp))
    for ( s in 1:length(gr)) {
      Xh <- X
      Xpred <- X
      for (i in 1:length(X)) {
        Xh[[i]] <- X[[i]][-gr[[s]],]
        Xpred[[i]] <- t(as.matrix(X[[i]][gr[[s]],]))
      }
      #print(as.matrix(Y[-gr[[s]],]))
      Y.p <- as.matrix(Y[-gr[[s]],]); rownames(Y.p) <- rownames(Xh[[1]])
      #print(Y.p)
      plsda <- suppressMessages(block.pls(Xh, Y.p, ncomp = ncomp, design=design, scale=F, mode = mode))
      Bt <- predict(plsda, newdata = Xpred)
      for ( n in 1:ncomp) {
        Y.hat[gr[[s]],,n] <- Bt$AveragedPredict[,,n]
        
        l <- ((s/length(gr))*(r/rep)+((r-1)/rep))
        setTxtProgressBar(pb, (n*l/ncomp))
      }
      setTxtProgressBar(pb, ((s/length(gr))*(r/rep)+((r-1)/rep)))
    }
    for ( n in 1:ncomp) {
      # Q2 optimization
      E <- Y - Y.hat[,,n]
      #rmsep <- c(rmsep, sqrt(sum(apply(E,2,function(x) sum(x**2)))))
      press <- sum(apply(E,2,function(x) sum(x**2)))
      sct <- sum(apply(Y,2,function(x) sum(x**2)))
      q2_rep <- c(q2_rep, (1-(press/sct)))
    }
    q2t [,r] <- q2_rep
  }
  
  # Cumulative R2
  Xpred = X
  Ypred = Y
  R2 = list()
  for (i in 1:ncomp){
    Y.hat <- predict(plsda, Xpred)$AveragedPredict[,,i]
    E <- Ypred - Y.hat
    scr <- sum(t(E)%*%E)
    sct <- sum(t(Y)%*%Y)
    R2[i] = (1 - (scr/sct))
  }
  return(list(R2=R2, Q2=q2t))
  
  # output final optimal 
  setTxtProgressBar(pb, 1)
  png("Q2.png")
  opt.comp <- optimal.comp(q2t, "opt", p = "Q^2", low.limit = low.limit,
                           optimal.threshold, plot = TRUE)
  opt.comp
  dev.off()
  
  png("Cumulative R2.png")
  barplot(unlist(R2), ylim=c(0,1), main = "Cumulative R2", names.arg = c(1:ncomp), xlab = "Nº components", col = "light blue")
  dev.off()
  
  
  
}

tune.var.block <- function(X, Y, varX, varY=NULL, factor, fold, ncomp = 10, option = "Q", rep = 50, 
                           low.limit = 0.5, sel.var = "opt", design, mode = "regression",
                           optimal.threshold = 0.0975) {
  # Test mode
  if (!is.element(mode, c("regression", "canonical"))){
    stop("Not valid mode: regression or canonical")
  }
  
  # Test rows
  if (nrow(X[[1]]) != nrow(Y)) {
    stop("nrow of X and Y must be equal")
  }
  
  if (length(X) != length(varX)) {
    stop("Nº of elements in X is not equal to nº of elements in varX")
  }
  
  # Include total number of Xs variables
  for ( l in 1:length(varX)) {
    if (!is.element(ncol(X[[l]]), varX[[l]])) {
      varX[[l]] <- c(varX[[l]], ncol(X[[l]]))
    }
  }
  
  # Create grid of combinations
  if (is.null(varY)) {
    grid <- expand.grid(varX)[order(expand.grid(varX)[,1]),]
  } else {
    varib <- varX; varib[["Y"]] <- varY
    grid <- expand.grid(varib)[order(expand.grid(varib)[,1]),]
  }
  
  
  # Resultado total
  q2t <- array(0, c(ncomp, dim(grid)[1], rep))
  r2t <- array(0, c(ncomp, dim(grid)[1], rep))
  
  #varX <- sort(varX, decreasing = F)
  
  # See fold suitability
  if (nrow(X[[1]]) <= 5) {
    message("* Low number of replicates. Leave-1-out cv will be performed.")
  } 
  if ((nlevels(factor)/length(factor))==1 & option == "Q") {
    message("* Only 1 sample for each condition is not suitable for prediction (Q^2 optimization).")
  }
  message("\n Running optimization: ")
  #pb <- txtProgressBar(style = 3)
  
  if (option == "Q") {
    for (r in 1:rep){
      # Opt. nº var(X)
      
      if (nrow(X[[1]]) <= 5) {
        gr <- create.fold(factor, length(factor))
      } 
      if ((nlevels(factor)/length(factor))==1) {
        gr <- create.fold(factor, length(factor))
      } else {
        gr <- create.fold(factor, fold)
      }
      var.list.n <- sapply(grid[1,],function(x)c(), simplify = F)
      q2_rep <- matrix(0, ncomp, length(varX))
      q2y_rep <- array(0, c(ncomp, length(varX), max(2,length(varY))))
      print("i1")
      for ( n in 1:ncomp) {
        nvar <- 1
        var.yforx <- c()
        Y.hat <- matrix(0, nrow(Y), ncol(Y))
        cont <- 0
        if ( is.null(varY)) {
          q2_per_n <- apply(grid, 1, function(gd) {
            var.list <- sapply(1:length(gd),function(x){c(var.list.n[[x]],gd[x])}, simplify = T)
            names(var.list) <- names(X)
            #var.y <- rep(ncol(Y), n)
            #Y.hat <- matrix(0, nrow(Y), ncol(Y))
            for ( s in 1:length(gr)) {
              Xh <- X
              Xpred <- X
              for (i in 1:length(X)) {
                Xh[[i]] <- X[[i]][-gr[[s]],]
                Xpred[[i]] <- X[[i]][gr[[s]],]
              }
              print("i2")
              x.list <- sapply(var.list[1:length(varX)], c, simplify = F)
              names(x.list) <- names(varX)
              spls <- suppressMessages(block.spls(Xh, as.matrix(Y[-gr[[s]],]), ncomp = n, 
                                                  design=design, scale=F,
                                                  keepX = x.list,
                                                  mode = mode))
              print("i3")
              Bt <- predict(spls, newdata = Xpred)
              Y.hat[gr[[s]],] <- Bt$WeightedPredict[,,n]
              print("i4")
            }
            print("i5")
            # Q2 optimization
            E <- Y - Y.hat
            #rmsep <- c(rmsep, sqrt(sum(apply(E,2,function(x) sum(x**2)))))
            press <- sum(apply(E,1,function(x) sum(x**2)))
            sct <- sum(apply(Y,2,function(x) sum(x**2)))
            (1 - (press/sct))
            #}
            #}
          }
          )
          print("i6")
          # Var opt.
          opt.per.n <- optimal.var.n.null.block(q2_per_n, sel.var,
                                                optimal.threshold)$var
          var.list.n <- sapply(names(var.list.n),function(x){
            c(var.list.n[[x]],
              sapply(grid[opt.per.n,],function(y)y, simplify = F)[[x]])
          }, simplify = F)
          q2t[n,,r] <- q2_per_n
        } else {
          q2_per_n <- apply(grid, 1, function(gd) {
            var.list <- sapply(1:length(gd),function(x){c(var.list.n[[x]],gd[x])}, simplify = T)
            names(var.list) <- names(X)
            #var.y <- rep(ncol(Y), n)
            #Y.hat <- matrix(0, nrow(Y), ncol(Y))
            for ( s in 1:length(gr)) {
              Xh <- X
              Xpred <- X
              for (i in 1:length(X)) {
                Xh[[i]] <- X[[i]][-gr[[s]],]
                Xpred[[i]] <- X[[i]][gr[[s]],]
              }
              x.list <- sapply(var.list[1:length(varX)], c, simplify = F)
              names(x.list) <- names(varX)
              spls <- suppressMessages(block.spls(Xh, as.matrix(Y[-gr[[s]],]), ncomp = n, 
                                                  design=design, scale=F,
                                                  keepX = x.list,
                                                  keepY = var.list[length(var.list)],
                                                  mode = mode))
              Bt <- predict(spls, newdata = Xpred)
              Y.hat[gr[[s]],] <- Bt$WeightedPredict[,,n]
            }
            # Q2 optimization
            E <- Y - Y.hat
            #rmsep <- c(rmsep, sqrt(sum(apply(E,2,function(x) sum(x**2)))))
            press <- sum(apply(E,1,function(x) sum(x**2)))
            sct <- sum(apply(Y,2,function(x) sum(x**2)))
            (1 - (press/sct))
            #}
            #}
          }
          )
          # Var opt.
          opt.per.n <- optimal.var.n.block(q2_per_n, sel.var,
                                           optimal.threshold = optimal.threshold)$var
          var.list.n <- sapply(names(var.list.n),function(x){
            c(var.list.n[[x]],
              sapply(grid[opt.per.n,],function(y)y, simplify = F)[[x]])
          }, simplify = F)
          q2t[n,,r] <- q2_per_n
        }
      }
      #     q2t[[toString(r)]] <- q2y_rep
    }
    # output final optimal 
    if ( is.null(varY)) {
      opt.var <- optimal.var.null.block(q2t, sel.var, p = "Q^2", grid, varX, low.limit = low.limit,
                                        optimal.threshold)
    } else {
      opt.var <- optimal.var.block(q2t, sel.var, p = "Q^2", grid, varX, varY, low.limit = low.limit,
                                   optimal.threshold)
    }
    
  } else {
    for (r in 1:rep){
      # Opt. nº var(X)
      
      if (nrow(X[[1]]) <= 5) {
        gr <- create.fold(factor, length(factor))
      } 
      if ((nlevels(factor)/length(factor))==1) {
        gr <- create.fold(factor, length(factor))
      } else {
        gr <- create.fold(factor, fold)
      }
      var.list.n <- sapply(grid[1,],function(x)c(), simplify = F)
      #q2_rep <- matrix(0, ncomp, length(varX))
      #q2y_rep <- array(0, c(ncomp, length(varX), max(2,length(varY))))
      for ( n in 1:ncomp) {
        nvar <- 1
        var.yforx <- c()
        Y.hat <- matrix(0, nrow(Y), ncol(Y))
        cont <- 0
        if ( is.null(varY)) {
          r2_per_n <- apply(grid, 1, function(gd) {
            var.list <- sapply(1:length(gd),function(x){c(var.list.n[[x]],gd[x])}, simplify = F)
            names(var.list) <- names(X)
            #var.y <- rep(ncol(Y), n)
            #Y.hat <- matrix(0, nrow(Y), ncol(Y))
            r2_parcial <- c()
            for ( s in 1:length(gr)) {
              Xh <- X
              Xpred <- X
              for (i in 1:length(X)) {
                Xh[[i]] <- X[[i]][-gr[[s]],]
                Xpred[[i]] <- X[[i]][gr[[s]],]
              }
              x.list <- sapply(var.list[1:length(varX)], c, simplify = F)
              names(x.list) <- names(varX) 
              spls <- suppressMessages(block.spls(Xh, as.matrix(Y[-gr[[s]],]), ncomp = n, 
                                                  design=design, scale=F,
                                                  keepX = x.list,
                                                  mode = mode))
              Y.model <- predict(spls, newdata = Xh)$WeightedPredict[,,n]
              E <- as.matrix(Y[-gr[[s]],]) - Y.model
              scr <- sum(apply(E,2,function(x) sum(x**2)))
              sct <- sum(apply(as.matrix(Y[-gr[[s]],]),2,function(x) sum(x**2)))
              r2 <- 1 - (scr/sct)
              r2_parcial <- c(r2_parcial, r2)
            }
            mean(r2_parcial)
          }
          )
          # Var opt.
          opt.per.n <- optimal.var.n.block(r2_per_n, sel.var,
                                           optimal.threshold = optimal.threshold)$var
          var.list.n <- sapply(names(var.list.n),function(x){
            c(var.list.n[[x]],
              sapply(grid[opt.per.n,],function(y)y, simplify = F)[[x]])
          }, simplify = F)
          r2t[n,,r] <- r2_per_n
        } else {
          r2_per_n <- apply(grid, 1, function(gd) {
            var.list <- sapply(1:length(gd),function(x){c(var.list.n[[x]],gd[x])}, simplify = T)
            names(var.list) <- names(X)
            #var.y <- rep(ncol(Y), n)
            #Y.hat <- matrix(0, nrow(Y), ncol(Y))
            r2_parcial <- c()
            for ( s in 1:length(gr)) {
              Xh <- X
              Xpred <- X
              for (i in 1:length(X)) {
                Xh[[i]] <- X[[i]][-gr[[s]],]
                Xpred[[i]] <- X[[i]][gr[[s]],]
              }
              x.list <- sapply(var.list[1:length(varX)], c, simplify = F)
              names(x.list) <- names(varX) 
              spls <- suppressMessages(block.spls(Xh, as.matrix(Y[-gr[[s]],]), ncomp = n, 
                                                  design=design, scale=F,
                                                  keepX = x.list,
                                                  keepY = var.list[length(var.list)],
                                                  mode = mode))
              Y.model <- predict(spls, newdata = Xh)$WeightedPredict[,,n]
              E <- as.matrix(Y[-gr[[s]],]) - Y.model
              scr <- sum(apply(E,2,function(x) sum(x**2)))
              sct <- sum(apply(as.matrix(Y[-gr[[s]],]),2,function(x) sum(x**2)))
              r2 <- 1 - (scr/sct)
              r2_parcial <- c(r2_parcial, r2)
            }
            mean(r2_parcial)
          }
          )
          # Var opt.
          opt.per.n <- optimal.var.n.block(r2_per_n, sel.var,
                                           optimal.threshold = optimal.threshold)$var
          var.list.n <- sapply(names(var.list.n),function(x){
            c(var.list.n[[x]],
              sapply(grid[opt.per.n,],function(y)y, simplify = F)[[x]])
          }, simplify = F)
          r2t[n,,r] <- r2_per_n
        }
      }
      #     q2t[[toString(r)]] <- q2y_rep
    }
    # output final optimal 
    if (is.null(varY)) {
      opt.var <- optimal.var.null.block(r2t, sel.var, p = "R^2", grid, varX, low.limit = low.limit,
                                        optimal.threshold)
    } else {
      opt.var <- optimal.var.block(r2t, sel.var, p = "R^2", grid, varX, varY, low.limit = low.limit,
                                   optimal.threshold)
    }
    # 
  }
  max <- sapply(grid[opt.var$max.var[1],],function(y)y, simplify = F)
  opt <- sapply(grid[opt.var$opt.var[1],],function(y)y, simplify = F)
  
  for (i in 2:ncomp) {
    max <- sapply(names(max), function(x){
      c(max[[x]],sapply(grid[opt.var$max.var[i],],function(y)y, simplify = F)[[x]])
    }
    , simplify = F)
    
    opt <- sapply(names(opt), function(x){
      c(opt[[x]],sapply(grid[opt.var$opt.var[i],],function(y)y, simplify = F)[[x]])
    }
    , simplify = F)
    
  }
  opt.list <- sapply(opt[1:length(varX)], function(x) c(x), simplify = F)
  names(opt.list) <- names(varX)
  
  max.list <- sapply(max[1:length(varX)], function(x) c(x), simplify = F)
  names(max.list) <- names(varX)
  rownames(opt.var$rq) <- sapply(c(1:ncomp), function(x) paste("Comp. ", x, sep = ""))
  
  list(var.opt.x = opt.list, var.max.x = max.list,
       var.opt.y = opt[length(opt)], var.max.y = max[length(max)],
       validation.data = data.frame(grid, t(opt.var$rq)))
}



####################################
####                             ###
####        Aux functions        ###
####                             ###
####################################

# Optimize nº of comp.
optimal.comp <- function (data, option, p="Q^2", plot =TRUE, low.limit = 0.5,
                          optimal.threshold = 0.0975) {
  rec=0
  #dev.off()
  mypar <- par(no.readonly = TRUE)
  mai <- c(par()$mai[1:3],par()$pin[2]-2)
  #par(mfrow= c(1,1), xpd = TRUE, mai = mai)
  y1 <- apply(data,1,mean)
  
  if (length(which(y1>low.limit))==0) {
    #message("\n * Low limit value is not reached in any component.")
    low.limit = 0
  }
  
  if (option == "max") {
    opt <- y1[which(y1 == max(y1))[1]]
  } else {
    rq <- which(y1 == max(y1))[1]
    for ( i in rq:1) {
      if (abs(y1[i]-y1[rq])<=optimal.threshold & y1[i] > low.limit) {
        rec <- y1[i]
      }
    }
    opt <- y1[i]
    for ( i in 2:length(y1)) {
      if (abs(y1[i]-rec)>0.0975 & y1[i] > low.limit) {
        opt <- y1[i]
      } else {
        break
      }
    }
  }
  if ( plot == TRUE) {
    message("\n Qualitative result:")
    if (p == "Q^2") {
      if (opt < 0.5) {
        message("\t- Bad prediction performance.")
      } else {
        if ( opt < 0.8) {
          message("\t- Good prediction performance.")
        } else {
          message("\t- Excellent prediction performance.")
        }
      }
    } else {
      if (opt < 0.5) {
        message("\t- Bad variability explanation of Y.")
      } else {
        if ( opt < 0.8) {
          message("\t- Good variability explanation of Y.")
        } else {
          message("\t- Excellent variability explanation of Y.")
        }
      }
    }
    
    # PLot matrix
    # points and lines
    y2 <- apply(data,1,mean)
    arrow_1 <- apply(data, 1, quantile, probs = 0.025)
    arrow_2 <- apply(data, 1, quantile, probs = 0.975)

    plot(seq(1,nrow(data),1),y2, type = "b", pch = 19,
         main = "Nº of components optimization",
         xlab = "Nº components", ylab = paste("Cumulative",p), col = "darkgreen",
         ylim = c((min(data)-abs(min(data)*0.2)),
                  (max(data)+abs(max(data)*0.2))), cex = 0.8, xaxt = "n", bty = "L")
    axis(1, at=seq(1,length(y1),1), labels= 1:length(y1), tck = -0.03)
    lines(rep(which(y1 == opt)[1],2), c((min(data)-abs(min(data)*0.2)),
                                        (max(data)+abs(max(data)*0.2))), col = "red")
    lines(rep(max(which(y1 == rec)[1],2),2), c((min(data)-abs(min(data)*0.2)),
                                               (max(data)+abs(max(data)*0.2))), col = "blue", lty = 5)
    legend((nrow(data)), (max(data)+abs(max(data))*0.1),
           legend = c("SIMCA-P rule",
                      "Optimal", "Cumulative Q^2"),
           col = c("red","blue", "darkgreen"), bty = "n", lty = c(1,2,1))

    if(sum(arrow_1-arrow_2)!=0){
      arrows(x0=seq(1,nrow(data),1),y0=arrow_1,
             x1=seq(1,nrow(data),1),y1=arrow_2,
             length=0.05, angle = 90, code = 3, col = 'darkgreen')
    }
  }
  rows <- c()
  for ( i in 1:length(y1)) {
    rows <- c(rows, paste("Comp.", i))
  }
  cols <- c()
  for ( i in 1:ncol(data)) {
    cols <- c(cols, paste("Rep.", i))
  }
  rownames(data) <- rows
  colnames(data) <- cols
  # par(mypar)
  
  list(simca.comp = which(y1 == opt)[1], 
       opt.comp = max(which(y1 == rec)[1],2),
       Q2 = data)
  
  
}

# Optimize var of X for comp.
optimal.var.n <- function (data, option,
                           optimal.threshold = 0.0975) {
  y1 <- data
  
  if (option == "max") {
    opt <- y1[which(y1 == min(y1))[1]]
  } else {
    rq <- which(y1 == max(y1))[1]
    for ( i in rq:1) {
      if (abs(y1[i]-y1[rq])<=optimal.threshold) {
        opt <- y1[i]
      }
    }
  }
  list(var = which(y1 == opt)[1], value = opt)
}
optimal.var.nx <- function (data, selection, option,
                            optimal.threshold = 0.0975) {
  
  y1 <- c()
  for ( i in 1:dim(data)[1]) {
    y1 <- c(y1, data[i, selection[i]])
  }
  
  if (option == "max") {
    opt <- y1[which(y1 == min(y1))[1]]
  } else {
    rq <- which(y1 == max(y1))[1]
    for ( i in rq:1) {
      if (abs(y1[i]-y1[rq])<=optimal.threshold) {
        opt <- y1[i]
      }
    }
  }
  list(var = which(y1 == opt)[1], value = opt)
}

# Optimize global nº of X variables (without Y)
optimal.var.null <- function (data.raw, option, p, varX, varY, low.limit = 0.5,
                              optimal.threshold = 0.0975) {
  #dev.off()
  mypar <- par(no.readonly = TRUE)
  mai <- c(par()$mai[1:3],par()$pin[2]-0.2)
  #par(mfrow= c(1,1), xpd = TRUE, mai = mai)
  data <- array(0, c(nrow(data.raw[[1]]), ncol(data.raw[[1]]), length(data.raw)))
  for ( i in 1:length(data.raw)) {
    data[,,i] <- data.raw[[i]][,,1]
  }
  
  matrix <- t(apply(data,1,function(x) apply(x,1,mean)))
  # PLot matrix
  # points and lines
  arrow_1 <- t(apply(data,1,function(x) apply(x,1,quantile, probs = 0.025)))
  arrow_2 <- t(apply(data,1,function(x) apply(x,1,quantile, probs = 0.975)))
  
  if (length(which(matrix>low.limit))==0) {
    #message("\n * Low limit value is not reached in any component.")
    low.limit = 0
  }
  
  #print(matrix)
  palette(rainbow(max(dim(matrix)[1],2), s = 0.7))
  plot(seq(1,(dim(matrix)[2]+1),1),c(matrix[1,],0), pch = 19, 
       main = "Nº of X variables optimization",
       xlab = "Nº of X variables", ylab = paste("Cumulative",p), col = 1,
       ylim = c((min(data)-abs(min(data)*0.2)),
                (max(data)+abs(max(data)*0.2))), 
       cex = 0.4, type = "n", xaxt = "n", bty = "L", las = 2)
  #if(sum(arrow_1[1,]-arrow_2[1,])!=0){
  #  arrows(x0=seq(1,ncol(matrix),1),y0=arrow_1[1,],
  #         x1=seq(1,ncol(matrix),1),y1=arrow_2[1,],
  #         length=0.05, angle = 90, code = 3, col = 1)
  #}
  
  for ( i in 1:dim(matrix)[1]) {
    lines(seq(1,dim(matrix)[2],1),matrix[i,], type = "b",
          pch = 19, col = i, cex = 0.4)
    
    if(sum(arrow_1[i,]-arrow_2[i,])!=0){
      arrows(x0=seq(1,ncol(matrix),1),y0=arrow_1[i,],
             x1=seq(1,ncol(matrix),1),y1=arrow_2[i,],
             length=0.05, angle = 90, code = 3, col = i)
    }
  }
  axis(1, at=seq(1,dim(matrix)[2],1), labels= varX, tck = -0.03)
  
  # Obt. opt
  var.x <- c()
  var.max <- c()
  for ( i in 1:dim(matrix)[1]) {
    if (length(which(matrix[i,]>low.limit))==0) {
      low.limit = min(matrix[i,])
    }
    var.max <- c(var.max, which(matrix[i,] == max(matrix[i,])))
    points(which(matrix[i,] == max(matrix[i,]))[1], max(matrix[i,]), col = i, pch = 8, cex = 3)
    rq <- which(matrix[i,] == max(matrix[i,]))[1]
    for ( j in rq:1) {
      if (abs(matrix[i,][j]-matrix[i,][rq])<=optimal.threshold & matrix[i,][j] > low.limit) {
        opt <- matrix[i,][j]
      }
    }
    points(which(matrix[i,] == opt[1])[1], opt[1], col = i, pch = 18, cex = 3)
    var.x <- c(var.x,  which(matrix[i,] ==opt[1])[1])
  }
  
  # Create legend
  text <- c()
  for ( i in 1:dim(matrix)[1]) {
    text <- c(text, paste("Comp.", i))
  }
  #legend("right", legend = text, 
  #      col = c(1:dim(matrix)[1]), inset=c(-0.55,0), bty = "n", lty = rep(1,dim(matrix)[1]))
  
  #text(rep(ncol(matrix), nrow(matrix)), matrix[,ncol(matrix)], labels = text, pos = 4)
  #legend("right", legend = c(paste("Max",p), "Optimal"),
  #       pch = c(8,18), bty = "n", inset = c(-0.5,0), pt.cex = 2)
  legend((dim(matrix)[2]+0.5), 
         (max(arrow_2)+abs(max(arrow_2)*0.2)), 
         legend = c("Max", "Optimal", text), col = c("black", "black", 1:dim(data)[1]),
         pch = c(8,18, rep(NA,dim(data)[1] )), lty = c(NA, NA, rep(1,dim(data)[1])),
         bty = "n", pt.cex = 2)
  rownames(matrix) <- text
  colnames(matrix) <- varX
  par(mypar)
  names(var.x) <- text; names(var.max) <- text
  list (optvar = var.x, maxvar = var.max, rq = t(matrix))
  ##############
}

# Optimize global nº of X variables (with Y)
optimal.var <- function (data.raw, option, p, varX, varY, low.limit = 0.5,
                         optimal.threshold = 0.0975) {
  #dev.off()
  mypar <- par(no.readonly = TRUE)
  mai <- c(par()$mai[1:3],par()$pin[2]-0.2)
  par(mfrow= c(1,1), xpd = TRUE, mai = mai)
  data <- array(0, c(nrow(data.raw[[1]]), ncol(data.raw[[1]]), dim(data.raw[[1]])[3]))
  arrow_1 <- array(0, c(nrow(data.raw[[1]]), ncol(data.raw[[1]]), dim(data.raw[[1]])[3]))
  arrow_2 <- array(0, c(nrow(data.raw[[1]]), ncol(data.raw[[1]]), dim(data.raw[[1]])[3]))
  for ( i in 1:dim(data.raw[[1]])[3]) {
    aux <- array(0, c(nrow(data.raw[[1]]), ncol(data.raw[[1]]), length(data.raw)))
    for ( j in 1:length(data.raw)) {
      aux[,,j] <- data.raw[[j]][,,i]
    }
    data[,,i] <- t(apply(aux,1,function(x) apply(x,1,mean)))
    arrow_1[,,i] <- t(apply(aux,1,function(x) apply(x,1,quantile, probs = 0.025)))
    arrow_2[,,i] <- t(apply(aux,1,function(x) apply(x,1,quantile, probs = 0.975)))
  }
  # Opt y and x
  # coord
  xv <- 0
  yv <- 0
  xq <- 0
  yq <- 0
  maxv <- c()
  maxq <- c()
  optv <- c()
  optq <- c()
  
  # Plot X
  rq_per_n <- data.frame("Comp.1" = c(t(data[1,,])))
  palette(rainbow(max(dim(data)[1],2), s = 0.7))
  plot(c(1:(dim(data)[2]*dim(data)[3])), c(t(data[1,,])), pch = 19,
       main = "Nº of variables optimization",
       xlab = "Nº variables of  Y (top) and X (bottom)", ylab = paste("Cumulative",p), col = 1,
       ylim = c((min(arrow_1)-abs(min(arrow_1)*0.2)),
                (max(arrow_2)+abs(max(arrow_2)*0.2))),
       cex = 0.5, type = "b", xaxt = "n", bty = "L", las = 2)
  if (sum(c(t(arrow_1[1,,]))-c(t(arrow_2[1,,]))) != 0) {
    arrows(x0=c(1:(dim(data)[2]*dim(data)[3])),y0=c(t(arrow_1[1,,])),
           x1=c(1:(dim(data)[2]*dim(data)[3])),y1=c(t(arrow_2[1,,])),
           length=0.05, angle = 90, code = 3, col = 1)
  }
  for ( l in 2:dim(data)[1]) {
    C <- paste("Comp.",l,sep="")
    rq_per_n <- data.frame(rq_per_n, C = c(t(data[l,,])))
    lines(c(1:(dim(data)[2]*dim(data)[3])), c(t(data[l,,])), type = "b",
          pch = 19, col = l, cex = 0.5)
    if (sum(c(t(arrow_1[l,,]))-c(t(arrow_2[l,,]))) != 0) {
      arrows(x0=c(1:(dim(data)[2]*dim(data)[3])),y0=c(t(arrow_1[l,,])),
             x1=c(1:(dim(data)[2]*dim(data)[3])),y1=c(t(arrow_2[l,,])),
             length=0.05, angle = 90, code = 3, col = l)
    }
  }
  y.axe <- rep(varY, length(varX))
  x.axe <- c(sapply(varX, function(x) rep(x, length(varY))))
  axis(1, at=c(1:(dim(data)[2]*dim(data)[3])), labels= y.axe, padj = 0)
  axis(1, at=c(1:(dim(data)[2]*dim(data)[3])), labels= x.axe, padj = 1.6)
  max <- c()
  opt <- c()
  
  # COmbinations
  lista <- list("X" = varX, "Y" = varY)
  comb <- expand.grid(lista)[order(expand.grid(lista)[,1]),]
  
  for (n in 1:dim(data)[1]) {
    y1 <- c(t(data[n,,]))
    if(max(y1)<low.limit) {
      low.limit = 0
    }
    max <- c(max, (which(y1 == max(y1))[1]))
    rq <- which(y1 == max(y1))[1]
    for ( i in rq:1) {
      if (abs(y1[i]-y1[rq])<=optimal.threshold & y1[i] >= lwo.limit) {
        rec <- which(y1 == y1[i])[1]
      }
    }
    opt <- c(opt, rec)
    points(max[n], y1[max[n]], pch = 8, cex = 3, col = n)
    points(opt[n], y1[opt[n]], pch = 18, cex = 3.1, col = n)
    maxv <- c(maxv, comb[max[n],"X"])
    maxq <- c(maxq, comb[max[n],"Y"])
    optv <- c(optv, comb[opt[n],"X"])
    optq <- c(optq, comb[opt[n],"X"])
  }
  
  # Create legend
  text <- c()
  for ( i in 1:dim(data)[1]) {
    text <- c(text, paste("Comp.", i))
  }
  colnames(rq_per_n) <- text

  # Plot legend
  legend(((dim(data)[2]*dim(data)[3])+0.5), 
         (max(arrow_2)+abs(max(arrow_2)*0.2)), 
         legend = c("Max", "Optimal", text), col = c("black", "black", 1:dim(data)[1]),
         pch = c(8,18, rep(NA,dim(data)[1] )), lty = c(NA, NA, rep(1,dim(data)[1])),
         bty = "n", pt.cex = 2)
  #legend(((dim(data)[2]*dim(data)[3])),(min(data)+0.2*min(data)), legend = "X axis:\n Nº variables of\n  Y (top) and X (bottom)", pch = NA, bty = "n")#,inset = c(-0.5,0))
  print(y.axe)
  print(x.axe)
  par(mypar)
  list (optvarx = optv, optvary= optq, 
        maxvarx = maxv, maxvary = maxq, 
        rq = data.frame("X"=x.axe, "Y"=y.axe, rq_per_n))
  ############
  
}


# Optimize nº of X variables for MB-S-PLS
optimal.var.n.null.block <- function (data, option, optimal.threshold = 0.0975) {
  y1 <- data
  if (option == "max") {
    opt <- y1[which(y1 == min(y1))[1]]
  } else {
    rq <- which(y1 == max(y1))[1]
    for ( i in rq:1) {
      if (abs(y1[i]-y1[rq])<=optimal.threshold) {
        opt <- y1[i]
      }
    }
  }
  list(var = which(y1 == opt)[1], value = opt)
}

optimal.var.n.block <- function (data, option, optimal.threshold = 0.0975) {
  y1 <- data
  if (option == "max") {
    opt <- y1[which(y1 == min(y1))[1]]
  } else {
    rq <- which(y1 == max(y1))[1]
    for ( i in rq:1) {
      if (abs(y1[i]-y1[rq])<=optimal.threshold) {
        opt <- y1[i]
      }
    }
  }
  list(var = which(y1 == opt)[1], value = opt)
}

# Optimize nº of X and Y variables for MB-S-PLS
optimal.var.null.block <- function (data, option, p, grid, varX, low.limit = 0.5, 
                                    optimal.threshold = 0.0975) {
  #dev.off()
  mypar <- par(no.readonly = TRUE)
  mar <- c(par()$pin[1]-0.5,par()$mar[2:3],par()$pin[2]+8)
  par(mfrow= c(1,1), xpd = TRUE, mar = mar)
  
  # Mean and arrows
  matrix <- t(apply(data,1,function(x) apply(x,1,mean)))
  arrow_1 <- t(apply(data,1,function(x) apply(x,1,quantile, probs = 0.025)))
  arrow_2 <- t(apply(data,1,function(x) apply(x,1,quantile, probs = 0.975)))
  
  if (length(which(matrix>low.limit))==0) {
    #message("\n * Low limit value is not reached in any component.")
    low.limit = 0
  }
  
  #print(matrix)
  palette(rainbow(max(dim(matrix)[1],2), s = 0.7))
  plot(seq(1,dim(matrix)[2],1),matrix[1,], pch = 19, 
       main = "Nº of X variables optimization",
       xlab = "", ylab = paste("Cumulative",p), col = 1,
       ylim = c((min(arrow_1)-abs(min(arrow_1)*0.2)),
                (max(arrow_2)+abs(max(arrow_2)*0.2))), 
       cex = 0.4, type = "b", xaxt = "n", bty = "L")
  if(sum(arrow_1[1,]-arrow_2[1,])!=0){
    arrows(x0=seq(1,ncol(matrix),1),y0=arrow_1[1,],
           x1=seq(1,ncol(matrix),1),y1=arrow_2[1,],
           length=0.05, angle = 90, code = 3, col = 1)
  }
  
  for ( i in 2:dim(matrix)[1]) {
    lines(seq(1,dim(matrix)[2],1),matrix[i,], type = "b",
          pch = 19, col = i, cex = 0.4)
    
    if(sum(arrow_1[i,]-arrow_2[i,])!=0){
      arrows(x0=seq(1,ncol(matrix),1),y0=arrow_1[i,],
             x1=seq(1,ncol(matrix),1),y1=arrow_2[i,],
             length=0.05, angle = 90, code = 3, col = i)
    }
  }
  for ( i in 1:length(varX)) {
    axis(1, at=seq(1,dim(matrix)[2],1), labels= grid[,i], padj = 0+(1.5*(i-1)), cex.axis = 1, outer = F)
  }
  
  # Obt. opt
  var.x <- c()
  var.max <- c()
  for ( i in 1:dim(matrix)[1]) {
    if (length(which(matrix[i,]>low.limit))==0) {
      low.limit = min(matrix[i,])
    }
    var.max <- c(var.max, which(matrix[i,] == max(matrix[i,])))
    points(which(matrix[i,] == max(matrix[i,]))[1], max(matrix[i,]), col = i, pch = 8, cex = 3)
    rq <- which(matrix[i,] == max(matrix[i,]))[1]
    for ( j in rq:1) {
      if (abs(matrix[i,][j]-matrix[i,][rq])<=optimal.threshold & matrix[i,][j] > low.limit) {
        opt <- matrix[i,][j]
      }
    }
    points(which(matrix[i,] == opt[1])[1], opt[1], col = i, pch = 18, cex = 3)
    var.x <- c(var.x,  which(matrix[i,] ==opt[1])[1])
  }
  
  # Create legend
  text <- c()
  for ( i in 1:dim(matrix)[1]) {
    text <- c(text, paste("Comp.", i))
  }
  #legend("right", legend = text, 
  #      col = c(1:dim(matrix)[1]), inset=c(-0.55,0), bty = "n", lty = rep(1,dim(matrix)[1]))
  
  #text(rep(ncol(matrix), nrow(matrix)), matrix[,ncol(matrix)], labels = text, pos = 4)
  #legend("right", legend = c(paste("Max",p), "Optimal"),
  #       pch = c(8,18), bty = "n", inset = c(-0.5,0), pt.cex = 2)
  legend((dim(matrix)[2]+0.2),(max(arrow_2)+abs(max(arrow_2)*0.2)), legend = c("Max", "Optimal", text), col = c("black", "black", 1:dim(matrix)[1]),
         pch = c(8,18, rep(NA,dim(matrix)[1])), lty = c(NA, NA, rep(1,dim(matrix)[1])),
         bty = "n", pt.cex = 2)
  legend((dim(matrix)[2]+0.2),(min(matrix)+0.2*min(matrix)), legend = "X axis:\n\n Nº variables of\n  X1 (top) to\n Xn (bottom)", pch = NA, bty = "n")#,inset = c(-0.5,0))
  #rownames(matrix) <- text
  #colnames(matrix) <- varX
  par(mypar)
  #names(var.x) <- text; names(var.max) <- text
  list (opt.var = var.x, max.var = var.max, rq = matrix)
  ##############
}
optimal.var.block <- function (data, option, p, grid, varX, varY, low.limit = 0.5, optimal.threshold = 0.0975) {
  #dev.off()
  mypar <- par(no.readonly = TRUE)
  mar <- c(par()$pin[1]+2,par()$mar[2:3],par()$pin[2]+8)
  par(mfrow= c(1,1), xpd = TRUE, mar = mar)
  
  # Mean and arrows
  matrix <- t(apply(data,1,function(x) apply(x,1,mean)))
  arrow_1 <- t(apply(data,1,function(x) apply(x,1,quantile, probs = 0.025)))
  arrow_2 <- t(apply(data,1,function(x) apply(x,1,quantile, probs = 0.975)))
  
  if (length(which(matrix>low.limit))==0) {
    #message("\n * Low limit value is not reached in any component.")
    low.limit = 0
  }
  
  #print(matrix)
  palette(rainbow(max(dim(matrix)[1],2), s = 0.7))
  plot(seq(1,dim(matrix)[2],1),matrix[1,], pch = 19, 
       main = "Nº of variables optimization",
       xlab = "", ylab = paste("Cumulative",p), col = 1,
       ylim = c((min(arrow_1)-abs(min(arrow_1)*0.2)),
                (max(arrow_2)+abs(max(arrow_2)*0.2))), 
       cex = 0.4, type = "b", xaxt = "n", bty = "L")
  if(sum(arrow_1[1,]-arrow_2[1,])!=0){
    arrows(x0=seq(1,ncol(matrix),1),y0=arrow_1[1,],
           x1=seq(1,ncol(matrix),1),y1=arrow_2[1,],
           length=0.05, angle = 90, code = 3, col = 1)
  }
  
  for ( i in 2:dim(matrix)[1]) {
    lines(seq(1,dim(matrix)[2],1),matrix[i,], type = "b",
          pch = 19, col = i, cex = 0.4)
    
    if(sum(arrow_1[i,]-arrow_2[i,])!=0){
      arrows(x0=seq(1,ncol(matrix),1),y0=arrow_1[i,],
             x1=seq(1,ncol(matrix),1),y1=arrow_2[i,],
             length=0.05, angle = 90, code = 3, col = i)
    }
  }
  varib <- varX; varib[["Y"]] <- varY
  for ( i in 1:length(varib)) {
    axis(1, at=seq(1,dim(matrix)[2],1), labels= grid[,i], padj = 0+(1.5*(i-1)), cex.axis = 1, outer = F)
  }
  
  # Obt. opt
  var.x <- c()
  var.max <- c()
  for ( i in 1:dim(matrix)[1]) {
    if (length(which(matrix[i,]>low.limit))==0) {
      low.limit = min(matrix[i,])
    }
    var.max <- c(var.max, which(matrix[i,] == max(matrix[i,])))
    points(which(matrix[i,] == max(matrix[i,]))[1], max(matrix[i,]), col = i, pch = 8, cex = 3)
    rq <- which(matrix[i,] == max(matrix[i,]))[1]
    for ( j in rq:1) {
      if (abs(matrix[i,][j]-matrix[i,][rq])<=optimal.threshold & matrix[i,][j] > low.limit) {
        opt <- matrix[i,][j]
      }
    }
    points(which(matrix[i,] == opt[1])[1], opt[1], col = i, pch = 18, cex = 3)
    var.x <- c(var.x,  which(matrix[i,] ==opt[1])[1])
  }
  
  # Create legend
  text <- c()
  for ( i in 1:dim(matrix)[1]) {
    text <- c(text, paste("Comp.", i))
  }
  #legend("right", legend = text, 
  #      col = c(1:dim(matrix)[1]), inset=c(-0.55,0), bty = "n", lty = rep(1,dim(matrix)[1]))
  
  #text(rep(ncol(matrix), nrow(matrix)), matrix[,ncol(matrix)], labels = text, pos = 4)
  #legend("right", legend = c(paste("Max",p), "Optimal"),
  #       pch = c(8,18), bty = "n", inset = c(-0.5,0), pt.cex = 2)
  legend((dim(matrix)[2]+0.3),(max(arrow_2)+abs(max(arrow_2)*0.2)), legend = c("Max", "Optimal", text), col = c("black", "black", 1:dim(matrix)[1]),
         pch = c(8,18, rep(NA,dim(matrix)[1])), lty = c(NA, NA, rep(1,dim(matrix)[1])),
         bty = "n", pt.cex = 2)
  legend((dim(matrix)[2]+0.3),(min(matrix)), legend = "X axis:\n\n Nº variables of\n  X1 (top) to\n Xn (bottom).\n Y: last row.", pch = NA, bty = "n")#,inset = c(-0.5,0))
  #rownames(matrix) <- text
  #colnames(matrix) <- varX
  par(mypar)
  #names(var.x) <- text; names(var.max) <- text
  list (opt.var = var.x, max.var = var.max, rq = matrix)
  ##############
  
}

# Create fold for CV
create.fold <- function(y, k) {
  #min_reps <- 2
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0) 
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k, 
                                                               size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)
  
  out <- split(seq(along = y), foldVector)
  names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
                      sep = "")
  out
}

