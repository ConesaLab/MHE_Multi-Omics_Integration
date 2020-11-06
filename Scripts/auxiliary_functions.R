#' Score plot of a NOISeq PCA object
#'
#' A function to score plot is already programmed in NOISeq package but this is a personalized version using ggplot2
#'
#' @param myPCA Object generated with dat(type = "PCA") function.
#' @param comp Principal components to plot in axis.
#' @param title An overall title for the plot.
#'
#' @return A score plot.
#' @examples plot.scores.gg(myPCA, comp = c(1,2), title = "My score plot")

plot.scores.gg <- function(myPCA, comp = c(1,2), title = NULL){
  # dependencies
  require(ggplot2)
  
  # score values and group attributes
  da <- data.frame(x=myPCA@dat$result$scores[,comp[1]], y=myPCA@dat$result$scores[,comp[2]], Group = myPCA@dat$factors[,2])
  
  # scatter plot
  ggplot(data = da, aes(x=x, y=y, color = Group)) +
    geom_point(size=5, shape=16) +
    geom_text(aes(label=myPCA@dat$factors[,1]), hjust=0.5, vjust=-0.75, size=5) +
    scale_x_continuous(name=sprintf("PC%i (%s%%)", comp[1], round(myPCA@dat$result$var.exp[comp[1]]*100, 2))) +
    scale_y_continuous(name=sprintf("PC%i (%s%%)", comp[2], round(myPCA@dat$result$var.exp[comp[2]]*100, 2))) +
    ggtitle(title) +
    theme_bw() +
    theme(legend.background = element_rect(fill="linen", size=1, linetype="solid"), text = element_text(size=20))
}


#' Loading plot of a NOISeq PCA object
#'
#' A function to loading plot is already programmed in NOISeq package but this is a personalized version using ggplot2
#'
#' @param myPCA Object generated with dat(type = "PCA") function.
#' @param comp Principal components to plot in axis.
#' @param title An overall title for the plot.
#' @param color data.frame with 2 columns: gene name and gene characteristic (i.e. %GC, length).
#'
#' @return A loading plot.
#' @examples plot.scores.gg(myPCA, comp = c(1,2), title = "My loading plot")
#' @examples plot.scores.gg(myPCA, comp = c(1,2), title = "My loading plot", color = data.frame(gene_name, gene_length))

plot.loadings.gg <- function(myPCA, comp = c(1,2), title = NULL, color = NULL){
  # dependencies
  require(ggplot2)
  
  # loading values and group attributes
  da <- data.frame(x=myPCA@dat$result$loadings[,comp[1]], y=myPCA@dat$result$loadings[,comp[2]])
  col <- color[names(myPCA@dat$result$loadings[,1]), 2]
  
  #plot
  ggplot(data = da, aes(x=x, y=y, color = col)) +
    geom_point(size=3, shape=16) +
    scale_x_continuous(name=sprintf("PC%i (%s%%)", comp[1], round(myPCA@dat$result$var.exp[comp[1]]*100, 2))) +
    scale_y_continuous(name=sprintf("PC%i (%s%%)", comp[2], round(myPCA@dat$result$var.exp[comp[2]]*100, 2))) +
    scale_color_gradientn(colours = rainbow(3))+
    ggtitle(title) +
    theme_bw() +
    theme(text = element_text(size=20))
}


#' Personalized version of uvGsa() function of mdgsa package (Montaner)
#'
#' Performs a Uni-Variate Gene Set Analysis using a logistic regression model like the original function.
#' We added to the output the genes names that are annotated inside each GO term.
#'
#' @param index ranking index, generally a numerical named vector.
#' @param annot an annotation list.
#' @param p.adjust.method p-value adjustment method for multiple testing.
#' @param family see glm.fit.
#' @param verbose verbose.
#' @param verbosity integer indicating which iterations should be indicated if verbose = TRUE.
#' @param fulltable if TRUE, 'sd', 't' and 'convergence' indicator from the glm fit are included in the output.
#' @param ... further arguments to be pasted to glm.fit, for instance 'weights'.
#'
#' @return A data.frame with a row for each Gene Set or block. Columns are:
#' N:  number of genes annotated to the Gene Set.
#' lor:  log Odds Ratio estimated for the Gene Set.
#' pval:  p-values associated to each log Odds Ratio.
#' padj:  adjusted p-values.
#' sd:  standard deviations associated to each log Odds Ratio.
#' t:  t statistic associated to each log Odds Ratio.
#'  
#' @examples plot.scores.gg(myPCA, comp = c(1,2), title = "My score plot")

uvGsa_withGeneNames <- function (index, annot, p.adjust.method = "BY", family = quasibinomial(), 
          verbose = TRUE, verbosity = 100, fulltable = FALSE, ...) 
{
  if (is.data.frame(index)) {
    index <- as.matrix(index)
  }
  if (is.matrix(index)) {
    genes <- rownames(index)
    rownames(index) <- NULL
  }
  if (is.vector(index)) {
    genes <- names(index)
    names(index) <- NULL
  }
  if (is.null(genes)) {
    stop("no geneIds found. Check names or rownames in 'index'")
  }
  blocks <- names(annot)
  if (is.null(blocks)) {
    stop("unnamed 'annot'")
  }
  res <- matrix(NA, nrow = length(blocks), ncol = 7)
  rownames(res) <- blocks
  colnames(res) <- c("N", "lor", "sd", "t", "pval", "conv", "genenames")
  t0 <- proc.time()
  counter <- 0
  if (verbose) {
    message("Analyzed blocks:")
  }

  X <- cbind(rep(1, times = length(genes)), index)
  colnames(X) <- NULL
  for (bl in blocks) {
    B <- as.numeric(genes %in% annot[[bl]])
    gene_names = paste(genes[genes %in% annot[[bl]]], collapse = "__") ## added code
    res.glm <- glm.fit(x = X, y = B, family = family, ...)
    res.sum <- summary.glm(res.glm)
    res[bl, ] <- c(sum(B), res.sum$coefficients[2, ], res.glm$converged, gene_names)
    if (verbose) {
      counter <- counter + 1
      if (counter%%verbosity == 0) {
        message(counter, ", ", appendLF = FALSE)
      }
      if (counter%%(10 * verbosity) == 0) {
        message("\n")
      }
    }
  }
  t1 <- proc.time()
  if (verbose) {
    message("time in seconds:")
    print(t1 - t0)
  }
  if (verbose) {
    if (any(res[, "conv"] == 0) & !fulltable) {
      tex <- "The analysis did not converge for some blocks.\n                    You may re-run uvGsa using 'fulltable = TRUE' to find them."
      warning(gsub("  +", "", tex))
    }
  }
  res <- cbind(res, padj = p.adjust(res[, "pval"], method = p.adjust.method))
  res <- res[, c("lor", "pval", "padj", "sd", "t", "conv", "N", "genenames")]
  if (!fulltable) {
    res <- res[, c("lor", "pval", "padj", "N", "genenames")]
  }
  res <- as.data.frame(res)
  res[,"lor"] <- as.numeric(res[,"lor"])
  res[,"pval"] <- as.numeric(res[,"pval"])
  res[,"padj"] <- as.numeric(res[,"padj"])
  res[,"N"] <- as.numeric(res[,"N"])
  res = res[order(res$pval),]
  res
}


#' Mean plot
#'
#' Plot showing values and mean per group for one variable
#'
#' @param data matrix or data.frame with rows = obs and columns = var.
#' @param groups a character string with group category.
#' @param plot a logical indicating whether show plot or not.
#'
#' @return A data.frame with mean values where rows are groups and columns are variables.
#' @examples meanplot(data = data, groups = groups, plot = TRUE)

meanplot <- function(data, groups, plot = TRUE){ 
  #data: rows = obs, columns = var
  mean = aggregate(data, by = list(groups), FUN = mean, na.rm = T)
  rownames(mean) = mean[,1]
  mean = mean[unique(groups),-1,drop=F]
  
  if (plot) {
    for (i in 1:ncol(mean)){
      plot(x=c(1:length(unique(groups))), y=mean[,i], type = "b", cex=1, pch=19, lwd=1, axes=F, col="red" ,
           xlab="", ylab="", main=colnames(mean)[i], ylim=c(min(data[,i]), max(data[,i])))
      axis(side=1, at=1:length(unique(groups)), labels=unique(groups))
      axis(2)
      #box()
      freqs <- table(groups)[unique(groups)]
      cols = c("blue", "orange", "darkgray")
      lim_down = 1
      lim_up = freqs[1]
      for (j in 1:length(freqs)){
        text(x=rep(j, freqs[j]), y=data[lim_down:lim_up, i], cex=0.75, pch=4, lwd = 2, col=cols[j], 
             labels = rownames(data[lim_down:lim_up, i, drop=F]))
        lim_down = lim_down + freqs[j]
        lim_up = lim_up + freqs[j+1]
      }
    }
  }
  
  return(mean)
}


#' eval.clusters
#'
#' Function to evaluate clusters of a correlation matrix by silhouette mean score.
#'
#' @param matrix correlation matrix.
#' @param clusters number of clusters.
#' @param algorithm string indicating clustering algorithm. Only 2 possible values: "hclust" or "kmeans". Default: "hclust".
#' @param method when "hclust" is selected, the agglomeration method to be used. Same as hclust() function.
#' @param corrplot a logical indicating whether show corrplot. Default:FALSE.
#'
#' @return List with components of each cluster, silhouette score and number of bad classified components. Optionally, a corrplot.
#' @examples eval.clusters(matrix = cormatrix, clusters = 5, algorithm = "hclust", method = "average")

eval.clusters <- function(matrix, clusters, algorithm = "hclust", method = "average", corrplot = FALSE){
  # dependencies
  require(cluster)
  require(corrplot)
  require(pals)
  
  if(algorithm == "hclust"){
    # hclust
    # methods = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
    distance <- as.dist(sqrt(2*(1-matrix^2))) #correlation must to be transformed to distance
    ClustOrig <- hclust(distance, method = method)
    cut <- cutree(ClustOrig, k = clusters)
    main <- "hclust"
  } else if(algorithm == "kmeans"){
    # kmeans
    # methods = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")
    distance <- as.dist(sqrt(2*(1-matrix^2))) #correlation must to be transformed to distance
    cut <- kmeans(distance, centers=clusters, iter.max=10, algorithm = method)$cluster
    main <- "kmeans"
  } else if(algorithm == "pam"){
    # Partitioning Around Medoids
    # methods = c("euclidean", "manhattan")
    distance <- as.dist(sqrt(2*(1-matrix^2))) #correlation must to be transformed to distance
    cut <- pam(distance, k=clusters, diss = TRUE, metric = method)$clustering
    main <- "pam"
  } else {
    print("Only 3 possible algorithms: hclust, kmeans or pam")
  }
  
  # silhouette parameters
  sL <- cluster::silhouette(cut, distance) 
  sil <- cbind(sL[,1],sL[,2],sL[,3])
  rownames(sil) <- names(cut)
  
  silmyorder <- sil[order(sil[,1], rownames(sil), sil[,3], decreasing = c(FALSE,TRUE)),]
  matrix <- matrix[rownames(silmyorder), rownames(silmyorder)]
  
  if(corrplot == TRUE){
  # corrplot
  mycolor <- pals::glasbey()[-c(1,2)]
  corrplot::corrplot(matrix,
                     order = "original",
                     tl.cex = 1.5,
                     tl.col = mycolor[silmyorder[,1]],
                     cl.lim = c(0,1),
                     cl.pos = "b",
                     cl.cex = 1,
                     cl.length = 2,
                     # cl.ratio = 0.5,
                     # title = method,
                     mar = c(0, 0, 0, 0))
  corrplot::corrRect(table(silmyorder[,1]), col = unique(mycolor[silmyorder[,1]]), lwd = 3)
  }
  
  # components of each cluster
  pr <- list()
  for (i in 1:clusters){
    pr[[i]] = rownames(silmyorder)[which(silmyorder[,1] == i)]
  }
  names(pr) <- c(1:clusters)
  
  return(list(components = pr, score = mean(sil[,3]), bad_classified = sum(sil[,3] < 0)))
}


#' Expression profile plot
#'
#' Function to plot expression profiles with positive correlation together and negative correlation flipped.
#'
#' @param list list of clusters (output of eval.clusters() function).
#' @param expr.matrix expression matrix (genes as rows).
#' @param corM correlation matrix.
#' @param path path to the directory where plots will be stored.
#'
#' @return .
#' @examples expr.profile.plot(list = clusters.list, expr.matrix = molecules.matrix, corM = corM, path = "01_clustering/")

expr.profile.plot <- function(list, expr.matrix, corM, path){
  for(i in 1:length(list)){
    png(paste0(path, "Group", names(list[i]), ".png"))
    par(mar=c(6,7,5,2)+.1, font.axis = 2, font.lab = 2) 
    
    # center and scale data
    x = scale(expr.matrix, center=T, scale=T)
    x = x[,list[[i]],drop=F]
    # palete color equal to corrplot in eval.clusters
    col = pals::glasbey()[-c(1,2)][i]
    
    # sort cluster by expression to keep positive correlations whith the more expressed as reference
    expr_gradient = names(sort(apply(expr.matrix[,list[[i]]], 2, mean), decreasing = T))
    # correlations
    number_of_negatives <- corM[expr_gradient, expr_gradient,drop=F] > 0
    pos_cor <- rownames(number_of_negatives)[number_of_negatives[,1,drop=F] == TRUE]
    neg_cor <- rownames(number_of_negatives)[number_of_negatives[,1,drop=F] == FALSE]
    
    # plot
    matplot(x[,pos_cor,drop=F], type = "b", pch=19, lwd = 3, col = col, lty = 1, 
            main = paste("Module ", names(list[i])), 
            ylab = "Scaled expression profile", 
            ylim = c(-3, 3),
            axes=F, cex.lab = 2, cex.axis =2, cex.main = 2)
    # add condition to plot profiles with negative correlation as flipped profile (oposite math symbol)
    if (!identical(neg_cor, character(0))){
      new_lines <- -x[,neg_cor, drop=F]
      apply(new_lines, 2, lines, type = "b", pch = 1, lty=2, lwd = 2, col = col)
    }
    
    # more 
    axis(side=2, lwd = 3, cex.axis =1.5)
    axis(side=1, lwd = 3, cex.axis =1.5, at=1:5, labels=rownames(expr.matrix)[1:5], las = 2, col.axis= "blue")
    axis(side=1, lwd = 3, cex.axis =1.5, at=6:11, labels=rownames(expr.matrix)[6:11], las = 2, col.axis= "orange")
    box()
    
    dev.off()
  }
}


## PLS
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

## S-PLS
tune.comp <- function(X, Y, factor, fold, ncomp = 10, rep = 50, low.limit = 0.5, mode = "regression",
                      optimal.threshold = 0.0975,  option = "Q", title = NULL) {
  
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
                             optimal.threshold = optimal.threshold, plot = TRUE, title = title)
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
                             optimal.threshold = optimal.threshold, plot = TRUE, title = title)
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

####################################
####                             ###
####        Aux functions        ###
####                             ###
####################################

# Optimize nº of comp.
optimal.comp <- function (data, option, p="Q^2", plot =TRUE, low.limit = 0.5,
                          optimal.threshold = 0.0975, title = NULL) {
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
         main = title, #"Nº of components optimization",
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









## Barplot (Mean, SD) ----
barplot.mean.sd <- function(dataset = t(median_norm), groups = patients$class, filename = "multipage.pdf"){
  require(tidyr)
  require(dplyr)
  require(ggplot2)
  #### Mean & SD ####
  dim(dataset)
  mean_table <- aggregate(dataset~groups, FUN = mean)
  sd_table <- aggregate(dataset~groups, FUN = sd)
  
  data <- data.frame(dataset, patient = groups, stringsAsFactors = F)
  
  # Calculates Mean, SD and N
  ggplot_df <- data %>%
    tibble::rownames_to_column("Case") %>%
    tidyr::gather(-c(patient, Case), key = "Metabolite", value = "Value") %>%
    dplyr::group_by(patient, Metabolite) %>%
    dplyr::summarise(Mean = mean(Value),
                     SD = sd(Value),
                     Number = n()) %>%
    dplyr::ungroup()
  
  plots_data <- lapply(unique(ggplot_df$Metabolite), function(i){
    met_plot <- ggplot(dplyr::filter(ggplot_df, Metabolite == i)) +
      geom_bar( aes(x = factor(patient, levels=unique(groups)), y=Mean ), 
                stat="identity", fill="forestgreen", alpha=0.5) +
      geom_errorbar( aes(x=patient, ymin=Mean-SD, ymax=Mean+SD), width=0.4, colour="orange", alpha=0.9, size=1.5) +
      theme_classic() + 
      theme(axis.text.x = element_text(size = 7)) +
      theme(axis.title.x=element_blank()) +
      ggtitle(i)
    
    return(met_plot)
  })
  
  ml <- gridExtra::marrangeGrob(plots_data, nrow=1, ncol=dim(dataset)[2], grobs = plots_data)
  ggsave(filename, ml, width = 11, device="pdf")
}

## PLS ----

################################
######     MixOmics     ########
################################
# score plot PCA mixOmics


#select loadings mixomics
toploadings <- function(pls, comp){
  sort(pls$loadings$X[,comp], decreasing = T)
  names(sort(plsda.res$loadings$X[,1], decreasing = T)[1:dim(up)[1]])
}


################################
######     DaVinci      ########
################################
# score plot X Manu
score.plot.X <- function(pls, factor, sub=NULL){
plot(pls$T[,1], pls$T[,2],xlab="X score Comp.1",ylab="X score Comp.2", col = c("tomato", "steelblue2")[factor], pch = 20, cex = 3.5,  main = "Score plot (X space)", sub=sub, asp = 1)
text(pls$T[,1], pls$T[,2], labels = rownames(pls$T), pos = 4, cex=1.5)
}
# loading plot X Manu
loading.plot.X <- function(pls, factor){
  plot(pls$W[,1], pls$W[,2], col = "springgreen3", pch = 20, cex = 2, main = "Loading plot (X space)", asp = 1, ylim=c(-0.2,0.2))
  top10_pc1 <- names(sort(abs(pls$W[,1]), decreasing=T)[1:5])
  sel_pc1 <- which(rownames(pls$W) %in% top10_pc1)
  text(pls$W[sel_pc1,1][1:2], pls$W[sel_pc1,2][1:2], labels = rownames(pls$W[sel_pc1,][1:2,]), pos = 2, col = "blue", cex= 1.5)
  text(pls$W[sel_pc1,1][3:5], pls$W[sel_pc1,2][3:5], labels = rownames(pls$W[sel_pc1,][3:5,]), pos = 4, col = "blue",cex=1.5)
  points(pls$W[sel_pc1,1], pls$W[sel_pc1,2], pch=1, col = "blue", cex = )
  # 
  top10_pc2 <- names(sort(abs(pls$W[,2]), decreasing=T)[1:5])
  sel_pc2 <- which(rownames(pls$W) %in% top10_pc2)
  text(pls$W[sel_pc2,1][1], pls$W[sel_pc2,2][1], labels = rownames(pls$W[sel_pc2,][1,,drop=F]), pos = 3, col = "orange",cex=1.5)
  text(pls$W[sel_pc2,1][2], pls$W[sel_pc2,2][2], labels = rownames(pls$W[sel_pc2,][2,,drop=F]), pos = 2, col = "orange",cex=1.5)
  text(pls$W[sel_pc2,1][3:4], pls$W[sel_pc2,2][3:4], labels = rownames(pls$W[sel_pc2,][3:4,,drop=F]), pos = 4, col = "orange",cex=1.5)
  text(pls$W[sel_pc2,1][5], pls$W[sel_pc2,2][5], labels = rownames(pls$W[sel_pc2,][5,,drop=F]), pos = 1, col = "orange",cex=1.5)
  points(pls$W[sel_pc2,1], pls$W[sel_pc2,2], pch=1, col = "orange")
  
  # top10B = names(sort(abs(pls$B[,1]), decreasing=T)[1:10])
  # selB <- which(rownames(pls$B) %in% top10B)
  # text(pls$W[selB,1], pls$W[selB,2], labels = rownames(pls$W[selB,]), pos = 3)
  # points(pls$W[selB,1], pls$W[selB,2], pch=1, col = "black", cex = )
}

# score plot Y Manu
score.plot.Y <- function(pls, factor, sub=NULL){
  plot(pls$U[,1], pls$U[,2], col = c("tomato", "steelblue2")[factor], pch = 20, cex = 3.5, main = "Score plot (Y space)", sub=sub, asp = 1)
  text(pls$U[,1], pls$U[,2], labels = rownames(pls$U), pos = 4, cex = 1.5)
}
# loading plot Y Manu
loading.plot.Y <- function(pls, factor){
  plot(pls$C[,1], pls$C[,2], col = "springgreen3", pch = 20, cex = 3.5, main = "Loading plot (Y space)", asp = 1)
  top10 <- names(sort(abs(pls$C[,1]), decreasing=T)[1:10])
  sel <- which(rownames(pls$C) %in% top10)
  text(pls$C[sel,1], pls$C[sel,2], labels = rownames(pls$C[sel,,drop=F]), pos = 4, cex = 1.5)
}



# # Modify margins to plot legend outside plot -----
# initpar <- par(c("mai", "pin", "xpd"))
# mai <- c(initpar$mai[1:3], initpar$pin[2]+3)
# # on.exit(par(mai = initpar$mai, xpd = initpar$xpd)) #if you use it in a function
# par(mai = mai)
# 
# # After plot...
# par(xpd = TRUE)
# legend((ncol(design.matrix)-0.17), 2,
#        legend = c("Hight",
#                   "Medium",
#                   "Low"),
#        col = c("indianred","gold", "lightgreen"),
#        bty = "n", lty = c(2,2,2),
#        title = "Batch effect magnitude")

## GOS ANOVA-----

# GOS: goodnes of separation score plot (ANOVA) -------
score.anova <- function(pls, group, space = "X"){
  # This function calculates the line between 2 centroids of groups and then project all values to this line. 
  # Then calculates non-parametric ANOVA between proyected values.
  # pls: object from pls2() function --> Manu's code
  # group: factor of conditions
  # space: "X" or "Y"
  
  group <- as.factor(group)
  
  if (space == "X"){
    mydata = pls$T
  } else if (space == "Y"){
    mydata = pls$U
  } else {
    print("Only X or Y values to 'space' are allowed")
  }

  ncomp <- ncol(mydata)
  if (ncomp == 1){
    # ANOVA
    aov = kruskal.test(mydata[,1] ~ group)
  } else{
  # PROYECTIONS
    # Centroid calculation
    centroid1 = apply(mydata[1:4,,drop=F], 2, median)
    centroid2 = apply(mydata[5:10,,drop=F], 2, median)
    
    x1 = centroid1[1]
    y1 = centroid1[2]
    x2 = centroid2[1]
    y2 = centroid2[2]
    
    slope <- (y2-y1)/(x2-x1) #m = (y2-y1)/(x2-x1)
    intercept <- y1-slope*x1 #b = y - ax
    
    # Proyection
    vector_dir <- c((x2-x1), (y2-y1)) #v = (x2-x1), (y2-y1) --> https://www.geogebra.org/m/UAafu2Qc
    proyection <- mydata %*% vector_dir
    
    
    # # Plots
    # score.plot.X(pls, factor)
    # abline(h=0, col = "grey", lty=2)
    # abline(v=0, col = "grey", lty=2)
    # points(x1,y1, col="darkolivegreen2", pch=18, cex=2)
    # points(x2,y2, col="darkolivegreen2", pch=18, cex=2)
    # abline(a=intercept, b=slope, col="darkolivegreen2", lwd = 2)
    # 
    # # Proyection plot
    # score.plot.X(list("T" = cbind(proyection, proyection*0)), factor)
    # abline(h=0, col = "grey", lty=2)
    # abline(v=0, col = "grey", lty=2)
    
    # ANOVA
    aov = kruskal.test(proyection[,1] ~ group)
  }
  return(aov)
}


## GOS F-stat-----
score.stat <- function(pls, groups, Fcrit = 4.7374){
  #Groups info
  n1=4
  n2=6
  # Centroid calculation
  U = pls$U
  c1 = apply(U[1:4,,drop=F], 2, mean)
  c2 = apply(U[5:10,,drop=F], 2, mean)
  # Euclidean distance: |AB| = sqrt( (x2-x1)^2 + (y2-y1)^2 )
  De = sqrt( (c2[1]-c1[1])^2 + (c2[2]-c1[2])^2 )
  # Mahalanobis distance: D^2 = (x - ??)' ??^-1 (x - ??)
  centroids = (rbind(c1,c2))
  d = centroids[2,,drop=F] - centroids[1,,drop=F] # d = [x2 - x1, y2 - y1]
  S = cov(U) #score matrix
  mah = mahalanobis(d, center=F, cov=S)
  # T2 Hotelling
  t2 = ((n1*n2)/(n1+n2)) * mah
  # F-statistic
  p=2
  f = ((n1+n2-p-1) / (p*(n1+n2-2))) * t2
  # P.value
  pval = pf(f, df1=2, df2=7,lower.tail = F)
  pval = as.numeric(pval)
  ### Compute table ###
  table = data.frame("Eu_dist"= De, "Mah_dist" = mah, "T2"=t2, "F.val"=f, "Crit_Fval"=Fcrit, "Pval"=pval)
  return(table)
}

## mdgsa (Montaner)----
# Para que aparezcan los genes implicados en cada pathway, hay que cambiar el c??digo de la funcion uvGsa
# fix(uvGsa)
# B <- as.numeric(genes %in% annot[[bl]])
# b = paste(genes[genes %in% annot[[bl]]], collapse = "__")
# res.glm <- glm.fit(x = X, y = B, family = family, ...)
# res.sum <- summary.glm(res.glm)
# res[bl, ] <- c((b), res.sum$coefficients[2, ], res.glm$converged)
#

## ORA (Over Representation Analysis)----
EnrichALLterms = function (test, notTest, annotation, p.adjust.method = "fdr") {
  
  annot2test = unique(annotation[,2])
  
  resultat2 = t(sapply(annot2test, Enrich1term, test = test, notTest = notTest, annotation = annotation))
  
  return (data.frame(resultat2, 
                     "adjPval" = p.adjust(as.numeric(resultat2[,"pval"]), method = p.adjust.method), 
                     stringsAsFactors = F))
  
}

Enrich1term = function (term, test, notTest, annotation) {
  annotTest = length(intersect(test, annotation[annotation[,2] == term,1]))
  annotTest_names = intersect(test, annotation[annotation[,2] == term,1])
  annotNOTtest = length(intersect(notTest, annotation[annotation[,2] == term,1]))
  # annotNOTtest_names = intersect(notTest, annotation[annotation[,2] == term,1])
  if (identical(annotTest_names, character(0))){
    annotTest_names = 0
  }
  
  if ((annotTest) > 0) {
    annotNOTtest = length(intersect(notTest, annotation[annotation[,2] == term,1]))
    mytest = matrix(c(annotTest, length(test)-annotTest, annotNOTtest, length(notTest)-annotNOTtest), ncol = 2)
    resultat = c(term, annotTest, length(test), annotNOTtest, length(notTest), 
                 fisher.test(mytest, alternative = "greater")$p.value, paste(sort(annotTest_names), collapse="__"))
    names(resultat) = c("term", "annotTest", "test", "annotNotTest", "notTest", "pval", "annotTest_names")
  } else {
    resultat = c(term, 0, 0, 0, 0, 100, 0)
    names(resultat) = c("term", "annotTest", "test", "annotNotTest", "notTest", "pval", "annotTest_names")
  }
  
  return(resultat)
  
}
