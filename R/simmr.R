#' Similarity-based multimodal regression
#' 
#' `simmr` fits a multimodal regression model using only distance matrices 
#' computed on the original data. Permutation testing is used to obtain 
#' *p*-values.
#' 
#' @param D List of distance matrices of class \link[base]{matrix} or 
#'   \link[stats]{dist}
#' @param X Design matrix, either output of \link[stats]{model.matrix} or valid
#'   input to the function.
#' @param variables Names of explanatory variables to regress on, must match
#'   with column names of `X`
#' @param tests Chosen test statistics from "Dempster", "PC". For "PC", defaults
#'   to using n-p-1 PCs, where `X` is `n x p`. If `pc.pv` and `pc.n` are both
#'   specified, `pc.pv` is prioritized.
#' @param D.scale Method for normalization of distance matrices
#' @param pc.pv Optional, calculate additional PC-based test statistics 
#'   for specified proportions of variation explained
#' @param pc.n Optional, calculate additional PC-based test statistics 
#'   for specific numbers of PCs. If 'all', tries every possible number of PCs.
#' @param n.perm Number of permutations to perform
#' @param eigen.tol Discard eigenvalues below this tolerance for both PCA and
#'   computing test statistics
#' @param debug Include internal objects in output (from unpermuted analysis)
#'
#' @return A list with the following components:
#'   \item{`stat`}{Test statistics computed on original data.}
#'   \item{`perms`}{Test statistics after permutations.}
#'   \item{`perm.p`}{Permutation p-values for each test.}
#'   
#' @import MDMR, stats
#' @export
#'
#' @examples
#' D <- list(dist(rnorm(10)), dist(rnorm(10)), dist(rnorm(10)))
#' X <- list("var" = runif(10))
#' simmr(D, X, "var")
simmr <- function(D, X, variables, tests = c("Dempster", "PC"), 
                  D.scale = c("max", "trace", "dvar", "none"), 
                  pc.pv = NULL,
                  pc.n = NULL,
                  n.perm = 999,
                  eigen.tol = 1e-10, 
                  debug = FALSE) {
  tests <- match.arg(tests, several.ok = TRUE)
  D.scale <- match.arg(D.scale)
  
  # If X is not a matrix, convert to model matrix and define reduced model
  if (!("matrix" %in% class(X))) {
    redcols <- unlist(sapply(variables, 
                      function(v) paste(v, levels(X[[v]]), sep = "")))
    X <- model.matrix(~., X)
    red <- !(colnames(X) %in% redcols)
  } else {
    red <- !(colnames(X) %in% variables)
  }
  
  design <- X
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #### Obtain MDS scores for each distance matrix ####
  Zlist <- lapply(D, function(x) {
    A <- -as.matrix(x)^2/2
    As <- A - rep(colMeans(A), rep.int(n, n))
    G <- t(As) - rep(rowMeans(As), rep.int(n, n))
    
    Gsvd <- svd(G)
    # Normalize Gower's matrices by chosen criteria
    switch(
      D.scale,
      "max" = {
        Gsvd$u %*% diag(sqrt(Gsvd$d/max(Gsvd$d)))
      },
      "trace" = {
        Gsvd$u %*% diag(sqrt(Gsvd$d/sum(Gsvd$d)))
      },
      "dvar" = {
        Gsvd$u %*% diag(sqrt(Gsvd$d/sqrt(sum(G^2/n^2))))
      },
      "none" = {
        Gsvd$u %*% diag(sqrt(Gsvd$d))
      })
  })
  Z <- do.call(cbind, Zlist)

  #### Perform PCA on MDS axes if specified ####
  Zpca <- NULL
  if ("PC" %in% tests) {
    if (identical(pc.n, "all")) {
      pc.n <- 1:(n-p-1)
    }
    
    Zpca <- prcomp(Z)
    ev <- Zpca$sdev^2
    pv <- ev/sum(ev)
    tpv <- cumsum(pv)
    
    # prioritize using pc.pv over pc.n
    if (!is.null(pc.pv)) {
      pc.n <- sapply(pc.pv, function(x) which(tpv > x)[1])
    }
    
    # keep maximum number of PCs (above eigen.tol) and specified numbers
    Zpc <- Zpca$x[,1:min(n-p-1, max(which(ev > eigen.tol)))]
    scores <- Zpca$x[,1:max(which(ev > eigen.tol))]
    
    Zpc_l <- lapply(pc.n, function(pcn) {
      scores[,1:min(dim(scores)[2], pcn), drop = FALSE]
    })
    
    # name PC output as pv if pc.pv is specified
    if (!is.null(pc.pv)) {
      names(Zpc_l) <- as.character(pc.pv*100)
    } else if (!is.null(pc.n)) {
      names(Zpc_l) <- paste0("_n", pc.n)
    }
  }
  
  #### Compute specified test statistics ####
  out <- NULL
  stats <- NULL
  for (i in 1:(n.perm+1)) {
    if (i > 1) {X <- design[sample(1:n),]}
    stat <- NULL
    
    #### Get SSCP error and regression ####
    h <- tcrossprod(tcrossprod(X, solve(crossprod(X))), X)
    E <- (t(Z) %*% (diag(n) - h) %*% Z)
    # reduced model
    Xred <- X[,red]
    hred <- tcrossprod(tcrossprod(Xred, solve(crossprod(Xred))), Xred)
    H <- (t(Z) %*% (diag(n) - hred) %*% Z) - E
    
    if (debug & (i == 1)) {
      out <- c(out, list(E = E, H = H))
    }
    
    sapply(
      tests, switch,
      "Dempster" = {
        stat <- c(stat,
                  "Dempster" = sum(diag(H))/sum(diag(E)))
      },
      "PC" = {
        # get SSCP matrices using PCs
        Epc <- (t(Zpc) %*% (diag(n) - h) %*% Zpc)
        Hpc <- (t(Zpc) %*% (diag(n) - hred) %*% Zpc) - Epc
        
        PCs_res <- NULL
        if (!is.null(pc.pv) | !is.null(pc.n)) {
          PCs_res <- lapply(Zpc_l, function(z) {
            Epc = (t(z) %*% (diag(n) - h) %*% z)
            Hpc = (t(z) %*% (diag(n) - hred) %*% z) - Epc
            HEpc = Hpc %*% solve(Epc)
            etapc = Re(eigen(HEpc)$values)
            etapc = etapc[etapc > eigen.tol]
            list(Epc = Epc, Hpc = Hpc, HEpc = HEpc, etapc = etapc)
          })
        }
        
        HEpc <- Hpc %*% solve(Epc)
        etapc <- Re(eigen(HEpc)$values)
        etapc <- etapc[etapc > eigen.tol]
        
        # get test statistics
        stat <- c(stat, 
                  "PC.Wilks" = prod(1/(1+etapc)),
                  "PC.Pillai" = sum(etapc/(1+etapc)),
                  "PC.Hotelling" = sum(etapc),
                  "PC.Roy" = etapc[1]/(1+etapc[1]))
        if (!is.null(pc.pv) | !is.null(pc.n)) {
          stat_l <- sapply(PCs_res, function(x) {
            etapc <- x$etapc
            c(
              "PC.Wilks" = prod(1/(1+etapc)),
              "PC.Pillai" = sum(etapc/(1+etapc)),
              "PC.Hotelling" = sum(etapc),
              "PC.Roy" = etapc[1]/(1+etapc[1])
            )
          })
          
          stats_l <- unlist(data.frame(stat_l))
          names(stats_l) <- c(outer(rownames(stat_l), colnames(stat_l), 
                                    FUN = paste0))
          
          stat <- c(stat, stats_l)
        }
        
        if (debug & (i == 1)) {
          out <- c(out,
                   list(Zpca = Zpca, PCs.res = PCs_res,
                        Epc = Epc, Hpc = Hpc))
        }
      }
    )
    
    stats <- rbind(stats, stat)
  }
  
  rownames(stats) <- NULL
  perm_p <- apply(stats, 2, function(x) mean(x >= x[1]))
  
  # Calculate Wilk's lambda p-values properly
  wilks <- grepl("Wilks", names(perm_p))
  perm_p[wilks] <- 1 - perm_p[wilks] + 1/(n.perm+1)
  
  out <- c(out,
           list(stat = stats[1,],
                perms = stats[-1,],
                perm.p = perm_p))
  return(out)
}

