#' Simulation code for similarity-based multimodal regression (SiMMR)
#' 
#' `sim_simmr` generates data and compares SiMMR test statistics to alternative
#' methods. Further details are available in the manuscript
#'
#' @param dir Directory to save results as a .Rdata file.
#' @param N Number of observations.
#' @param M Number of modalities.
#' @param Q Number of features per modality.
#' @param pc PCs to shift, can be selected at "all", "none", "##pct", or the
#'   number of PCs as a numeric value.
#' @param x.shift Amount to multiply the PC shifts by
#' @param rho Correlation values used in `cov`, a vector of length two if `cov`
#'   is chosen as "block"
#' @param cov Correlation structure, can be "exchangeable", "ar1" for AR(1), or
#'   "block" for separating within and between modality correlations.
#' @param one.modality If TRUE, only the first modality is related to the
#'   simulated outcome.
#' @param continuous If TRUE, the simulated covariate is drawn from a standard
#'   normal distribution, otherwise drawn from a Bernoulli distribution with 
#'   probability 0.5.
#' @param debug Whether to print debug messages.
#'
#' @return
#' 
#' @import MDMR, mvtnorm, stats
#' 
#' @export
#'
#' @examples
#' # Generate 10 observations, each with 2 modalities of 5 features
#' # Exchangable correlation with rho=0.5
#' # Binary covariate with shifts along the first PC
#' simulation("", 10, 3, 2, 5, 1, 0.3, 0.5)
simulation <- function(dir, N, M, Q, pc, x.shift, rho, 
                       cov = "exchangeable", one.modality = FALSE,
                       continuous = TRUE, debug = FALSE) {
  # break if file exists
  path <- paste0(dir, "/sim_mmdmr")
  fname <- paste(path, "N", N, "M", M, "Q", Q, "PC", pc, "res.Rdata", sep = "_")
  if (file.exists(fname)) {
    break
  }
  
  # PCs to shift, can be selected at "all", "none", "##pct", or a number of PCs
  if (is.na(as.numeric(pc))) {
    if (pc == "all") {
      x.pc <- 1:(M*Q)
    } else if (pc == "none") {
      x.pc <-  1
      x.shift <- 0
    } else {
      pc.pct <- gsub("pct", "", pc)
      x.pc <- 1:(floor(as.numeric(pc.pct)/100*M*Q))
    }
  } else {
    x.pc <- 1:as.numeric(pc)
  }
  
  # choose a correlation structure
  # if block, rho is a vector of within correlation and between correlation
  if (cov == "exchangeable") {
    Sigma <- diag(1-rho, M*Q) + rho
  } else if (cov == "ar1") {
    expt <- abs(matrix(1:(M*Q) - 1, nrow = M*Q, ncol = M*Q, byrow = TRUE) - 
                  (1:(M*Q) - 1))
    Sigma <- rho^expt
  } else if (cov == "block") {
    Sigma <- matrix(rho[2], M*Q, M*Q) # between correlations
    for (m in 1:M) {
      Sigma[((m-1)*Q+1):(m*Q), ((m-1)*Q+1):(m*Q)] <- diag(1-rho[1], Q) + rho[1]
    }
  }
  
  phi <- eigen(Sigma)$vectors[,x.pc, drop = FALSE]
  pc.shift <- x.shift*rowSums(phi)
  
  if (one.modality) {pc.shift[-(1:Q)] <- 0}
  
  # store input parameters and initial sim params
  sim_params <- as.list(environment())
  
  sim_out <- NULL
  mdmr_res <- NULL
  mdmr_p_all <- NULL
  mdmr_p <- vector("numeric", nsim)
  evs <- NULL
  sim_p <- NULL
  
  lm_p <- vector("numeric", nsim)
  for (s in 1:nsim) {
    print(paste0("Simulation: ", s))
    
    # generate covariate
    if (continuous) {
      x <- rnorm(N)
    } else {
      x <- rbinom(N, 1, x.prop)
    }
    
    X <- cbind(`Intercept` = 1, x)
    
    D <- vector("list", M)
    
    # generate correlated sets of features
    Yall <- rmvnorm(N, sigma = Sigma)
    Yall <- Yall + x %*% t(pc.shift)
    Y <- NULL
    for (m in 1:M) {
      # split into separate feature sets and get distance matrices
      Y[[m]] <- Yall[,((m-1)*Q+1):(m*Q)]
      D[[m]] <- dist(Y[[m]])
    }
    
    # separate MDMR with bonferroni
    invisible(capture.output(
      mdmr_res <- lapply(D, function(d) mdmr(data.frame(x), d, nperm = 999))
      ))
    mdmr_ps <- sapply(mdmr_res, function(r) r$pv[2,1])
    # save MDMR p-values
    names(mdmr_ps) <- paste0("MDMR.", 1:M)
    mdmr_p_all <- rbind(mdmr_p_all, mdmr_ps)
    mdmr_p[s] <- min(p.adjust(mdmr_ps, method = "bonferroni"))
    
    # MMR concatenating Y
    if (M*Q < N-2) {
      lm_res <- lm(Yall ~ X)
      lm_p[s] <- anova(lm_res)$`Pr(>F)`[2]
    } else {
      lm_res <- NULL
      lm_p[s] <- NA
    }
    
    # multiple MDMR
    res <- simmr(D, X, "x", n.perm = 999, 
                 pc.n = "all")
    sim_p <- rbind(sim_p, res$perm.p)
  }
  sim_p <- cbind(sim_p, "MC_MDMR" = mdmr_p, "MMR" = lm_p, mdmr_p_all)
  
  if (debug) {
    c(sim_out, list(params = sim_params, sim.p = sim_p))
  } else {
    save(sim_p, sim_params,
         file = paste(path, "N", N, "M", M, "Q", Q, "PC", pc, 
                      "res.Rdata", sep = "_"))
  }
}