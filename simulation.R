# =============================================================================
# Simulation study for "Incremental Causal Effect for Time to Treatment
# Initialization" (Ying, Zhao & Xu, ICLR 2025), Section 4.1.
#
# Evaluates the finite-sample performance of the IPW estimator of the
# incremental causal effect psi(theta) = E[Y_{T(theta)}], where the treatment-
# initiation hazard lambda(t|L) is scaled by a constant factor theta. Reports
# bias, empirical SE (SEE), bootstrap-based SE (SD), and 95% coverage (CP)
# across a grid of theta values and sample sizes n = 200, 1000, 5000.
#
# NOTE: rep_num (Monte Carlo replications) and B (multiplier-bootstrap draws)
# are set to small values below for a quick test run. The paper uses
# R = 1000 replications and B = 200 bootstrap draws (Table 1); set them
# accordingly to reproduce the published numbers (this is computationally heavy).
# =============================================================================

library(survival)   # coxph(), Surv() for the nuisance hazard model

# ---- Data-generating parameters (see Section 4.1 of the paper) --------------
para_set <- list(mu_L = 0.5,
                 sigma_L = 0.5,
                 beta_L = 0.25,   # covariate effect on the treatment-time hazard
                 mu_LY = 1.5,     # covariate effect on the outcome
                 mu_AY = 1,       # treatment-time effect on the outcome
                 sigma_Y = 0.5    # outcome noise SD
)

# ---- Generate one simulated dataset -----------------------------------------
# L: baseline covariate ~ Unif(0,1)
# A: time to treatment initiation, exponential with hazard scaled by hz_ratio (theta)
# Delta: event indicator 1(A < tau), with administrative censoring at tau = 2
# Y: continuous outcome depending on L and the (censored) treatment time A
data_gen <- function(N, para_set, hz_ratio) {
  L <- runif(N, 0, 1)
  A <- rexp(N, hz_ratio * exp(para_set$beta_L * L))
  Delta <- A < 2
  A <- pmin(A, 2)   # administrative censoring at tau = 2
  Y <- rnorm(N, exp(1 - para_set$mu_LY * L - para_set$mu_AY * (2 - A)), para_set$sigma_Y)

  df <- data.frame(L, A, Y, Delta)
  return(df)
}

# ---- Ground-truth psi(theta) via a very large Monte Carlo sample ------------
# Approximate the true incremental effect for each theta by directly generating
# data under the intervened hazard (hz_ratio = theta) and averaging Y.
set.seed(12345)
N <- 10^7
cand_set <- c(1/3, 1/2.5, 1/2, 1/1.5, 1.5, 2, 2.5, 3)  # grid of theta values
truth_set <- c()

for (hz_ratio in cand_set) {
  df <- data_gen(N, para_set, hz_ratio)
  truth_set <- c(truth_set, mean(df$Y))
}

rep_num <- 1   # R: number of Monte Carlo replications  (paper: 1000)
B <- 2         # B: number of multiplier-bootstrap draws (paper: 200)

# ---- One simulation run: estimate psi(theta), its bootstrap SE, and CP -------
simu <- function() {
  est_set <- c()
  sd_set <- c()
  CP_set <- c()
  rep <- 1
  while(rep <= rep_num) {
    # Generate observed data under the *factual* hazard (theta = 1)
    df <- data_gen(N, para_set, 1)
    est <- c()
    # Nuisance: fit Cox model for the treatment-initiation hazard lambda(t|L)
    fit <- coxph(Surv(time = A, event = Delta) ~ L, data = df)
    est_sd <- c()
    cumhaz <- predict(fit, type = "expected")  # estimated cumulative hazard per subject
    # IPW estimator (Theorem 1): reweight Y by the Radon-Nikodym weights
    #   theta^Delta * exp(-(theta - 1) * cumhaz)
    for (hz_ratio in cand_set) {
      est <- c(est, mean(hz_ratio^(df$Delta) * exp(-(hz_ratio - 1) * cumhaz) * df$Y))
    }
    est_set <- rbind(est_set, est)

    # Multiplier (weighted) bootstrap for inference:
    # draw exp(1) weights, refit Cox with those weights, recompute the estimator
    boot_est_set <- c()
    for (b in 1:B) {
      boot_weights <- rexp(N)
      boot_weights <- boot_weights/mean(boot_weights)
      boot_est <- c()
      fit <- coxph(Surv(time = A, event = Delta) ~ L, data = df, weights = boot_weights)
      cumhaz <- predict(fit, type = "expected")
      for (hz_ratio in cand_set) {
        boot_est <- c(boot_est, mean(boot_weights * hz_ratio^(df$Delta) * exp(-(hz_ratio - 1) * cumhaz) * df$Y))
      }
      boot_est_set <- rbind(boot_est_set, boot_est)
    }
    # Bootstrap SE and Wald-type 95% coverage indicator (per theta)
    sd_set <- rbind(sd_set, apply(boot_est_set, 2, sd))
    CP_set <- rbind(CP_set, truth_set <= est + qnorm(0.975) * apply(boot_est_set, 2, sd) & truth_set >= est + qnorm(0.025) * apply(boot_est_set, 2, sd))
    rep = rep + 1
  }
  return(list(est_set = est_set,
              sd_set = sd_set,
              CP_set = CP_set))
}

# ---- Summarize a simulation result into the reported rows -------------------
# Rows: Bias, empirical SE (SEE), average bootstrap SE (SD), 95% coverage (CP)
result_gen <- function(result) {
  return(with(result,
              rbind(colMeans(est_set) - truth_set,
                    apply(est_set, 2, sd),
                    colMeans(sd_set),
                    colMeans(CP_set))))
}

# ---- Format a result matrix as LaTeX table rows (as in Table 1) -------------
print_table <- function(my_matrix) {
  add_ampersand <- Vectorize(function(x) paste("&", x, sep = ""))
  my_matrix <- round(my_matrix * 100, digits = 3)
  my_matrix <- cbind(c("Bias", "SEE", "SD", "95\\% CP"), my_matrix)
  my_matrix <- apply(my_matrix, c(1, 2), add_ampersand)
  my_matrix <- cbind(my_matrix, c("\\\\", "\\\\", "\\\\", "\\\\"))
  #my_matrix <- rbind(c("theta =", cand_set), my_matrix)
  return(my_matrix)
}

# ---- Run for a sample size and print the LaTeX table rows ----------------
N <- 200
result <- simu()
my_matrix <- result_gen(result)
print_table(my_matrix)
