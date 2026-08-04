# =============================================================================
# Data application: rheumatoid arthritis (Methotrexate) — analysis only
#
# Paper: "Incremental Causal Effect for Time to Treatment Initialization"
#        (Ying, Zhao & Xu, ICLR 2025), Section 4.2.
#
# This script starts from an ALREADY-PREPARED, one-row-per-patient data frame
# `df` and estimates the incremental causal effect of the time to initiate
# Methotrexate (MTX) on the number of tender joints at one year, with
# multiplier-bootstrap pointwise 95% confidence intervals.
#
# -----------------------------------------------------------------------------
# REQUIRED INPUT: a data frame `df` with one row per patient and the columns:
#   A       : observed time to MTX initiation, administratively censored at tau
#             (i.e. A = min(time-to-MTX, tau); numeric)
#   Delta   : event indicator, 1 if MTX was initiated before tau, else 0
#   Y        : outcome — number of tender joints at month tau (numeric)
#   <covariates in L.ext> : baseline covariates (see L.ext below), numeric
#
# tau is the follow-up horizon (tau = 12 months in the paper).
# `Y` must be non-missing (patients with missing Y should already be removed).
#
# The construction of `df` from raw per-visit records (and the raw file
# "ra.xlsx") is intentionally omitted here; that data is patient-level and not
# distributed. Provide your own prepared `df` before running this script.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
library(survival)   # coxph, Surv, cox.zph — Cox model and PH test
library(timereg)    # cox.aalen — cumulative martingale residual check
library(ggplot2)    # plotting psi_hat(theta)

# ---- Analysis settings ------------------------------------------------------
tau <- 12  # follow-up horizon (months); Y is measured at tau, A censored at tau

# Baseline covariates L used to model the MTX-initiation hazard.
# These column names must exist in `df`.
L.ext <- c("age_0", "sex", "smoke_0", "edu_0", "duration_0", "year_0",
           "rapos", "haqc_0", "gsc_0", "esrc_0", "dmrd_0", "onprd2_0", "jc_0")

# (Optional) quick sanity checks on the prepared data:
# stopifnot(all(c("A", "Delta", "Y", L.ext) %in% names(df)))
# df <- df[!is.na(df$Y), ]
# print(nrow(df))   # 1010 patients in the paper

# -----------------------------------------------------------------------------
# 1) Cox proportional-hazards (PH) goodness-of-fit checks
# -----------------------------------------------------------------------------
surv_obj <- Surv(time = df$A, event = df$Delta)

# Schoenfeld residuals test of the PH assumption.
# H0: residuals are independent of time (PH holds) for each covariate.
formula <- as.formula(paste("surv_obj ~ ", paste(L.ext, collapse = " + ")))
multi_fit <- coxph(formula, data = df)
print(cox.zph(multi_fit))   # in the paper, year_0 and dmrd_0 violate PH

# Cumulative martingale residuals (a second PH check) via cox.aalen.
formula2 <- as.formula(paste("surv_obj ~", paste0("prop(", L.ext, ")", collapse = " + ")))
multi_fit2 <- cox.aalen(formula2, data = df)
summary(multi_fit2)         # p-values for the test of proportionality
plot(multi_fit2, score = TRUE)   # cumulative martingale residual plots

# -----------------------------------------------------------------------------
# 2) Drop covariates that violate the PH assumption
# -----------------------------------------------------------------------------
L.ext <- setdiff(L.ext, c("year_0", "dmrd_0"))
df <- df[, setdiff(names(df), c("year_0", "dmrd_0"))]

# -----------------------------------------------------------------------------
# 3) IPW estimation of psi(theta) with multiplier-bootstrap pointwise 95% CIs
# -----------------------------------------------------------------------------
cand_set <- seq(0.2, 5, 0.1)  # grid of hazard-ratio interventions theta
B <- 200                       # number of multiplier-bootstrap draws
N <- nrow(df)                  # sample size

# IPW estimator (Theorem 1): fit a Cox model for the MTX-initiation hazard,
# then reweight Y by  theta^Delta * exp(-(theta - 1) * cumhaz).
# The multiplier bootstrap re-fits the Cox model with exp(1) weights to get SEs.
real_apply <- function() {
  est_set <- c()

  surv_obj <- Surv(time = df$A, event = df$Delta)
  formula <- as.formula(paste("surv_obj ~ ", paste(L.ext, collapse = " + ")))

  # Point estimate over the theta grid
  fit <- coxph(formula, data = df)
  cumhaz <- predict(fit, type = "expected")   # per-subject cumulative hazard
  for (hz_ratio in cand_set) {
    est_set <- c(est_set, mean(hz_ratio^(df$Delta) * exp(-(hz_ratio - 1) * cumhaz) * df$Y))
  }

  # Multiplier (weighted) bootstrap for standard errors
  boot_est_set <- c()
  for (b in 1:B) {
    boot_weights <- rexp(N)
    boot_weights <- boot_weights / mean(boot_weights)
    boot_est <- c()
    fit <- coxph(formula, data = df, weights = boot_weights)
    cumhaz <- predict(fit, type = "expected")
    for (hz_ratio in cand_set) {
      boot_est <- c(boot_est, mean(boot_weights * hz_ratio^(df$Delta) * exp(-(hz_ratio - 1) * cumhaz) * df$Y))
    }
    boot_est_set <- rbind(boot_est_set, boot_est)
  }

  sd_set <- apply(boot_est_set, 2, sd)
  CI_U_set <- est_set + qnorm(0.975) * sd_set   # pointwise 95% CI upper
  CI_L_set <- est_set + qnorm(0.025) * sd_set   # pointwise 95% CI lower

  return(data.frame(theta = cand_set,
                    est = est_set,
                    sd = sd_set,
                    CI_U = CI_U_set,
                    CI_L = CI_L_set))
}

set.seed(12345)
result <- real_apply()

# -----------------------------------------------------------------------------
# 4) Plot psi_hat(theta) with the pointwise 95% CI band
# -----------------------------------------------------------------------------
ggplot(result, aes(x = theta, y = est)) +
  geom_line(aes(color = "Estimate"), linewidth = 0.5) +
  geom_ribbon(aes(ymin = CI_L, ymax = CI_U, fill = "Pointwise 95% CI"), alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  labs(y = expression(paste("Average number of tender joints ", hat(psi)(theta))),
       x = expression(paste("Hazards ratio ", theta))) +
  theme_classic() +
  scale_x_log10(limits = range(result$theta), breaks = c(0.2, 0.5, 1, 2, 5)) +
  scale_color_manual(name = "", values = c("Estimate" = "black")) +
  scale_fill_manual(name = "", values = c("Pointwise 95% CI" = "lightblue")) +
  theme(legend.position = "top",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_line(color = "black"))

# -----------------------------------------------------------------------------
# 5) Numerical summaries reported in the paper
#    (percent change in tender joints relative to the observed regime, theta = 1)
# -----------------------------------------------------------------------------
round(100 * (1 - result$est[cand_set == 2] / result$est[cand_set == 1]), 2)   # theta = 2
round(100 * (1 - result$est[cand_set == 5] / result$est[cand_set == 1]), 2)   # theta = 5
round(100 * (result$est[cand_set == .5] / result$est[cand_set == 1] - 1), 2)  # theta = 0.5
round(100 * (result$est[cand_set == .2] / result$est[cand_set == 1] - 1), 2)  # theta = 0.2

round(result$est[c(19, 49, 4, 1)], 2)    # point estimates at selected theta
round(result$CI_U[c(19, 49, 4, 1)], 2)   # their pointwise CI upper bounds
round(result$CI_L[c(19, 49, 4, 1)], 2)   # their pointwise CI lower bounds
