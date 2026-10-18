# =============================================================================
# Bootstrap confidence intervals for cross-validated AUC
# Two approaches:
# bootstrap_pooled_auc()   - post hoc, resamples the stored out-of-fold
#  predictions. Seconds to run. Conditional on the fitted models.


bootstrap_auc_multinom <- function(fit, n_boot = 2000, seed = 42,
                                   stratified = TRUE) {

  preds <- fit$predictions
  lev   <- levels(preds$obs_class)
  if (is.null(lev)) lev <- sort(unique(as.character(preds$obs_class)))

  auc <- function(p, y) {
    n1 <- sum(y); n0 <- length(y) - n1
    if (n1 == 0 || n0 == 0) return(NA_real_)
    (sum(rank(p)[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }

  auc_all <- function(P, obs) {
    vapply(seq_along(lev), function(k) auc(P[, k], obs == lev[k]), numeric(1))
  }

  set.seed(seed)

  out <- lapply(fit$comparison$Model, function(m) {
    d <- preds[preds$Model == m, ]

    d <- aggregate(d[, lev, drop = FALSE],
                   by = list(rowIndex = d$rowIndex, obs = d$obs_class),
                   FUN = mean)

    P   <- as.matrix(d[, lev, drop = FALSE])
    obs <- as.character(d$obs)

    strata <- split(seq_len(nrow(d)), obs)

    b <- replicate(n_boot, {
      i <- if (stratified) {
        unlist(lapply(strata,
                      function(v) v[sample.int(length(v), length(v), replace = TRUE)]),
               use.names = FALSE)
      } else {
        sample.int(nrow(d), nrow(d), replace = TRUE)
      }
      a <- auc_all(P[i, , drop = FALSE], obs[i])
      c(macro = if (anyNA(a)) NA_real_ else mean(a), a)
    })

    ci <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2) return(c(NA_real_, NA_real_, NA_real_))
      c(sd(x), quantile(x, c(0.025, 0.975), names = FALSE))
    }

    obs_auc <- auc_all(P, obs)
    st      <- t(apply(b, 1, ci))

    res <- data.frame(
      Model            = m,
      AUC_Macro_pooled = if (anyNA(obs_auc)) NA_real_ else mean(obs_auc),
      boot_se          = st[1, 1],
      ci_low           = st[1, 2],
      ci_high          = st[1, 3],
      n_valid          = sum(is.finite(b[1, ])),
      n_patients       = nrow(d),
      stringsAsFactors = FALSE
    )

    for (k in seq_along(lev)) {
      res[[paste0("AUC_", lev[k])]]          <- obs_auc[k]
      res[[paste0("AUC_", lev[k], "_low")]]  <- st[k + 1, 2]
      res[[paste0("AUC_", lev[k], "_high")]] <- st[k + 1, 3]
    }

    res
  })

  merge(fit$comparison, do.call(rbind, out), by = "Model",
        all.x = TRUE, sort = FALSE)
}