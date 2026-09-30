bmzdp_log_h_ijr <- function(t, zr_ij, gamma_ij, beta_r, mu_j, rec) {
  if (rec$piecewise) {
    tr <- 0
    for (q in seq_along(rec$lambda_n)) {
      tr <- tr + rec$lambda_n[q] * ifelse(t > rec$time_q[q] & t <= rec$time_q[q + 1], 1, 0)
    }
    val <- log(tr) + log(gamma_ij) + beta_r %*% zr_ij + mu_j
  } else {
    val <- log(rec$kesi) + (rec$kesi - 1) * log(t) + log(gamma_ij) + beta_r %*% zr_ij + mu_j
  }
  val
}

bmzdp_log_s_ijr <- function(t, zr_ij, gamma_ij, beta_r, mu_j, rec) {
  lambda_ijr <- gamma_ij * exp(c(beta_r %*% zr_ij) + mu_j)
  if (rec$piecewise) {
    sr <- 0
    for (q in seq_along(rec$lambda_n)) {
      sr <- sr + rec$lambda_n[q] * max(0, min(rec$time_q[q + 1] - rec$time_q[q], t - rec$time_q[q]))
    }
    val <- (-sr * lambda_ijr)
  } else {
    val <- (-lambda_ijr * (t^rec$kesi))
  }
  val
}

bmzdp_s_ijr <- function(t, zr_ij, gamma_ij, beta_r, mu_j, rec) {
  lambda_ijr <- gamma_ij * exp(c(beta_r %*% zr_ij) + mu_j)
  if (rec$piecewise) {
    sr <- 0
    for (q in seq_along(rec$lambda_n)) {
      sr <- sr + rec$lambda_n[q] * max(0, min(rec$time_q[q + 1] - rec$time_q[q], t - rec$time_q[q]))
    }
    val <- exp(-lambda_ijr * sr)
  } else {
    val <- exp(-lambda_ijr * t^rec$kesi)
  }
  val
}

bmzdp_log_h_ijd <- function(t, zd_ij, sigma, gamma_ij, zeta, mu_intercept, beta_d, tao, mu_j) {
  log(sigma) + ((sigma) - 1) * log(t) - (zeta * sigma) * log(gamma_ij) - ((mu_intercept + beta_d %*% zd_ij + tao * mu_j) * sigma)
}

bmzdp_log_s_ijd <- function(t, zd_ij, sigma, gamma_ij, zeta, mu_intercept, beta_d, tao, mu_j) {
  lambda_ijd <- (gamma_ij)^(-zeta * sigma) * exp(-(mu_intercept + c(beta_d %*% zd_ij) + tao * mu_j) * sigma)
  (-lambda_ijd * (t^(sigma)))
}

bmzdp_s_ijd <- function(t, zd_ij, sigma, gamma_ij, zeta, mu_intercept, beta_d, tao, mu_j) {
  lambda_ijd <- (gamma_ij)^(-zeta * sigma) * exp(-(mu_intercept + c(beta_d %*% zd_ij) + tao * mu_j) * sigma)
  exp(-lambda_ijd * t^(sigma))
}

bmzdp_log_l_ij <- function(u_ij, delta_ij, d_ij_min, t_vector, y_ij, zd_ij, sigma, gamma_ij, zeta, mu_intercept, beta_d, tao, mu_j, zr_ij, beta_r, rec) {
  if (u_ij == 1 & delta_ij == 1) {
    sum_result <- 0
    for (i in seq_along(t_vector)) {
      sum_result <- sum_result + bmzdp_log_h_ijr(t_vector[i], zr_ij, gamma_ij, beta_r, mu_j, rec)
    }
    val <- sum_result + bmzdp_log_s_ijr(d_ij_min, zr_ij, gamma_ij, beta_r, mu_j, rec) + bmzdp_log_s_ijd(d_ij_min, zd_ij, sigma, gamma_ij, zeta, mu_intercept, beta_d, tao, mu_j)
  } else if (u_ij == 1 & delta_ij == 0) {
    sum_result <- 0
    for (i in seq_along(t_vector)) {
      sum_result <- sum_result + bmzdp_log_h_ijr(t_vector[i], zr_ij, gamma_ij, beta_r, mu_j, rec)
    }
    val <- sum_result + bmzdp_log_s_ijr(d_ij_min, zr_ij, gamma_ij, beta_r, mu_j, rec) + bmzdp_log_s_ijd(d_ij_min, zd_ij, sigma, gamma_ij, zeta, mu_intercept, beta_d, tao, mu_j) + bmzdp_log_h_ijd(d_ij_min, zd_ij, sigma, gamma_ij, zeta, mu_intercept, beta_d, tao, mu_j)
  } else if (u_ij == 0 & delta_ij == 0) {
    val <- (bmzdp_log_s_ijr(d_ij_min, zr_ij, gamma_ij, beta_r, mu_j, rec) + bmzdp_log_s_ijd(d_ij_min, zd_ij, sigma, gamma_ij, zeta, mu_intercept, beta_d, tao, mu_j) + bmzdp_log_h_ijd(d_ij_min, zd_ij, sigma, gamma_ij, zeta, mu_intercept, beta_d, tao, mu_j)) * (1 - y_ij)
  } else if (u_ij == 0 & delta_ij == 1) {
    val <- log(y_ij + (1 - y_ij) * bmzdp_s_ijd(d_ij_min, zd_ij, sigma, gamma_ij, zeta, mu_intercept, beta_d, tao, mu_j) * bmzdp_s_ijr(d_ij_min, zr_ij, gamma_ij, beta_r, mu_j, rec))
  }
  val
}

bmzdp_loglik_l <- function(l, dat, st) {
  bmzdp_log_l_ij(dat$u[l], dat$delta[l], dat$time[l], dat$tv[[l]], st$y[l], dat$zd[[l]], st$sigma[l], st$gamma[l], st$zeta, st$mu_intercept, st$betad, st$tao, st$mu_hos[dat$cluster[l]], dat$zr[[l]], st$betar, st$rec)
}

bmzdp_loglik_sum <- function(idx, dat, st, nan_stop) {
  sloglik <- 0
  for (l in idx) {
    sloglik <- sloglik + bmzdp_loglik_l(l, dat, st)
    if (nan_stop && is.nan(sloglik)) {
      break
    }
  }
  sloglik
}

bmzdp_lh_eta <- function(eta_val, k, eta_temp, z_matrix, y_prev, sigma2eta) {
  eta_temp[k] <- eta_val
  logit_prob <- z_matrix %*% eta_temp
  prob <- 1/(1 + exp(-logit_prob))
  sum(dbinom(x = y_prev, size = rep(1, times = length(y_prev)), prob = prob, log = TRUE)) - (eta_val^2)/(2 * sigma2eta)
}

bmzdp_lh_betar <- function(betar_val, k, dat, st, sigma2r) {
  st$betar[k] <- betar_val
  bmzdp_loglik_sum(seq_len(dat$n), dat, st, FALSE) - betar_val^2/(2 * sigma2r)
}

bmzdp_lh_lambda_n <- function(lambda_n_val, s, dat, st, sigma2l) {
  st$rec$lambda_n[s] <- lambda_n_val
  bmzdp_loglik_sum(seq_len(dat$n), dat, st, FALSE) - lambda_n_val^2/(2 * sigma2l)
}

bmzdp_lh_kesi <- function(kesi_val, dat, st, akesi, bkesi) {
  st$rec$kesi <- kesi_val
  bmzdp_loglik_sum(seq_len(dat$n), dat, st, FALSE) + (akesi - 1) * log(kesi_val) - (kesi_val/bkesi)
}

bmzdp_lh_gamma_sub <- function(gamma_sub_val, k, dat, st, epislon_k) {
  st$gamma[k] <- gamma_sub_val
  val <- bmzdp_loglik_l(k, dat, st)
  val - ((log(gamma_sub_val))^2)/(2 * epislon_k)
}

bmzdp_lh_zeta <- function(zeta_val, dat, st, sigma2zeta) {
  st$zeta <- zeta_val
  bmzdp_loglik_sum(seq_len(dat$n), dat, st, FALSE) - (zeta_val)^2/(2 * sigma2zeta)
}

bmzdp_lh_mu_intercept <- function(mu_intercept_val, dat, st, sigma2mu_intercept) {
  st$mu_intercept <- mu_intercept_val
  bmzdp_loglik_sum(seq_len(dat$n), dat, st, FALSE) - (mu_intercept_val)^2/(2 * sigma2mu_intercept)
}

bmzdp_lh_betad <- function(betad_val, k, dat, st, sigma2d) {
  st$betad[k] <- betad_val
  bmzdp_loglik_sum(seq_len(dat$n), dat, st, FALSE) - betad_val^2/(2 * sigma2d)
}

bmzdp_lh_tao <- function(tao_val, dat, st, sigma2tao) {
  st$tao <- tao_val
  bmzdp_loglik_sum(seq_len(dat$n), dat, st, FALSE) - (tao_val)^2/(2 * sigma2tao)
}

bmzdp_lh_b_dp_cluster <- function(b_dp_cluster_val, in_k, dat, st, v2) {
  st$mu_hos[in_k] <- b_dp_cluster_val
  bmzdp_loglik_sum(seq_len(dat$n - 1), dat, st, TRUE) - b_dp_cluster_val^2/(2 * v2)
}

bmzdp_lh_b_dp_cluster_n <- function(b_dp_cluster_val_n, in_k, dat, st, v2) {
  st$sigma[in_k] <- b_dp_cluster_val_n
  bmzdp_loglik_sum(seq_len(dat$n - 1), dat, st, FALSE) - b_dp_cluster_val_n^2/(2 * v2)
}

bmzdp_draw_label <- function(prob) {
  which(rmultinom(1, size = 1, prob = prob)[, 1] == 1)
}

bmzdp_onehot_index <- function(m) {
  as.integer(colSums(m * seq_len(nrow(m))))
}

bmzdp_stick <- function(lab, K, alpha) {
  w_inter <- rep(0, times = K)
  w_inter[K] <- 1
  for (k in 1:(K - 1)) {
    s_larger_cluster <- 0
    for (s in (k + 1):K) {
      s_larger_cluster <- s_larger_cluster + sum(lab == s)
    }
    w_inter[k] <- rbeta(1, (1 + sum(lab == k)), alpha + s_larger_cluster)
  }
  w <- rep(0, times = K)
  w[1] <- w_inter[1]
  for (k in 2:K) {
    w[k] <- w_inter[k]
    for (s in 1:(k - 1)) {
      w[k] <- w[k] * (1 - w_inter[s])
    }
  }
  w
}

bmzdp_rep <- function(x, m, name) {
  if (!(length(x) %in% c(1L, m))) {
    stop(sprintf("'%s' must have length 1 or %d", name, m), call. = FALSE)
  }
  rep_len(as.numeric(x), m)
}
