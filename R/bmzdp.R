bmzdp <- function(event_times, time, delta, cluster, zr, zd, z,
                  baseline, time_q,
                  mcmc_samples, num_dp_cluster, num_dp_cluster_n, alpha, v2,
                  sigma2r, sigma2d, sigma2l, akesi, bkesi,
                  sigma2zeta, sigma2mu_intercept, sigma2tao, sigma2eta, a0, b0,
                  metrop_var_betar, metrop_var_betad, metrop_var_lambda_n, metrop_var_kesi,
                  metrop_var_gamma_sub, metrop_var_eta, metrop_var_zeta, metrop_var_mu_intercept,
                  metrop_var_tao, metrop_var_b_dp_cluster, metrop_sd_b_dp_cluster_n,
                  betar_init, betad_init, lambda_n_init, kesi_init, gamma_init, epislon_init,
                  eta_init, zeta_init, mu_intercept_init, tao_init, mu_hos_init, sigma_init,
                  b_dp_cluster_init, b_dp_cluster_n_init) {
  baseline <- match.arg(baseline, c("piecewise", "weibull"))
  piecewise <- baseline == "piecewise"

  zr <- as.matrix(zr)
  zd <- as.matrix(zd)
  z <- as.matrix(z)
  storage.mode(zr) <- "double"
  storage.mode(zd) <- "double"
  storage.mode(z) <- "double"
  time <- as.numeric(time)
  delta <- as.numeric(delta)
  n <- length(time)
  pr <- ncol(zr)
  pd <- ncol(zd)
  pc <- ncol(z)

  if (!is.list(event_times) || length(event_times) != n) stop("'event_times' must be a list with one element per subject", call. = FALSE)
  if (length(delta) != n || length(cluster) != n) stop("'time', 'delta' and 'cluster' must have the same length", call. = FALSE)
  if (nrow(zr) != n || nrow(zd) != n || nrow(z) != n) stop("'zr', 'zd' and 'z' must have one row per subject", call. = FALSE)
  if (!all(delta %in% c(0, 1))) stop("'delta' must be coded 0/1", call. = FALSE)
  if (n < 2) stop("at least two subjects are required", call. = FALSE)
  if (any(is.na(cluster)) || any(cluster != round(cluster)) || any(cluster < 1)) stop("'cluster' must contain integer codes 1, ..., J", call. = FALSE)
  if (mcmc_samples < 2) stop("'mcmc_samples' must be at least 2", call. = FALSE)
  if (num_dp_cluster < 2 || num_dp_cluster_n < 2) stop("'num_dp_cluster' and 'num_dp_cluster_n' must be at least 2", call. = FALSE)

  cluster <- as.integer(cluster)
  num_hos <- max(cluster)
  num_rec <- lengths(event_times)
  u <- ifelse(num_rec > 0, 1, 0)
  tv <- lapply(seq_len(n), function(l) if (num_rec[l] > 0) as.numeric(event_times[[l]]) else 0)

  if (piecewise) {
    time_q <- as.numeric(time_q)
    if (length(time_q) < 2 || is.unsorted(time_q, strictly = TRUE)) stop("'time_q' must be strictly increasing with at least two knots", call. = FALSE)
    num_piece <- length(time_q) - 1
    all_t <- unlist(tv[num_rec > 0])
    if (length(all_t) > 0 && (min(all_t) <= time_q[1] || max(all_t) > time_q[num_piece + 1])) stop("all recurrent event times must lie in (time_q[1], time_q[length(time_q)]]", call. = FALSE)
  }

  dat <- list(
    n = n,
    u = u,
    delta = delta,
    time = time,
    tv = tv,
    zd = lapply(seq_len(n), function(l) zd[l, ]),
    zr = lapply(seq_len(n), function(l) zr[l, ]),
    cluster = cluster
  )
  members <- lapply(seq_len(num_hos), function(k) which(cluster == k))

  metrop_var_betar <- bmzdp_rep(metrop_var_betar, pr, "metrop_var_betar")
  metrop_var_betad <- bmzdp_rep(metrop_var_betad, pd, "metrop_var_betad")
  metrop_var_gamma_sub <- bmzdp_rep(metrop_var_gamma_sub, n, "metrop_var_gamma_sub")
  metrop_var_eta <- bmzdp_rep(metrop_var_eta, pc, "metrop_var_eta")
  metrop_var_b_dp_cluster <- bmzdp_rep(metrop_var_b_dp_cluster, num_dp_cluster, "metrop_var_b_dp_cluster")
  metrop_sd_b_dp_cluster_n <- bmzdp_rep(metrop_sd_b_dp_cluster_n, num_dp_cluster_n, "metrop_sd_b_dp_cluster_n")

  betar <- matrix(0, nrow = mcmc_samples, ncol = pr)
  betar[1, ] <- bmzdp_rep(betar_init, pr, "betar_init")
  betar_accept <- rep(0, times = pr)

  if (piecewise) {
    metrop_var_lambda_n <- bmzdp_rep(metrop_var_lambda_n, num_piece, "metrop_var_lambda_n")
    lambda_n <- matrix(0, nrow = mcmc_samples, ncol = num_piece)
    lambda_n[1, ] <- bmzdp_rep(lambda_n_init, num_piece, "lambda_n_init")
    lambda_n_accept <- rep(0, times = num_piece)
  } else {
    kesi <- rep(0, times = mcmc_samples)
    kesi[1] <- kesi_init
    kesi_accept <- 0
  }

  gamma_sub <- matrix(0, nrow = mcmc_samples, ncol = n)
  gamma_sub[1, ] <- bmzdp_rep(gamma_init, n, "gamma_init")
  gamma_sub_accept <- rep(0, times = n)

  epislon <- matrix(0, nrow = mcmc_samples, ncol = num_hos)
  epislon[1, ] <- bmzdp_rep(epislon_init, num_hos, "epislon_init")

  eta <- matrix(0, nrow = mcmc_samples, ncol = pc)
  eta[1, ] <- bmzdp_rep(eta_init, pc, "eta_init")
  eta_temp <- eta[1, ]
  eta_accept <- rep(0, times = pc)

  y_index <- which(u == 0)
  p <- matrix(0, nrow = mcmc_samples, ncol = length(y_index))
  z_matrix <- z[y_index, , drop = FALSE]
  logit_prob_initial <- z_matrix %*% eta[1, ]
  p[1, ] <- 1/(1 + exp(-logit_prob_initial))

  y <- matrix(0, nrow = mcmc_samples, ncol = n)
  y[1, y_index] <- rbinom(length(y_index), 1, p[1, ])

  zeta <- rep(0, times = mcmc_samples)
  zeta[1] <- zeta_init
  zeta_accept <- 0

  mu_intercept <- rep(0, times = mcmc_samples)
  mu_intercept[1] <- mu_intercept_init
  mu_intercept_accept <- 0

  betad <- matrix(0, nrow = mcmc_samples, ncol = pd)
  betad[1, ] <- bmzdp_rep(betad_init, pd, "betad_init")
  betad_accept <- rep(0, times = pd)

  tao <- rep(0, times = mcmc_samples)
  tao[1] <- tao_init
  tao_accept <- 0

  mu_hos <- matrix(0, nrow = mcmc_samples, ncol = num_hos)
  mu_hos[1, ] <- bmzdp_rep(mu_hos_init, num_hos, "mu_hos_init")
  sigma <- matrix(0, nrow = mcmc_samples, ncol = n)
  sigma[1, ] <- bmzdp_rep(sigma_init, n, "sigma_init")

  s_label <- matrix(0L, nrow = mcmc_samples, ncol = num_hos)
  s_label_n <- matrix(0L, nrow = mcmc_samples, ncol = n)
  s_label[1, ] <- bmzdp_onehot_index(rmultinom(num_hos, size = 1, prob = rep((1/num_dp_cluster), times = num_dp_cluster)))
  s_label_n[1, ] <- bmzdp_onehot_index(rmultinom(n, size = 1, prob = rep((1/num_dp_cluster_n), times = num_dp_cluster_n)))
  pai <- rep(0, times = num_dp_cluster)
  pai_n <- rep(0, times = num_dp_cluster_n)

  b_dp_cluster <- matrix(0, nrow = mcmc_samples, ncol = num_dp_cluster)
  b_dp_cluster[1, ] <- bmzdp_rep(b_dp_cluster_init, num_dp_cluster, "b_dp_cluster_init")
  b_dp_cluster_n <- matrix(0, nrow = mcmc_samples, ncol = num_dp_cluster_n)
  b_dp_cluster_n[1, ] <- bmzdp_rep(b_dp_cluster_n_init, num_dp_cluster_n, "b_dp_cluster_n_init")
  b_dp_cluster_accept <- rep(0, times = num_dp_cluster)
  b_dp_cluster_accept_n <- rep(0, times = num_dp_cluster_n)

  w_dp_cluster <- matrix(0, nrow = mcmc_samples, ncol = num_dp_cluster)
  w_dp_cluster_n <- matrix(0, nrow = mcmc_samples, ncol = num_dp_cluster_n)
  w_dp_cluster[1, ] <- rep((1/num_dp_cluster), times = num_dp_cluster)
  w_dp_cluster_n[1, ] <- rep((1/num_dp_cluster_n), times = num_dp_cluster_n)

  st <- list(
    betar = betar[1, ],
    gamma = gamma_sub[1, ],
    y = y[1, ],
    sigma = sigma[1, ],
    zeta = zeta[1],
    mu_intercept = mu_intercept[1],
    betad = betad[1, ],
    tao = tao[1],
    mu_hos = mu_hos[1, ],
    rec = if (piecewise) list(piecewise = TRUE, lambda_n = lambda_n[1, ], time_q = time_q) else list(piecewise = FALSE, kesi = kesi[1])
  )

  for (i in 2:mcmc_samples) {
    betar[i, ] <- betar[(i - 1), ]
    for (k in seq_len(pr)) {
      betar_proposed <- rnorm(n = 1, mean = betar[(i - 1), k], sd = sqrt(metrop_var_betar[k]))
      ratio_1 <- exp(bmzdp_lh_betar(betar_proposed, k, dat, st, sigma2r) - bmzdp_lh_betar(betar[(i - 1), k], k, dat, st, sigma2r))
      if (ratio_1 >= runif(n = 1, min = 0, max = 1)) {
        betar[i, k] <- betar_proposed
        betar_accept[k] <- betar_accept[k] + 1
      }
      st$betar[k] <- betar[i, k]
    }

    if (piecewise) {
      lambda_n[i, ] <- lambda_n[(i - 1), ]
      for (s in seq_len(num_piece)) {
        lambda_n_proposed <- abs(rnorm(n = 1, mean = lambda_n[(i - 1), s], sd = sqrt(metrop_var_lambda_n[s])))
        if (lambda_n_proposed > 0) {
          ratio <- exp(bmzdp_lh_lambda_n(lambda_n_proposed, s, dat, st, sigma2l) - bmzdp_lh_lambda_n(lambda_n[(i - 1), s], s, dat, st, sigma2l))
          if (ratio >= runif(n = 1, min = 0, max = 1)) {
            lambda_n[i, s] <- lambda_n_proposed
            lambda_n_accept[s] <- lambda_n_accept[s] + 1
          }
        }
        st$rec$lambda_n[s] <- lambda_n[i, s]
      }
    } else {
      kesi[i] <- kesi[(i - 1)]
      kesi_proposed <- rnorm(n = 1, mean = kesi[(i - 1)], sd = sqrt(metrop_var_kesi))
      if (kesi_proposed > 0) {
        ratio <- exp(bmzdp_lh_kesi(kesi_proposed, dat, st, akesi, bkesi) - bmzdp_lh_kesi(kesi[(i - 1)], dat, st, akesi, bkesi))
        if (ratio >= runif(n = 1, min = 0, max = 1)) {
          kesi[i] <- kesi_proposed
          kesi_accept <- kesi_accept + 1
        }
      }
      st$rec$kesi <- kesi[i]
    }

    gamma_sub[i, ] <- gamma_sub[(i - 1), ]
    for (k in seq_len(n)) {
      gamma_sub_proposed <- rnorm(n = 1, mean = gamma_sub[(i - 1), k], sd = sqrt(metrop_var_gamma_sub[k]))
      if (gamma_sub_proposed > 0) {
        epislon_k <- epislon[(i - 1), cluster[k]]
        ratio <- exp(bmzdp_lh_gamma_sub(gamma_sub_proposed, k, dat, st, epislon_k) - bmzdp_lh_gamma_sub(gamma_sub[(i - 1), k], k, dat, st, epislon_k))
        if (!is.nan(ratio)) {
          if (ratio >= runif(n = 1, min = 0, max = 1)) {
            gamma_sub[i, k] <- gamma_sub_proposed
            gamma_sub_accept[k] <- gamma_sub_accept[k] + 1
          }
        }
      }
      st$gamma[k] <- gamma_sub[i, k]
      if ((gamma_sub_accept[k]/i) > 0.5) {
        metrop_var_gamma_sub[k] <- metrop_var_gamma_sub[k] + 5
      }
      if ((gamma_sub_accept[k]/i) < 0.2) {
        metrop_var_gamma_sub[k] <- metrop_var_gamma_sub[k] - 10
        if (metrop_var_gamma_sub[k] <= 0) {
          metrop_var_gamma_sub[k] <- 0.00000001
        }
      }
    }

    for (k in seq_len(num_hos)) {
      hos_ind_chosen <- members[[k]]
      epislon[i, k] <- 1/rgamma(1, (a0 + (length(hos_ind_chosen)/2)), (b0 + ((sum((log(gamma_sub[i, hos_ind_chosen]))^2))/2)))
    }

    eta[i, ] <- eta[(i - 1), ]
    y_prev <- y[(i - 1), y_index]
    for (k in seq_len(pc)) {
      eta_proposed <- rnorm(n = 1, mean = eta[(i - 1), k], sd = sqrt(metrop_var_eta[k]))
      ratio <- exp(bmzdp_lh_eta(eta_proposed, k, eta_temp, z_matrix, y_prev, sigma2eta) - bmzdp_lh_eta(eta[(i - 1), k], k, eta_temp, z_matrix, y_prev, sigma2eta))
      if (ratio >= runif(n = 1, min = 0, max = 1)) {
        eta[i, k] <- eta_proposed
        eta_accept[k] <- eta_accept[k] + 1
      }
      eta_temp[k] <- eta[i, k]
    }
    p[i, ] <- 1/(1 + exp(-z_matrix %*% eta[i, ]))

    for (j in seq_along(y_index)) {
      k <- y_index[j]
      st$y[k] <- 1
      val_1 <- bmzdp_loglik_l(k, dat, st)
      val_1 <- val_1 + log(p[i, j])
      st$y[k] <- 0
      val_0 <- bmzdp_loglik_l(k, dat, st)
      val_0 <- val_0 + log(1 - p[i, j])
      max_num <- max(val_1, val_0)
      val_1 <- val_1 - max_num
      val_0 <- val_0 - max_num
      val_1 <- exp(val_1)/sum(exp(val_1) + exp(val_0))
      y[i, k] <- rbinom(1, 1, val_1)
      st$y[k] <- y[i, k]
    }

    for (l in seq_len(n)) {
      for (s in seq_len(num_dp_cluster_n)) {
        st$sigma[l] <- b_dp_cluster_n[(i - 1), s]
        pai_n[s] <- log(w_dp_cluster_n[(i - 1), s]) + bmzdp_loglik_l(l, dat, st)
      }
      pai_n <- pai_n - max(pai_n)
      pai_n <- exp(pai_n)/sum(exp(pai_n))
      s_label_n[i, l] <- bmzdp_draw_label(pai_n)
      st$sigma[l] <- b_dp_cluster_n[(i - 1), s_label_n[i, l]]
    }

    b_dp_cluster_n[i, ] <- b_dp_cluster_n[(i - 1), ]
    for (k in seq_len(num_dp_cluster_n)) {
      in_k <- s_label_n[i, ] == k
      b_dp_propose_n <- abs(rnorm(n = 1, mean = b_dp_cluster_n[(i - 1), k], sd = metrop_sd_b_dp_cluster_n[k]))
      ratio_2 <- exp(bmzdp_lh_b_dp_cluster_n(b_dp_propose_n, in_k, dat, st, v2) - bmzdp_lh_b_dp_cluster_n(b_dp_cluster_n[(i - 1), k], in_k, dat, st, v2))
      if (ratio_2 >= runif(n = 1, min = 0, max = 1)) {
        b_dp_cluster_n[i, k] <- b_dp_propose_n
        b_dp_cluster_accept_n[k] <- b_dp_cluster_accept_n[k] + 1
      }
      sigma[i, in_k] <- b_dp_cluster_n[i, k]
      st$sigma[in_k] <- b_dp_cluster_n[i, k]
    }

    w_dp_cluster_n[i, ] <- bmzdp_stick(s_label_n[i, ], num_dp_cluster_n, alpha)

    zeta[i] <- zeta[(i - 1)]
    zeta_proposed <- rnorm(n = 1, mean = zeta[(i - 1)], sd = sqrt(metrop_var_zeta))
    ratio <- exp(bmzdp_lh_zeta(zeta_proposed, dat, st, sigma2zeta) - bmzdp_lh_zeta(zeta[(i - 1)], dat, st, sigma2zeta))
    if (ratio >= runif(n = 1, min = 0, max = 1)) {
      zeta[i] <- zeta_proposed
      zeta_accept <- zeta_accept + 1
    }
    st$zeta <- zeta[i]

    mu_intercept[i] <- mu_intercept[(i - 1)]
    mu_intercept_proposed <- rnorm(n = 1, mean = mu_intercept[(i - 1)], sd = sqrt(metrop_var_mu_intercept))
    ratio <- exp(bmzdp_lh_mu_intercept(mu_intercept_proposed, dat, st, sigma2mu_intercept) - bmzdp_lh_mu_intercept(mu_intercept[(i - 1)], dat, st, sigma2mu_intercept))
    if (ratio >= runif(n = 1, min = 0, max = 1)) {
      mu_intercept[i] <- mu_intercept_proposed
      mu_intercept_accept <- mu_intercept_accept + 1
    }
    st$mu_intercept <- mu_intercept[i]

    betad[i, ] <- betad[(i - 1), ]
    for (k in seq_len(pd)) {
      betad_proposed <- rnorm(n = 1, mean = betad[(i - 1), k], sd = sqrt(metrop_var_betad[k]))
      ratio <- exp(bmzdp_lh_betad(betad_proposed, k, dat, st, sigma2d) - bmzdp_lh_betad(betad[(i - 1), k], k, dat, st, sigma2d))
      if (ratio >= runif(n = 1, min = 0, max = 1)) {
        betad[i, k] <- betad_proposed
        betad_accept[k] <- betad_accept[k] + 1
      }
      st$betad[k] <- betad[i, k]
    }

    tao[i] <- tao[(i - 1)]
    tao_proposed <- rnorm(n = 1, mean = tao[(i - 1)], sd = sqrt(metrop_var_tao))
    ratio <- exp(bmzdp_lh_tao(tao_proposed, dat, st, sigma2tao) - bmzdp_lh_tao(tao[(i - 1)], dat, st, sigma2tao))
    if (ratio >= runif(n = 1, min = 0, max = 1)) {
      tao[i] <- tao_proposed
      tao_accept <- tao_accept + 1
    }
    st$tao <- tao[i]

    for (k in seq_len(num_hos)) {
      for (s in seq_len(num_dp_cluster)) {
        st$mu_hos[k] <- b_dp_cluster[(i - 1), s]
        ploglik <- 0
        for (l in members[[k]]) {
          ploglik <- ploglik + bmzdp_loglik_l(l, dat, st)
        }
        pai[s] <- (s_label[(i - 1), k] == s) * log(w_dp_cluster[(i - 1), s]) + ploglik
      }
      pai <- pai - max(pai)
      pai <- exp(pai)/sum(exp(pai))
      s_label[i, k] <- bmzdp_draw_label(pai)
      st$mu_hos[k] <- b_dp_cluster[(i - 1), s_label[i, k]]
    }

    b_dp_cluster[i, ] <- b_dp_cluster[(i - 1), ]
    for (k in seq_len(num_dp_cluster)) {
      in_k <- s_label[i, ] == k
      b_dp_proposed <- rnorm(n = 1, mean = b_dp_cluster[(i - 1), k], sd = sqrt(metrop_var_b_dp_cluster[k]))
      ratio <- exp(bmzdp_lh_b_dp_cluster(b_dp_proposed, in_k, dat, st, v2) - bmzdp_lh_b_dp_cluster(b_dp_cluster[(i - 1), k], in_k, dat, st, v2))
      if (ratio >= runif(n = 1, min = 0, max = 1)) {
        b_dp_cluster[i, k] <- b_dp_proposed
        b_dp_cluster_accept[k] <- b_dp_cluster_accept[k] + 1
      }
      mu_hos[i, in_k] <- b_dp_cluster[i, k]
      st$mu_hos[in_k] <- b_dp_cluster[i, k]
    }

    w_dp_cluster[i, ] <- bmzdp_stick(s_label[i, ], num_dp_cluster, alpha)
  }

  accept <- list(
    betar = betar_accept,
    gamma_sub = gamma_sub_accept,
    eta = eta_accept,
    zeta = zeta_accept,
    mu_intercept = mu_intercept_accept,
    betad = betad_accept,
    tao = tao_accept,
    b_dp_cluster = b_dp_cluster_accept,
    b_dp_cluster_n = b_dp_cluster_accept_n
  )
  if (piecewise) {
    accept$lambda_n <- lambda_n_accept
    base_draws <- list(lambda_n = lambda_n)
  } else {
    accept$kesi <- kesi_accept
    base_draws <- list(kesi = kesi)
  }

  c(
    list(baseline = baseline, betar = betar),
    base_draws,
    list(
      gamma_sub = gamma_sub,
      epislon = epislon,
      eta = eta,
      p = p,
      y = y,
      y_index = y_index,
      zeta = zeta,
      mu_intercept = mu_intercept,
      betad = betad,
      tao = tao,
      mu_hos = mu_hos,
      sigma = sigma,
      s_label = s_label,
      s_label_n = s_label_n,
      b_dp_cluster = b_dp_cluster,
      b_dp_cluster_n = b_dp_cluster_n,
      w_dp_cluster = w_dp_cluster,
      w_dp_cluster_n = w_dp_cluster_n,
      accept = accept,
      metrop_var_gamma_sub = metrop_var_gamma_sub
    )
  )
}
