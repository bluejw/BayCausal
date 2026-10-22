# Run the real RJ functions with a fixed dimension-matching draw and capture
# the acceptance calculation. Choosing the last singleton reverses the split.
probe_zero_rj <- function(fun, args, reject = FALSE) {
  terms <- NULL
  probe_env <- new.env(parent = environment(fun))
  probe_env$sample <- function(x, size, ...) tail(x, size)
  probe_env$runif <- function(n, min = 0, max = 1) {
    if (min == -1) return(0.25)
    frame <- parent.frame()
    terms <<- mget(c("log_proposal_ratio", "log_jacobian", "log_accept_prob"),
                    envir = frame)
    if (reject) 0.5 else 0
  }
  if (reject) {
    calls <- 0L
    bad_likelihood <- function(...) {
      calls <<- calls + 1L
      if (calls == 1L) 0 else -1e6
    }
    probe_env$compute_loglik_rjmcmc_cpp <- bad_likelihood
    probe_env$compute_loglik_rjmcmc_t_rcpp <- bad_likelihood
  }
  environment(fun) <- probe_env
  result <- do.call(fun, args)
  list(result = result, terms = terms)
}

zero_rj_args <- function(k = 0L, H = 3L, ps = 0.5, active = FALSE) {
  Q <- H + 1L
  n <- 6L
  delta <- matrix(0, Q, max(1L, k))
  if (k > 0) delta[cbind(seq_len(k), seq_len(k))] <- 1
  if (active) delta[Q, 1] <- 1
  list(mu = rep(0, Q), A = matrix(0, Q, 1), B = matrix(0, Q, Q),
       L = 0.4 * delta, C = matrix(if (k == 0) 0 else 0.2, n, max(1L, k)),
       sigma2 = rep(0.8, Q), a1 = 1, a2 = 2, sparsity_matrix = delta,
       pivot = seq_len(k), P_star = k, kappa = 1,
       mh_setup_lst = list(ps = ps),
       prior_lst = list(H = H, a_sigma = 1, b_sigma = 1),
       data = list(Y = matrix(seq_len(n * Q) / (n * Q), n, Q),
                   X = matrix(0, n, 1)))
}

zero_rj_next <- function(args, result) {
  for (name in c("L", "C", "sigma2", "sparsity_matrix", "pivot", "P_star")) {
    args[[name]] <- result[[paste0("mcmc_", name)]]
  }
  args
}

test_that("move selection permits the last singleton to merge", {
  zero <- matrix(0, 3, 1)
  singleton <- matrix(c(1, 0, 0), 3, 1)
  active <- matrix(c(1, 1, 0), 3, 1)
  expect_equal(compute_psplit(integer(0), zero, 0.3, 2), 1)
  expect_equal(compute_psplit(1L, singleton, 0.3, 2), 0.3)
  expect_equal(compute_psplit(1L, singleton, 0.3, 1), 0)
  expect_equal(compute_psplit(1L, active, 0.3, 2), 1)
  expect_true(is.na(compute_psplit(1L, active, 0.3, 1)))
})

test_that("paired RJ moves use the actual boundary probabilities", {
  # Expected reverse/forward proposal ratios include the Uniform(-1, 1)
  # density and the existing H-slot and singleton selection factors.
  cases <- list(
    list(args = zero_rj_args(H = 3), ratio = 3),
    list(args = zero_rj_args(H = 1), ratio = 2),
    list(args = zero_rj_args(k = 1, H = 3), ratio = 2),
    list(args = zero_rj_args(k = 1, H = 2), ratio = 2),
    list(args = zero_rj_args(H = 3, ps = 0.3), ratio = 4.2),
    list(args = zero_rj_args(k = 1, H = 3, ps = 0.3), ratio = 14 / 3),
    list(args = zero_rj_args(k = 1, H = 3, ps = 0.3, active = TRUE), ratio = 2.8),
    list(args = zero_rj_args(k = 2, H = 3, ps = 0.3), ratio = 20 / 9)
  )
  for (suffix in c("", "_t")) {
    split <- get(paste0("update_P_star_split", suffix), asNamespace("BayCausal"))
    merge <- get(paste0("update_P_star_merge", suffix), asNamespace("BayCausal"))
    for (case in cases) {
      args <- case$args
      if (suffix == "_t") args$nu_t <- 5
      forward <- probe_zero_rj(split, args)
      reverse <- probe_zero_rj(merge, zero_rj_next(args, forward$result))
      expect_equal(forward$result$mcmc_P_star, args$P_star + 1)
      expect_equal(ncol(forward$result$mcmc_L), args$P_star + 1)
      expect_equal(forward$terms$log_proposal_ratio, log(case$ratio))
      expect_equal(reverse$terms$log_proposal_ratio, -log(case$ratio))
      expect_equal(forward$terms$log_jacobian + reverse$terms$log_jacobian, 0)
      expect_equal(unname(forward$terms$log_accept_prob + reverse$terms$log_accept_prob), 0)
      for (name in c("P_star", "L", "C", "sigma2", "sparsity_matrix", "pivot")) {
        expect_equal(unname(reverse$result[[paste0("mcmc_", name)]]),
                     unname(args[[name]]))
      }
    }
  }
})

test_that("rejected boundary proposals leave the state intact", {
  for (suffix in c("", "_t")) {
    for (move in c("split", "merge")) {
      args <- zero_rj_args(k = if (move == "split") 0L else 1L)
      if (suffix == "_t") args$nu_t <- 5
      fun <- get(paste0("update_P_star_", move, suffix), asNamespace("BayCausal"))
      result <- probe_zero_rj(fun, args, reject = TRUE)$result
      for (name in c("P_star", "L", "C", "sigma2", "sparsity_matrix", "pivot")) {
        expect_identical(result[[paste0("mcmc_", name)]], args[[name]])
      }
    }
  }
})

test_that("post-processing preserves zero-factor bookkeeping", {
  zero <- list(P_star = 0L, L = matrix(0, 3, 1),
               sparsity_matrix = matrix(0, 3, 1), pivot = integer(0))
  singleton <- list(P_star = 1L, L = matrix(c(0.05, 0, 0), 3, 1),
                    sparsity_matrix = matrix(c(1, 0, 0), 3, 1), pivot = 1L)
  expect_identical(post_process_L(list(zero), TRUE, TRUE)[[1]], zero)
  thresholded <- post_process_L(list(singleton), TRUE, TRUE, cutoff = 0.1)
  removed <- post_process_L(list(singleton), TRUE, TRUE, remove_singleton = TRUE)
  expect_equal(thresholded[[1]], zero)
  expect_equal(removed[[1]], zero)
  expect_equal(post_process_L(thresholded, TRUE, TRUE), thresholded)
  expect_equal(extract_P_star(c(list(zero), thresholded)), c(0, 0))
  expect_equal(extract_post_mean(c(list(zero), thresholded), "L", 0), zero$L)
  expect_equal(post_process_L_mean(zero$L, TRUE, TRUE), zero$L)
  expect_equal(post_process_L_mean(singleton$L, TRUE, TRUE, cutoff = 0.1), zero$L)
  expect_equal(post_process_L_mean(singleton$L, TRUE, TRUE,
                                   remove_singleton = TRUE), zero$L)
})

test_that("all RJ sampler entry points run from zero factors", {
  set.seed(92)
  data <- list(Y = matrix(rnorm(60), 20, 3), X = matrix(rnorm(40), 20, 2))
  prior <- set_prior_default(data)
  prior$nu_t <- 5
  mh <- set_mh_default()
  chain <- set_chain(Nit = 30, burn = 0, thin = 1, seed = 12)
  init <- set_init_default(1, 0L, data, prior)
  init_wo_x <- set_init_default_wo_x(1, 0L, data, prior)
  for (state in list(init, init_wo_x)) {
    expect_identical(state$P_star, 0L)
    expect_equal(state$L, matrix(0, 3, 1))
    expect_equal(state$C, matrix(0, 20, 1))
    expect_equal(state$sparsity_matrix, matrix(0, 3, 1))
    expect_identical(state$pivot, integer(0))
  }
  for (name in c("glvcausal", "glvcausal_t", "glvcausal_wo_x",
                 "glvcausal_only_UGLT", "glvcausal_only_UGLT_t",
                 "glvcausal_only_UGLT_wo_x", "lingfa", "lingfa_wo_x")) {
    fun <- get(name, asNamespace("BayCausal"))
    start <- if (grepl("wo_x", name)) init_wo_x else init
    result <- fun(data, mh, start, prior, chain, FALSE)
    samples <- if (grepl("only_UGLT", name)) list(result) else result
    expect_gt(length(samples), 0)
    for (state in samples) {
      expect_true(all(is.finite(state$sigma2)) && all(state$sigma2 > 0))
      expect_equal(ncol(state$L), max(1, state$P_star))
      expect_length(state$pivot, state$P_star)
      if (state$P_star == 0) {
        expect_equal(state$L, matrix(0, 3, 1))
        expect_equal(state$sparsity_matrix, matrix(0, 3, 1))
        if (!is.null(state$C)) expect_equal(state$C, matrix(0, 20, 1))
      }
    }
  }
})

test_that("factor-analysis samplers continue after merging the last factor", {
  set.seed(92)
  data <- list(Y = matrix(rnorm(60), 20, 3), X = matrix(rnorm(40), 20, 2))
  prior <- set_prior_default(data)
  mh <- set_mh_default()
  chain <- set_chain(Nit = 100, burn = 0, thin = 1, seed = 12)
  for (name in c("lingfa", "lingfa_wo_x")) {
    init <- if (name == "lingfa") set_init_default(1, 1L, data, prior) else
      set_init_default_wo_x(1, 1L, data, prior)
    fun <- get(name, asNamespace("BayCausal"))
    samples <- fun(data, mh, init, prior, chain, FALSE)
    expect_length(samples, chain$Nit - 1L)
    counts <- c(init$P_star, extract_P_star(samples))
    expect_true(any(head(counts, -1) == 1 & tail(counts, -1) == 0))
    expect_true(any(head(counts, -1) == 0 & tail(counts, -1) == 1))
    zero_samples <- Filter(function(state) state$P_star == 0, samples)
    expect_true(all(vapply(zero_samples, function(state) {
      identical(dim(state$L), c(3L, 1L)) && all(state$L == 0) &&
        all(state$sparsity_matrix == 0) && length(state$pivot) == 0
    }, logical(1))))
  }
})
