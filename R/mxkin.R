#' Estimation of parameter distributions from mmkin objects
#'
#' @param object An mmkin row object containing several fits of the same model to different datasets
#' @return A fitted object of class 'mrkin'
#' @examples
#' sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)

m_SFO <- mkinmod(parent = mkinsub("SFO"))
d_SFO_1 <- mkinpredict(m_SFO,
  c(k_parent_sink = 0.1),
  c(parent = 98), sampling_times)
d_SFO_1_long <- mkin_wide_to_long(d_SFO_1, time = "time")
d_SFO_2 <- mkinpredict(m_SFO,
  c(k_parent_sink = 0.05),
  c(parent = 102), sampling_times)
d_SFO_2_long <- mkin_wide_to_long(d_SFO_2, time = "time")
d_SFO_3 <- mkinpredict(m_SFO,
  c(k_parent_sink = 0.02),
  c(parent = 103), sampling_times)
d_SFO_3_long <- mkin_wide_to_long(d_SFO_3, time = "time")

d1 <- add_err(d_SFO_1, function(value) 3, n = 1)
d2 <- add_err(d_SFO_2, function(value) 2, n = 1)
d3 <- add_err(d_SFO_3, function(value) 4, n = 1)
ds <- c(d1 = d1, d2 = d2, d3 = d3)

f <- mmkin("SFO", ds)
x <- mrkin(f)
as.numeric(x)

#'
#' @export
mrkin <- function(object) {
  if (nrow(object) > 1) stop("Only row objects allowed")
  n_d <- ncol(object)
  p_names <- names(parms(object[[1, 1]]))
  p_names_trans <- names(parms(object[[1, 1]]))

  p_mat_start_trans <- sapply(object, parms, transformed = TRUE)
  colnames(p_mat_start_trans) <- colnames(object)
  p_mat_attr_trans <- attributes(p_mat_start_trans)

  p_dist_names <- p_names[grepl("^log_", p_names)]
  p_free_names <- p_names[!grepl("^log_", p_names)]

  cost <- function(P) {
    p_cost_mat <- P
    attributes(P) <- p_mat_attr_trans

    ll_ds <- 0
    for (i_d in 1:n_d) {
      
    }

  }

  p_mat_start
}
