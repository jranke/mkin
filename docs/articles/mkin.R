## ---- echo = TRUE, cache = TRUE, fig = TRUE, fig.width = 8, fig.height = 7----
require(mkin)
m_SFO_SFO_SFO <- mkinmod(parent = mkinsub("SFO", "M1"),
                         M1 = mkinsub("SFO", "M2"),
                         M2 = mkinsub("SFO"), 
                         use_of_ff = "max", quiet = TRUE)

sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)

d_SFO_SFO_SFO <- mkinpredict(m_SFO_SFO_SFO,
  c(k_parent = 0.03, 
    f_parent_to_M1 = 0.5, k_M1 = log(2)/100, 
    f_M1_to_M2 = 0.9, k_M2 = log(2)/50),
  c(parent = 100, M1 = 0, M2 = 0),
  sampling_times)

d_SFO_SFO_SFO_err <- add_err(d_SFO_SFO_SFO, function(x) 3, n = 1, seed = 123456789 )

f_SFO_SFO_SFO <- mkinfit(m_SFO_SFO_SFO, d_SFO_SFO_SFO_err[[1]], quiet = TRUE)

plot_sep(f_SFO_SFO_SFO, lpos = c("topright", "bottomright", "bottomright"))

