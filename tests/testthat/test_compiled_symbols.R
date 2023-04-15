context("Use of precompiled symbols in mkinpredict")

test_that("We can safely use compiled code", {

  # Generate temporary DLL
  sfo_sfo_tmp <- mkinmod(DMTA = mkinsub("SFO", to = "M23"),
    M23 = mkinsub("SFO"))

  # Generate temporary DLL and move to user specified location
  if (!dir.exists("test_dlls")) dir.create("test_dlls")
  sfo_sfo_dll <- mkinmod(DMTA = mkinsub("SFO", to = "M23"),
    M23 = mkinsub("SFO"),
    dll_dir = "test_dlls",
    name = "sfo_sfo",
    unload = TRUE, overwrite = TRUE
  )

  if (Sys.info()["sysname"] != "Windows") {
    # mclapply using forks
    expect_known_output(
      mmkin(list(sfo_sfo_dll), dmta_ds, cores = n_cores, quiet = TRUE),
      "print_mmkin_sfo_sfo_dmta.txt"
    )

    # cluster describing itself as socket cluster
    cl_fork <- parallel::makeForkCluster(n_cores)
    expect_known_output(
      mmkin(list(sfo_sfo_tmp), dmta_ds, cluster = cl_fork, quiet = TRUE),
      "print_mmkin_sfo_sfo_dmta.txt"
    )
    expect_known_output(
      mmkin(list(sfo_sfo_dll), dmta_ds, cluster = cl_fork, quiet = TRUE),
      "print_mmkin_sfo_sfo_dmta.txt"
    )
    parallel::stopCluster(cl_fork)
  }

  # PSOCK cluster
  cl_psock <- parallel::makePSOCKcluster(n_cores)
  expect_known_output(
    mmkin(list(sfo_sfo_tmp), dmta_ds, cluster = cl_psock, quiet = TRUE),
    "print_mmkin_sfo_sfo_dmta.txt"
  )
  expect_known_output(
    mmkin(list(sfo_sfo_dll), dmta_ds, cluster = cl_psock, quiet = TRUE),
    "print_mmkin_sfo_sfo_dmta.txt"
  )
  parallel::stopCluster(cl_psock)

  # Clean up
  expect_true(file.remove("test_dlls/sfo_sfo.so"))
  expect_true(file.remove("test_dlls"))
})

