library(pkgdown)
library(testthat)

debug(pkgdown:::topic_index)

build_reference_index(".")

extract_tag(rd[[1]], "tag_title") # gives "\n"
extract_tag(rd[[4]], "tag_title") # gives the title

extract_tag(rd[[1]], "tag_alias") # one aliase
extract_tag(rd[[13]], "tag_alias") # two aliases

library(magrittr)

extract_tag_1 <- function(x, tag) {
  x %>%
    purrr::keep(inherits, tag) %>%
    unlist %>%
    paste(collapse = "") %>%
    trimws
}

extract_tag <- function(x, tag) {
  x %>%
    purrr::keep(inherits, tag) %>%
    purrr::map_chr(c(1, 1))
}

extract_tag_1 <- function(x, tag) {
  x %>%
    purrr::keep(inherits, tag) %>%
    purrr::map_chr(function(x) trimws(paste(x, collapse = " ")))
}

subset(rd[[1]][[1]], 

sapply(extract_tag_1(rd[[1]], "tag_title"), function(x) trimws(paste(x, collapse = " ")))

extract_tag_1(rd[[1]], "tag_title")
extract_tag_1(rd[[4]], "tag_title")

extract_tag_1(rd[[1]], "tag_alias") # one aliase
extract_tag_1(rd[[13]], "tag_alias") # two aliases

