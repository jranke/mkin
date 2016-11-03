library(pkgdown)

debug(pkgdown:::as_pkgdown)
debug(pkgdown:::topic_index)
debug(pkgdown:::package_rd) # creates rd object (list of Rd file representations)


undebug(pkgdown:::rd_file) # creates rd file representations

build_reference_index(".")

extract_tag(rd[[1]], "tag_title") # gives "\n"
extract_tag(rd[[4]], "tag_title") # gives the title

library(magrittr)
library(stringr)
extract_tag_local <- function(x, tag) {
  x %>%
    purrr::keep(inherits, tag) %>%
    unlist %>%
    paste(collapse = "") %>%
    str_trim 
}

extract_tag_local(rd[[1]], "tag_title")
extract_tag_local(rd[[4]], "tag_title")


titles <- purrr::map_chr(rd, extract_tag, "tag_title")

titles <- purrr::map_chr(rd, extract_tag_local, "tag_title")

rd[[2]]
