PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename $(PWD))
TGZ     := $(PKGSRC)_$(PKGVERS).tar.gz
TGZVNR  := $(PKGSRC)_$(PKGVERS)-vignettes-not-rebuilt.tar.gz
WINBIN  := $(PKGSRC)_$(PKGVERS).zip

# Specify the directory holding R binaries. To use an alternate R build (say a
# pre-prelease version) use `make RBIN=/path/to/other/R/` or `export RBIN=...`
# If no alternate bin folder is specified, the default is to use the folder
# containing the first instance of R on the PATH.
RBIN ?= $(shell dirname "`which R`")

# Specify package and static documentation directories for subversion on r-forge
RFSVN ?= $(HOME)/svn/kinfit.r-forge
RFDIR ?= $(RFSVN)/pkg/mkin
SDDIR ?= $(RFSVN)/www/mkin_static

# Vignettes are listed in the build target
pkgfiles = \
	.Rbuildignore \
	data/* \
	DESCRIPTION \
	man/* \
	NAMESPACE \
	NEWS.md \
	R/* \
	tests/* \
	tests/testthat* \
	TODO

all: build

$(TGZ): $(pkgfiles) vignettes
	$(RM) -r vignettes/*_cache
	$(RM) -r vignettes/*_files
	$(RM) -r vignettes/*.R
	$(RM) -r vignettes/web_only/*.R
	$(RM) Rplots.pdf
	"$(RBIN)/R" CMD build . 2>&1 | tee build.log

roxygen: 
	"$(RBIN)/Rscript" -e 'devtools::document()'

$(TGZVNR): $(pkgfiles)
	"$(RBIN)/R" CMD build . --no-build-vignettes;\
	mv $(TGZ) $(TGZVNR)

build: roxygen $(TGZ)

build-no-vignettes: $(TGZVNR)

install: build
	"$(RBIN)/R" CMD INSTALL $(TGZ)

quickinstall: build-no-vignettes
	"$(RBIN)/R" CMD INSTALL $(TGZVNR)

check: roxygen build
	_R_CHECK_CRAN_INCOMING_REMOTE_=false "$(RBIN)/R" CMD check --as-cran --no-tests $(TGZ) 2>&1 | tee check.log

quickcheck: roxygen build-no-vignettes
	mv $(TGZVNR) $(TGZ)
	"$(RBIN)/R" CMD check --no-tests --no-build-vignettes --no-vignettes $(TGZ)
	mv $(TGZ) $(TGZVNR)

clean:
	$(RM) -r vignettes/*_cache
	$(RM) -r vignettes/*_files
	$(RM) -r vignettes/*.R
	$(RM) -r vignettes/web_only/*.R
	$(RM) Rplots.pdf

test: install
	"$(RBIN)/Rscript" -e 'devtools::test()' 2>&1 | tee test.log
	sed -i -e "s/\r.*\r//" test.log

slowtests: install
	NOT_CRAN=true "$(RBIN)/Rscript" -e 'library(mkin); testthat::test_dir("tests/testthat/slow")' 2>&1 | tee tests_slow.log
	sed -i -e "s/\r.*\r//" tests_slow.log

vdiffr:
	"$(RBIN)/Rscript" -e 'vdiffr::manage_cases(filter = "plots|nafta")'

testcheck: test check

README.html: README.md
	"$(RBIN)/Rscript" -e "rmarkdown::render('README.md', output_format = 'html_document', output_options = list(mathjax = NULL))"

vignettes/%.html: vignettes/mkin_vignettes.css vignettes/references.bib vignettes/%.Rmd
	"$(RBIN)/Rscript" -e "tools::buildVignette(file = 'vignettes/$*.Rmd', dir = 'vignettes')"

vignettes: vignettes/mkin.html vignettes/FOCUS_D.html vignettes/FOCUS_L.html vignettes/twa.html

vignettes/web_only/%.html: vignettes/references.bib vignettes/web_only/%.Rmd
	"$(RBIN)/Rscript" -e "tools::buildVignette(file = 'vignettes/web_only/$*.Rmd', dir = 'vignettes/web_only')"

articles: vignettes/web_only/FOCUS_Z.html vignettes/web_only/compiled_models.html

pd:
	"$(RBIN)/Rscript" -e "pkgdown::build_site(run_dont_run = TRUE, lazy = TRUE)"
	git add -A
	git commit -m 'Static documentation rebuilt by pkgdown' -e

r-forge:
	git archive master > $(HOME)/git/mkin/mkin.tar;\
	cd $(RFDIR) && rm -r `ls` && tar -xf $(HOME)/git/mkin/mkin.tar;\
	rm -r $(SDDIR)/*;\
	cp -a docs/* $(SDDIR);\
	svn add --force .; svn rm --force `svn status | grep "\!" | cut -d " " -f 8`; cd $(RFSVN) && svn commit -m 'sync with git'

winbuilder: build
	date
	@echo "Uploading to R-release on win-builder"
	curl -T $(TGZ) ftp://anonymous@win-builder.r-project.org/R-release/
	@echo "Uploading to R-devel on win-builder"
	curl -T $(TGZ) ftp://anonymous@win-builder.r-project.org/R-devel/

drat: build
	"$(RBIN)/Rscript" -e "drat::insertPackage('$(TGZ)', commit = TRUE)"

$(WINBIN): build
	@echo "Building windows binary package..."
	"$(RBIN)/R" CMD INSTALL $(TGZ) --build
	@echo "DONE."

winbin: $(WINBIN)

dratwin: winbin
	"$(RBIN)/Rscript" -e "drat::insertPackage('$(WINBIN)', '~/git/drat/', commit = TRUE)"

submit:
	@echo "\nHow about make test, make check, make pd, make winbuilder"
	@echo "\nIs the DESCRIPTION file up to date?"
	@echo "\nIs the NEWS.md file up to date?"
	@echo "\nAre you sure you want to release to CRAN?"
	@echo "\nThen make r-forge, commit leftover changes if any, tag the release and use the form at http://cran.r-project.org/submit.html"
