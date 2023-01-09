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
RDEVBIN=/home/jranke/svn/R/r-devel/build/bin

# Specify package and static documentation directories for subversion on r-forge
RFSVN ?= $(HOME)/svn/r-forge/kinfit
RFDIR ?= $(RFSVN)/pkg/mkin
SDDIR ?= $(RFSVN)/www/mkin_static

# Vignettes are listed in the build target
pkgfiles = \
	.Rbuildignore \
	data/* \
	DESCRIPTION \
	inst/WORDLIST \
	inst/dataset_generation/* \
	inst/rmarkdown/templates/hierarchical_kinetics/template.yaml \
	inst/rmarkdown/templates/hierarchical_kinetics/skeleton/skeleton.Rmd \
	inst/testdata/* \
	man/* \
	NAMESPACE \
	NEWS.md \
	R/* \
	tests/* \
	tests/testthat*

all: build

$(TGZ): $(pkgfiles) vignettes
	$(RM) Rplots.pdf
	"$(RBIN)/R" CMD build . 2>&1 | tee log/build.log

roxygen:
	"$(RBIN)/Rscript" -e 'devtools::document()'

$(TGZVNR): $(pkgfiles)
	"$(RBIN)/R" CMD build . --no-build-vignettes;\
	mv $(TGZ) $(TGZVNR)

build: roxygen $(TGZ)

build-no-vignettes: $(TGZVNR)

install: build
	"$(RBIN)/R" CMD INSTALL $(TGZ)

devinstall: build
	"$(RDEVBIN)/R" CMD INSTALL $(TGZ)

quickinstall: build-no-vignettes
	"$(RBIN)/R" CMD INSTALL $(TGZVNR)

check: roxygen build
	_R_CHECK_CRAN_INCOMING_REMOTE_=false "$(RBIN)/R" CMD check --as-cran --no-tests $(TGZ) 2>&1 | tee log/check.log

devcheck: roxygen build
	_R_CHECK_CRAN_INCOMING_REMOTE_=false "$(RDEVBIN)/R" CMD check --as-cran --no-tests $(TGZ) 2>&1 | tee log/check_dev.log

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

# We set PROCESSX_NOTIFY_OLD_SIGCHILD in order to avoid the message
# "Error while shutting down parallel: unable to terminate some child processes",
# which is said to be harmless, see https://processx.r-lib.org/#mixing-processx-and-the-parallel-base-r-package
# and https://github.com/r-lib/processx/issues/236
test: install
	PROCESSX_NOTIFY_OLD_SIGCHLD=true "$(RBIN)/Rscript" -e 'options(cli.dynamic = TRUE); devtools::test()' 2>&1 | tee log/test.log
	sed -i -e "s/.*\r.*\r//" log/test.log

devtest: devinstall
	PROCESSX_NOTIFY_OLD_SIGCHLD=true "$(RDEVBIN)/Rscript" -e 'options(cli.dynamic = TRUE); devtools::test()' 2>&1 | tee log/test_dev.log
	sed -i -e "s/\r.*\r//" log/test_dev.log

slowtests: install
	NOT_CRAN=true "$(RBIN)/Rscript" -e 'cli.dynamic = TRUE); library(mkin); testthat::test_dir("tests/testthat/slow")' 2>&1 | tee log/tests_slow.log
	sed -i -e "s/\r.*\r//" log/tests_slow.log

testcheck: roxygen test check

README.html: README.md
	"$(RBIN)/Rscript" -e "rmarkdown::render('README.md', output_format = 'html_document', output_options = list(mathjax = NULL))"

vignettes/%.html: vignettes/mkin_vignettes.css vignettes/references.bib vignettes/%.rmd
	"$(RBIN)/Rscript" -e "tools::buildVignette(file = 'vignettes/$*.rmd', dir = 'vignettes')"

vignettes: vignettes/mkin.html vignettes/FOCUS_D.html vignettes/FOCUS_L.html vignettes/twa.html

vignettes/web_only/%.html: vignettes/references.bib vignettes/web_only/%.rmd
	"$(RBIN)/Rscript" -e "tools::buildVignette(file = 'vignettes/web_only/$*.rmd', dir = 'vignettes/web_only', keep=c('mkin_benchmarks.rda', 'saem_benchmarks.rda'))"

vignettes/prebuilt/%.pdf: vignettes/prebuilt/references.bib vignettes/prebuilt/%.rmd
	"$(RBIN)/Rscript" -e "rmarkdown::render('vignettes/prebuilt/$*.rmd')"

pd: roxygen
	"$(RBIN)/Rscript" -e "pkgdown::build_site(run_dont_run = TRUE, lazy = TRUE)"
	git add -A

pd_all: roxygen
	"$(RBIN)/Rscript" -e "pkgdown::build_site(run_dont_run = TRUE)"
	git add -A

r-forge:
	git archive main > $(HOME)/git/mkin/mkin.tar;\
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
