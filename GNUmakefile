PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename $(PWD))
TGZ     := $(PKGSRC)_$(PKGVERS).tar.gz
TGZVNR  := $(PKGSRC)_$(PKGVERS)-vignettes-not-rebuilt.tar.gz

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
	"$(RBIN)/R" CMD build . 2>&1 | tee build.log

$(TGZVNR): $(pkgfiles)
	"$(RBIN)/R" CMD build . --no-build-vignettes;\
	mv $(TGZ) $(TGZVNR)

build: $(TGZ)

build-no-vignettes: $(TGZVNR)

install: build
	"$(RBIN)/R" CMD INSTALL $(TGZ)

quickinstall: build-no-vignettes
	"$(RBIN)/R" CMD INSTALL $(TGZVNR)

check: build
	"$(RBIN)/R" CMD check --as-cran --no-tests $(TGZ) 2>&1 | tee check.log

quickcheck: build-no-vignettes
	mv $(TGZVNR) $(TGZ)
	"$(RBIN)/R" CMD check --no-tests --no-build-vignettes --no-vignettes $(TGZ)
	mv $(TGZ) $(TGZVNR)

clean:
	$(RM) -r vignettes/*_cache
	$(RM) -r vignettes/*_files
	$(RM) -r vignettes/*.R
	$(RM) Rplots.pdf

test: quickinstall
	NOT_CRAN=true "$(RBIN)/Rscript" -e 'devtools::test()' 2>&1 | tee test.log

README.html: README.md
	"$(RBIN)/Rscript" -e "rmarkdown::render('README.md', output_format = 'html_document', output_options = list(mathjax = NULL))"

vignettes/%.html: vignettes/mkin_vignettes.css vignettes/references.bib vignettes/%.Rmd
	"$(RBIN)/Rscript" -e "tools::buildVignette(file = 'vignettes/$*.Rmd', dir = 'vignettes')"

vignettes: vignettes/mkin.html vignettes/FOCUS_D.html vignettes/FOCUS_L.html vignettes/FOCUS_Z.html vignettes/compiled_models.html

pd:
	"$(RBIN)/Rscript" -e "pkgdown::build_site()"
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

submit:
	@echo "\nHow about make test, make check, make pd, make winbuilder"
	@echo "\nIs the DESCRIPTION file up to date?"
	@echo "\nIs the NEWS.md file up to date?"
	@echo "\nAre you sure you want to release to CRAN?"
	@echo "\nThen make r-forge, commit leftover changes if any, tag the release and use the form at http://cran.r-project.org/submit.html"
