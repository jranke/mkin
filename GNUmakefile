PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename $(PWD))
TGZ     := ../$(PKGSRC)_$(PKGVERS).tar.gz
TGZVNR  := ../$(PKGSRC)_$(PKGVERS)-vignettes-not-rebuilt.tar.gz

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
	README.md \
	tests/* \
	tests/testthat* \
	TODO

all: build

$(TGZ): $(pkgfiles)
	cd ..;\
		"$(RBIN)/R" CMD build $(PKGSRC) 2>&1 | tee $(PKGNAME)/build.log

$(TGZVNR): $(pkgfiles)
	cd ..;\
		"$(RBIN)/R" CMD build $(PKGSRC) --no-build-vignettes;\
		cd $(PKGSRC);\
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
	$(RM) -r $(PKGNAME).Rcheck/
	$(RM) -r vignettes/*.bbl
	$(RM) -r vignettes/*.blg
	$(RM) -r vignettes/*.fls
	$(RM) -r vignettes/*.fdb_latexmk
	$(RM) -r vignettes/cache
	$(RM) -r vignettes/figure
	$(RM) -r vignettes/*_cache
	$(RM) -r vignettes/*_files
	$(RM) -r vignettes/*-concordance.tex
	$(RM) -r vignettes/*.synctex.gz
	$(RM) Rplots.pdf

test: quickinstall
	cd tests;\
		"$(RBIN)/Rscript" testthat.R 2>&1 | tee ../test.log

vignettes/%.pdf: vignettes/header.tex vignettes/references.bib vignettes/%.Rnw
	"$(RBIN)/Rscript" -e "tools::buildVignette(file = 'vignettes/$*.Rnw', dir = 'vignettes')"

vignettes/%.html: vignettes/mkin_vignettes.css vignettes/%.Rmd
	"$(RBIN)/Rscript" -e "tools::buildVignette(file = 'vignettes/$*.Rmd', dir = 'vignettes')"

vignettes: vignettes/mkin.pdf vignettes/FOCUS_D.html vignettes/FOCUS_L.html vignettes/FOCUS_Z.pdf vignettes/compiled_models.html

sd: install
	rm -rf $(SDDIR)/*
	@echo Now execute
	@echo "\n  library(staticdocs); build_site(site_path = '$(SDDIR)')\n"
	$(RBIN)/R

sd2: install
	rm -rf $(SDDIR)/*
	xvfb-run $(RBIN)/R -e "library(staticdocs); build_site(site_path = '$(SDDIR)')"

r-forge: sd
	cd $(SDDIR) && svn add --force .
	git archive master > $(HOME)/mkin.tar;\
	cd $(RFDIR) && rm -r `ls` && tar -xf $(HOME)/mkin.tar;\
	svn add --force .; svn rm --force `svn status | grep "\!" | cut -d " " -f 8`; cd $(RFSVN) && svn commit -m 'sync with git'
	git add -A
	git commit -m 'Vignettes rebuilt by staticdocs::build_site() for static documentation on r-forge'

winbuilder: build
	date
	@echo "Uploading to R-release on win-builder"
	curl -T $(TGZ) ftp://anonymous@win-builder.r-project.org/R-release/
	@echo "Uploading to R-devel on win-builder"
	curl -T $(TGZ) ftp://anonymous@win-builder.r-project.org/R-devel/

submit:
	@echo "\nHow about make test, make check, make winbuilder"
	@echo "\nIs the DESCRIPTION file up to date?"
	@echo "\nIs the NEWS.md file up to date?"
	@echo "\nAre you sure you want to release to CRAN?"
	@echo "\nThen make r-forge, commit leftover changes if any, tag the release and use the form at http://cran.r-project.org/submit.html"
