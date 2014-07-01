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

# Specify the directory where the static documentation belongs
SDDIR ?= $(HOME)/svn/kinfit.r-forge/www/mkin_static

pkgfiles = NEWS \
	   data/* \
	   DESCRIPTION \
	   inst/unitTests* \
	   inst/staticdocs/README \
	   man/* \
	   NAMESPACE \
	   R/* \
	   README.md \
	   tests/* \
	   TODO \
	   vignettes/*

all: NEWS check clean

# convert markdown to R's NEWS format (from knitr package)
NEWS: NEWS.md
	sed -e 's/^-/ -/' -e 's/^## *//' -e 's/^#/\t\t/' <NEWS.md | fmt -80 >NEWS

$(TGZ): $(pkgfiles)
	cd ..;\
		"$(RBIN)/R" CMD build $(PKGSRC)

$(TGZVNR): $(pkgfiles)
	cd ..;\
		"$(RBIN)/R" CMD build $(PKGSRC) --no-build-vignettes;\
		cd $(PKGSRC);\
	mv $(TGZ) $(TGZVNR)
                
build: $(TGZ)

build-no-vignettes: $(TGZVNR)

install: build
	"$(RBIN)/R" CMD INSTALL $(TGZ)

install-no-vignettes: build-no-vignettes
	"$(RBIN)/R" CMD INSTALL $(TGZVNR)

check: build
	# Vignettes have been rebuilt by the build target
	"$(RBIN)/R" CMD check --as-cran --no-tests --no-build-vignettes $(TGZ)

check-no-vignettes: build-no-vignettes
	mv $(TGZVNR) $(TGZ)
	"$(RBIN)/R" CMD check --as-cran --no-tests $(TGZ)
	mv $(TGZ) $(TGZVNR)

clean: 
	$(RM) -r $(PKGNAME).Rcheck/

test: install-no-vignettes
	cd tests;\
		"$(RBIN)/Rscript" doRUnit.R

.PHONY: vignettes
vignettes:
	"$(RBIN)/Rscript" -e "tools::buildVignettes(dir = '.')"
		
sd:
	"$(RBIN)/Rscript" -e "library(staticdocs); build_site()"

move-sd:
	rm -rf $(SDDIR)/*;\
		cp -r inst/web/* $(SDDIR)

winbuilder: build
	date
	@echo "Uploading to R-release on win-builder"
	curl -T $(TGZ) ftp://anonymous@win-builder.r-project.org/R-release/
	@echo "Uploading to R-devel on win-builder"
	curl -T $(TGZ) ftp://anonymous@win-builder.r-project.org/R-devel/

r-forge:
	@echo "\nHow about make test and make check?"
	@echo "\nIs the DESCRIPTION file up to date?"
	@echo "\nTo update the svn repository tied to the local r-forge branch with"
	@echo "changes in the local master branch, run:"
	@echo "'git checkout r-forge'"
	@echo "'git merge --squash -srecursive -Xtheirs --no-commit master'"
	@echo "'git commit'"
	@echo "'git svn dcommit'"
	@echo "\nThen change back to the master branch:"
	@echo "'git checkout master'"

submit:
	@echo "\nAre you sure you want to release to CRAN?"
	@echo "\nThen use the form at http://cran.r-project.org/submit.html"
