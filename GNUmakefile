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

.PHONY: help

help:
	@echo "\nExecute development tasks for $(PKGNAME)\n"
	@echo "Usage: \`make <task>\` where <task> is one of:"
	@echo ""
	@echo "Development Tasks"
	@echo "-----------------"
	@echo "  build                   Create the package"
	@echo "  build-no-vignettes      Create the package without rebuilding vignettes"
	@echo "  check                   Invoke build and then check the package"
	@echo "  check-no-vignettes      Invoke build without rebuilding vignettes, and then check"
	@echo "  install                 Invoke build and then install the result"
	@echo "  install-no-vignettes    Invoke build without rebuilding vignettes and then install the result"
	@echo "  test                    Install a new copy of the package without vignette rebuilding"
	@echo "                          and run it through the testsuite"
	@echo "  sd                      Build the static documentation"
	@echo "  move-sd                 Move the static documentation where it belongs"
	@echo ""
	@echo "Packaging Tasks"
	@echo "---------------"
	@echo "  winbuilder              Check building on Windows using the winbuilder service"
	@echo "  r-forge                 Give reminders how to sync the r-forge repo"
	@echo "  submit                  Submit to CRAN"
	@echo ""
	@echo "Using R in: $(RBIN)"
	@echo "Set the RBIN environment variable to change this."
	@echo ""

#------------------------------------------------------------------------------
# Development Tasks
#------------------------------------------------------------------------------

# These must be manually kept up to date 
pkgfiles = ChangeLog \
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
	# Vignettes have been rebuilt by the build target so do not repeat that here
	"$(RBIN)/R" CMD check --as-cran --no-tests --no-build-vignettes $(TGZ)

check-no-vignettes: build-no-vignettes
	"$(RBIN)/R" CMD check --as-cran --no-tests $(TGZVNR)

test: install-no-vignettes
	cd tests;\
		"$(RBIN)/Rscript" doRUnit.R

sd:
	"$(RBIN)/Rscript" -e "library(staticdocs); build_site()"

move-sd:
	rm -rf $(SDDIR)/*;\
		cp -r inst/web/* $(SDDIR)

#------------------------------------------------------------------------------
# Packaging Tasks
#------------------------------------------------------------------------------
winbuilder: build
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
