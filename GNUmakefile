PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename $(PWD))

# Specify the directory holding R binaries. To use an alternate R build (say a
# pre-prelease version) use `make RBIN=/path/to/other/R/` or `export RBIN=...`
# If no alternate bin folder is specified, the default is to use the folder
# containing the first instance of R on the PATH.
RBIN ?= $(shell dirname "`which R`")

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
	@echo "  release    Give some reminders"
	@echo ""
	@echo "Using R in: $(RBIN)"
	@echo "Set the RBIN environment variable to change this."
	@echo ""


#------------------------------------------------------------------------------
# Development Tasks
#------------------------------------------------------------------------------

build:
	cd ..;\
		"$(RBIN)/R" CMD build $(PKGSRC)

build-no-vignettes:
	cd ..;\
		"$(RBIN)/R" CMD build $(PKGSRC) --no-build-vignettes

install: build
	cd ..;\
		"$(RBIN)/R" CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

install-no-vignettes: build-no-vignettes
	cd ..;\
		"$(RBIN)/R" CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
		"$(RBIN)/R" CMD check --as-cran --no-tests $(PKGNAME)_$(PKGVERS).tar.gz

check-no-vignettes: build-no-vignettes
	cd ..;\
		"$(RBIN)/R" CMD check --as-cran --no-tests $(PKGNAME)_$(PKGVERS).tar.gz

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
release:
	@echo "\nHow about make test and make check?"
	@echo "\nIs the DESCRIPTION file up to date?"
	@echo "\nTo update the svn repository tied to the local r-forge branch with"
	@echo "changes in the local master branch, run:"
	@echo "'git checkout r-forge'"
	@echo "'git merge --squash master'"
	@echo "'git commit'"
	@echo "'git svn dcommit'"
	@echo "\nThen change back to the master branch:"
	@echo "'git checkout master'"
