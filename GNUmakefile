PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename $(PWD))

# Specify the directory holding R binaries. To use an alternate R build (say a
# pre-prelease version) use `make RBIN=/path/to/other/R/` or `export RBIN=...`
# If no alternate bin folder is specified, the default is to use the folder
# containing the first instance of R on the PATH.
RBIN ?= $(shell dirname "`which R`")


.PHONY: help

help:
	@echo "\nExecute development tasks for $(PKGNAME)\n"
	@echo "Usage: \`make <task>\` where <task> is one of:"
	@echo ""
	@echo "Development Tasks"
	@echo "-----------------"
	@echo "  build      Invoke docs and then create a package"
	@echo "  check      Invoke build and then check the package"
	@echo "  install    Invoke build and then install the result"
	@echo "  test       Install a new copy of the package and run it "
	@echo "             through the testsuite"
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

install: build
	cd ..;\
		"$(RBIN)/R" CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
		"$(RBIN)/R" CMD check --as-cran --no-tests $(PKGNAME)_$(PKGVERS).tar.gz

test: install
	cd tests;\
		"$(RBIN)/Rscript" doRUnit.R

#------------------------------------------------------------------------------
# Packaging Tasks
#------------------------------------------------------------------------------
release:
	@echo "\nPull in changes from svn and merge local commits"
	@git svn rebase
	@echo "\nHow about make test, make check and make vignette?"
	@echo "\nIs the DESCRIPTION file up to date?"
	@echo "\nPerform final changes and commit with 'git commit --amend'."
	@echo "\nIf the above is taken care of, run 'git svn dcommit'"
	@echo "and then 'git push origin master'"
