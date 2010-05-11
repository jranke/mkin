R.exe -e "Sweave('mkin.Rnw', stylepath=FALSE)"
pdflatex.exe mkin
bibtex.exe mkin
pdflatex.exe mkin
pdflatex.exe mkin
