

filename=docu_ssimp

pdflatex ${filename}.tex
bibtex ${filename}.aux
pdflatex ${filename}.tex
pdflatex ${filename}.tex

rm -f ${filename}.{ps,log,aux,out,dvi,bbl,blg,toc}
