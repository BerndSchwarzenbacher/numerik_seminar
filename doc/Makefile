all: doc.pdf presentation.pdf

doc.pdf: doc.tex doc.bib
	pdflatex doc.tex doc.pdf
	bibtex doc
	pdflatex doc.tex doc.pdf
	pdflatex doc.tex doc.pdf

presentation.pdf: presentation.tex
	pdflatex presentation.tex presentation.pdf
	pdflatex presentation.tex presentation.pdf

clean:
	rm *.log *.aux *.pdf *.toc *.vrb *.nav *.snm *.out *.bbl *.blg

