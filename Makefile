all: derivations.pdf

derivations.pdf: derivations.tex
	pdflatex derivations.tex
	pdflatex derivations.tex

clean:
	rm *.aux *.log *.tex~