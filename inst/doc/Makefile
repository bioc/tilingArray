all: findsegments costMatrix assessNorm segmentation clean

findsegments: findsegments.tex
	pdflatex findsegments
	pdflatex findsegments

costMatrix: costMatrix.tex
	pdflatex costMatrix
	pdflatex costMatrix

assessNorm: assessNorm.tex
	cp -p ../scripts/assessNorm.pdf .

segmentation: segmentation.tex
	cp -p ../scripts/segmentation.pdf .

clean:
	rm -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	rm -f *.dvi *.ps findsegments-* costMatrix-* Rplots.pdf


