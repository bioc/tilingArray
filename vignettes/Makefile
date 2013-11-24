# vignettes actually build
BUILD = findsegments costMatrix plotAlongChrom

# prebuild vignettes copied from the SRCS_DIR directory
STUBS = assessNorm segmentation
SRCS_DIR = ../inst/vignettes

all: pdf copy

pdf: $(addsuffix .pdf,$(BUILD))

copy: $(STUBS)

%.pdf: %.tex
	pdflatex $*
	pdflatex $*

%: %.tex
	cp -p $(SRCS_DIR)/$*.pdf .

clean:
	$(RM) *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	$(RM) $(addsuffix -*,$(BUILD))
