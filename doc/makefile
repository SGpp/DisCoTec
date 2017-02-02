all:	gene.pdf style_guide.pdf tutorial.pdf clean

.SUFFIXES:
.SUFFIXES: .tex .dvi .pdf .ps .txt

.tex.pdf:
	((pdflatex $< && pdflatex $<) || \
	(latex $< && latex $<))

.dvi.ps:
	dvips $<

gene.pdf: gene.tex gene-faq.tex

style_guide.pdf: style_guide.tex

tutorial.pdf: tutorial.tex

clean:
	rm -f *.dvi *.aux *.log *~ *.toc *.nav *.vrb *.snm

distclean:  clean
	rm -f gene.pdf style_guide.pdf tutorial.pdf
