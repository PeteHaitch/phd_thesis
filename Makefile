outline.pdf: bibtex/phd_thesis.bib markdown/outline.md markdown/bibliography.md
	mkdir -p pdf
	pandoc -s --number-sections --bibliography=bibtex/phd_thesis.bib -o pdf/outline.pdf markdown/outline.md markdown/bibliography.md

preamble.pdf: bibtex/phd_thesis.bib markdown/preamble.md markdown/bibliography.md
	mkdir -p pdf
	pandoc -s --table-of-contents --number-sections --bibliography=bibtex/phd_thesis.bib -o pdf/preamble.pdf markdown/preamble.md markdown/bibliography.md

notation.md: Rmarkdown/notation.Rmd
	cd Rmarkdown; Rscript -e "knitr::knit(input = 'notation.Rmd', output = '../markdown/notation.md')"

notation.pdf: notation.md
	pandoc -s --number-sections -o pdf/notation.pdf markdown/notation.md

introduction.pdf: bibtex/phd_thesis.bib markdown/introduction.md markdown/bibliography.md
	mkdir -p pdf
	pandoc -s --table-of-contents --number-sections --bibliography=bibtex/phd_thesis.bib -o pdf/introduction.pdf markdown/introduction.md markdown/bibliography.md

statmodel.pdf: bibtex/phd_thesis.bib markdown/a_statistical_model_of_methlyC-seq.md markdown/bibliography.md
	mkdir -p pdf
	pandoc -s --table-of-contents --number-sections --bibliography=bibtex/phd_thesis.bib -o pdf/a_statistical_model_of_methlyC-seq.pdf markdown/a_statistical_model_of_methlyC-seq.md markdown/bibliography.md markdown/bibliography.md

phd_thesis.pdf: bibtex/phd_thesis.bib markdown/preamble.md notation.md markdown/introduction.md markdown/bibliography.md
	mkdir -p pdf
	pandoc --table-of-contents --number-sections --bibliography=bibtex/phd_thesis.bib -o pdf/phd_thesis.pdf markdown/preamble.md markdown/notation.md markdown/introduction.md markdown/bibliography.md

phd_thesis.html: bibtex/phd_thesis.bib markdown/preamble.md notation.md markdown/introduction.md markdown/bibliography.md
	mkdir -p html
	pandoc -s --mathjax --table-of-contents --number-sections --bibliography=bibtex/phd_thesis.bib -o html/phd_thesis.html markdown/preamble.md markdown/notation.md markdown/introduction.md markdown/bibliography.md

all: phd_thesis.pdf phd_thesis.html