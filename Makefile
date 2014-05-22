outline.pdf: latex/phd_thesis.bib markdown/outline.md markdown/bibliography.md
	mkdir -p pdf
	pandoc -s --number-sections --bibliography=latex/phd_thesis.bib -o pdf/outline.pdf markdown/outline.md markdown/bibliography.md

notation.pdf: Rmarkdown/notation.Rmd
	mkdir -p pdf
	Rscript -e "rmarkdown::render(input = 'Rmarkdown/notation.Rmd', output_file = '../pdf/notation.pdf')"

introduction.latex: markdown/introduction.md
	mkdir -p latex
	pandoc --chapters -o latex/introduction.tex markdown/introduction.md

introduction.pdf: latex/phd_thesis.bib markdown/introduction.md markdown/bibliography.md
	mkdir -p pdf
	pandoc --natbib -s --table-of-contents --number-sections --bibliography=latex/phd_thesis.bib -o latex/introduction_sc.tex markdown/introduction.md markdown/bibliography.md
	pdflatex -output-directory latex latex/introduction_sc
	bibtex latex/introduction_sc
	pdflatex -output-directory latex latex/introduction_sc
	pdflatex -output-directory latex latex/introduction_sc
	mv latex/introduction_sc.pdf pdf/introduction.pdf
	rm latex/introduction_sc*

statmodel.latex: latex/phd_thesis.bib markdown/a_statistical_model_of_methlyC-seq.md
	mkdir -p latex
	pandoc --chapters -o latex/a_statistical_model_of_methlyC-seq.tex markdown/a_statistical_model_of_methlyC-seq.md

statmodel.pdf: latex/phd_thesis.bib markdown/a_statistical_model_of_methlyC-seq.md markdown/bibliography.md
	mkdir -p pdf
	pandoc --natbib -s --table-of-contents --number-sections --bibliography=latex/phd_thesis.bib -o latex/a_statistical_model_of_methlyC-seq_sc.tex markdown/a_statistical_model_of_methlyC-seq.md markdown/bibliography.md
	pdflatex -output-directory latex latex/a_statistical_model_of_methlyC-seq_sc
	bibtex latex/a_statistical_model_of_methlyC-seq_sc
	pdflatex -output-directory latex latex/a_statistical_model_of_methlyC-seq_sc
	pdflatex -output-directory latex latex/a_statistical_model_of_methlyC-seq_sc
	mv latex/a_statistical_model_of_methlyC-seq_sc.pdf pdf/a_statistical_model_of_methlyC-seq.pdf
	rm latex/a_statistical_model_of_methlyC-seq_sc*

wgbs_analysis.latex: latex/phd_thesis.bib markdown/wgbs_analysis.md
	mkdir -p latex
	pandoc --chapters -o latex/wgbs_analysis.tex markdown/wgbs_analysis.md

wgbs_analysis.pdf: latex/phd_thesis.bib markdown/wgbs_analysis.md markdown/bibliography.md
	mkdir -p pdf
	pandoc --natbib -s --table-of-contents --number-sections --bibliography=latex/phd_thesis.bib -o latex/wgbs_analysis_sc.tex markdown/wgbs_analysis.md markdown/bibliography.md
	pdflatex -output-directory latex latex/wgbs_analysis_sc
	bibtex latex/wgbs_analysis_sc
	pdflatex -output-directory latex latex/wgbs_analysis_sc
	pdflatex -output-directory latex latex/wgbs_analysis_sc
	mv latex/wgbs_analysis_sc.pdf pdf/wgbs_analysis.pdf
	rm latex/wgbs_analysis_sc*

phd_thesis.pdf: introduction.latex statmodel.latex wgbs_analysis.latex
	mkdir -p pdf
	cd latex; pdflatex phd_thesis; \
	bibtex phd_thesis; \
	pdflatex phd_thesis; \
	pdflatex phd_thesis; \
	mv phd_thesis.pdf ../pdf/

paper_reviews.pdf: latex/phd_thesis.bib markdown/paper_reviews.md markdown/bibliography.md
	mkdir -p pdf
	pandoc --natbib -s --table-of-contents --number-sections --bibliography=latex/phd_thesis.bib -o latex/paper_reviews_sc.tex markdown/paper_reviews.md markdown/bibliography.md
	pdflatex -output-directory latex latex/paper_reviews_sc
	bibtex latex/paper_reviews_sc
	pdflatex -output-directory latex latex/paper_reviews_sc
	pdflatex -output-directory latex latex/paper_reviews_sc
	mv latex/paper_reviews_sc.pdf pdf/paper_reviews.pdf
	rm latex/paper_reviews_sc*

phd_thesis.html:
	echo "Citations aren't yet supported!"
	mkdir html
	pandoc -s --mathjax --table-of-contents --number-sections --bibliography=latex/phd_thesis.bib -o html/phd_thesis.html markdown/preamble.md markdown/introduction.md markdown/a_statistical_model_of_methlyC-seq.md markdown/bibliography.md

all: phd_thesis.pdf phd_thesis.html
