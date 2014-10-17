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

comethylation.md: latex/phd_thesis.bib Rmarkdown/comethylation.Rmd
	Rscript -e "rmarkdown::render('Rmarkdown/comethylation.Rmd', output_format = 'html_document', output_dir = '../markdown')"
	rm markdown/comethylation.html

comethylation.latex: comethylation.md latex/phd_thesis.bib Rmarkdown/comethylation.Rmd
	mkdir -p latex
	pandoc --chapters -o latex/comethylation.tex markdown/comethylation.md

comethylation.pdf: latex/phd_thesis.bib Rmarkdown/comethylation.Rmd markdown/bibliography.md
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/comethylation.Rmd', output_format = 'pdf_document', output_file = '../comethylation.pdf', output_dir = '../pdf')"

comethylation_review.md: latex/phd_thesis.bib Rmarkdown/comethylation_review.Rmd
	Rscript -e "rmarkdown::render('Rmarkdown/comethylation_review.Rmd', output_format = 'html_document', output_dir = '../markdown')"
	rm markdown/comethylation_review.html

comethylation_review.latex: comethylation_review.md latex/phd_thesis.bib Rmarkdown/comethylation_review.Rmd
	mkdir -p latex
	pandoc --chapters -o latex/comethylation_review.tex markdown/comethylation_review.md

comethylation_review.pdf: latex/phd_thesis.bib Rmarkdown/comethylation_review.Rmd markdown/bibliography.md
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/comethylation_review.Rmd', output_format = 'pdf_document', output_file = '../comethylation_review.pdf', output_dir = '../pdf')"

methsim.md: latex/phd_thesis.bib Rmarkdown/methsim.Rmd
	Rscript -e "rmarkdown::render('Rmarkdown/methsim.Rmd', output_format = 'html_document', output_dir = '../markdown')"
	rm markdown/methsim.html

methsim.latex: methsim.md latex/phd_thesis.bib Rmarkdown/methsim.Rmd
	mkdir -p latex
	pandoc --chapters -o latex/methsim.tex markdown/methsim.md

methsim.pdf: latex/phd_thesis.bib Rmarkdown/methsim.Rmd markdown/bibliography.md
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/methsim.Rmd', output_format = 'pdf_document', output_file = '../methsim.pdf', output_dir = '../pdf')"

phd_thesis.pdf: introduction.latex statmodel.latex wgbs_analysis.latex comethylation_review.latex comethylation.latex
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

phd_thesis.html: comethylation.md
	echo "Citations aren't yet supported!"
	mkdir html
	pandoc -s --mathjax --table-of-contents --number-sections --bibliography=latex/phd_thesis.bib -o html/phd_thesis.html markdown/preamble.md markdown/introduction.md markdown/a_statistical_model_of_methlyC-seq.md markdown/comethylation.md markdown/bibliography.md

all: phd_thesis.pdf phd_thesis.html
