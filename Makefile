# TODO: Standalone pdfs are missing references.

phd_thesis.pdf: introduction.tex biology_background.tex wgbs_bioinformatics_analysis.tex wgbs_statistical_framework.tex datasets.tex wgbs_downstream_analyses.tex comethylation_review.tex comethylation.tex methsim.tex concluding_remarks.tex appendix.tex
	mkdir -p pdf
	cd latex; pdflatex phd_thesis; \
	bibtex phd_thesis; \
	pdflatex phd_thesis; \
	pdflatex phd_thesis; \
	mv phd_thesis.pdf ../pdf/
# Annoyingly the printers can't seem to print the pdf! So this version replaces
# all pdf figures with high-quality pngs.
png_phd_thesis.pdf: introduction.tex biology_background.tex wgbs_bioinformatics_analysis.tex wgbs_statistical_framework.tex datasets.tex wgbs_downstream_analyses.tex comethylation_review.tex comethylation.tex methsim.tex concluding_remarks.tex appendix.tex
		mogrify -trim -density 300 -format png figures/*.pdf
		sed -i "" "s/pdf}/png}/g" latex/*.tex
		mkdir -p pdf
		cd latex; pdflatex phd_thesis; \
		bibtex phd_thesis; \
		pdflatex phd_thesis; \
		pdflatex phd_thesis; \
		mv phd_thesis.pdf ../pdf/peter-hickey-phd-thesis-png.pdf

introduction.tex: introduction.md
	mkdir -p latex
	pandoc -o latex/introduction.tex markdown/introduction.md

introduction.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/introduction.Rmd', 'markdown/introduction.md')"

introduction.pdf:
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/introduction.Rmd', output_format = 'pdf_document', output_file = '../introduction.pdf', output_dir = '../pdf')"

biology_background.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/biology_background.Rmd', 'markdown/biology_background.md')"

biology_background.tex: biology_background.md
	mkdir -p latex
	pandoc -o latex/biology_background.tex markdown/biology_background.md

biology_background.pdf:
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/biology_background.Rmd', output_format = 'pdf_document', output_file = '../biology_background.pdf', output_dir = '../pdf')"

wgbs_bioinformatics_analysis.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/wgbs_bioinformatics_analysis.Rmd', 'markdown/wgbs_bioinformatics_analysis.md')"

wgbs_bioinformatics_analysis.tex: wgbs_bioinformatics_analysis.md
	mkdir -p latex
	pandoc -o latex/wgbs_bioinformatics_analysis.tex markdown/wgbs_bioinformatics_analysis.md

wgbs_bioinformatics_analysis.pdf:
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/wgbs_bioinformatics_analysis.Rmd', output_format = 'pdf_document', output_file = '../wgbs_bioinformatics_analysis.pdf', output_dir = '../pdf')"

wgbs_statistical_framework.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/wgbs_statistical_framework.Rmd', 'markdown/wgbs_statistical_framework.md')"

wgbs_statistical_framework.tex: wgbs_statistical_framework.md
	mkdir -p latex
	pandoc -o latex/wgbs_statistical_framework.tex markdown/wgbs_statistical_framework.md

wgbs_statistical_framework.pdf:
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/wgbs_statistical_framework.Rmd', output_format = 'pdf_document', output_file = '../wgbs_statistical_framework.pdf', output_dir = '../pdf')"

wgbs_downstream_analyses.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/wgbs_downstream_analyses.Rmd', 'markdown/wgbs_downstream_analyses.md')"

wgbs_downstream_analyses.tex: wgbs_downstream_analyses.md
	mkdir -p latex
	pandoc -o latex/wgbs_downstream_analyses.tex markdown/wgbs_downstream_analyses.md

wgbs_downstream_analyses.pdf:
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/wgbs_downstream_analyses.Rmd', output_format = 'pdf_document', output_file = '../wgbs_downstream_analyses.pdf', output_dir = '../pdf')"

datasets.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/datasets.Rmd', 'markdown/datasets.md')"

datasets.tex: datasets.md
	mkdir -p latex
	pandoc -o latex/datasets.tex markdown/datasets.md

datasets.pdf:
	Rscript -e "rmarkdown::render('Rmarkdown/datasets.Rmd', output_format = 'pdf_document', output_file = '../datasets.pdf', output_dir = '../pdf')"

Avy_epialleles.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/Avy_epialleles.Rmd', 'markdown/Avy_epialleles.md')"

Avy_epialleles.tex:
	mkdir -p latex
	pandoc -o latex/Avy_epialleles.tex markdown/Avy_epialleles.md

Avy_epialleles.pdf:
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/Avy_epialleles.Rmd', output_format = 'pdf_document', output_file = '../Avy_epialleles.pdf', output_dir = '../pdf')"

comethylation_review.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/comethylation_review.Rmd', 'markdown/comethylation_review.md')"

comethylation_review.tex: comethylation_review.md
	mkdir -p latex
	pandoc -o latex/comethylation_review.tex markdown/comethylation_review.md

comethylation_review.pdf:
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/comethylation_review.Rmd', output_format = 'pdf_document', output_file = '../comethylation_review.pdf', output_dir = '../pdf')"

comethylation.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/comethylation.Rmd', 'markdown/comethylation.md')"

comethylation.tex: comethylation.md
	mkdir -p latex
	pandoc -o latex/comethylation.tex markdown/comethylation.md

comethylation.pdf:
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/comethylation.Rmd', output_format = 'pdf_document', output_file = '../comethylation.pdf', output_dir = '../pdf')"

methsim.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/methsim.Rmd', 'markdown/methsim.md')"

methsim.tex: methsim.md
	mkdir -p latex
	pandoc -o latex/methsim.tex markdown/methsim.md

methsim.pdf:
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/methsim.Rmd', output_format = 'pdf_document', output_file = '../methsim.pdf', output_dir = '../pdf')"

concluding_remarks.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/concluding_remarks.Rmd', 'markdown/concluding_remarks.md')"

concluding_remarks.tex: concluding_remarks.md
	mkdir -p latex
	pandoc -o latex/concluding_remarks.tex markdown/concluding_remarks.md

concluding_remarks.pdf:
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/concluding_remarks.Rmd', output_format = 'pdf_document', output_file = '../methsim.pdf', output_dir = '../pdf')"

appendix.md:
	mkdir -p markdown
	Rscript -e "knitr::knit('Rmarkdown/appendix.Rmd', 'markdown/appendix.md')"

appendix.tex: appendix.md
	mkdir -p latex
	pandoc -o latex/appendix.tex markdown/appendix.md

appendix.pdf:
	mkdir -p pdf
	Rscript -e "rmarkdown::render('Rmarkdown/appendix.Rmd', output_format = 'pdf_document', output_file = '../appendix.pdf', output_dir = '../pdf')"

paper_reviews.pdf: phd_thesis.bib paper_reviews.md markdown/bibliography.md
	mkdir -p pdf
	pandoc --natbib -s --table-of-contents --number-sections --bibliography=latex/phd_thesis.bib -o latex/paper_reviews_sc.tex markdown/paper_reviews.md markdown/bibliography.md
	pdflatex -output-directory latex latex/paper_reviews_sc
	bibtex latex/paper_reviews_sc
	pdflatex -output-directory latex latex/paper_reviews_sc
	pdflatex -output-directory latex latex/paper_reviews_sc
	mv latex/paper_reviews_sc.pdf pdf/paper_reviews.pdf
	rm latex/paper_reviews_sc*

phd_thesis.html: introduction.md biology_background.md wgbs_bioinformatics_analysis.md wgbs_statistical_framework.md datasets.md wgbs_downstream_analyses.md comethylation_review.md comethylation.md methsim.md concluding_remarks.md appendix.md
	echo "Citations aren't yet supported!"
	mkdir -p html
	pandoc -s --mathjax --table-of-contents --number-sections --bibliography=latex/phd_thesis.bib -o html/phd_thesis.html markdown/introduction.md markdown/biology_background.md markdown/wgbs_bioinformatics_analysis.md markdown/wgbs_statistical_framework.md markdown/datasets.md markdown/wgbs_downstream_analyses.md markdown/comethylation_review.md markdown/comethylation.md markdown/methsim.md markdown/appendix.md
all: phd_thesis.pdf phd_thesis.html
