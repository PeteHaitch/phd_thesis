
# Short term

- [x] Avoid first person in thesis? Ask Terry for opinion.
  - Modify if I grimace while reading the sentence.
- [x] Give each chapter a summary.
  - Checked each chapter
- [x] Put table captions at the top, figure captions at the bottom.
  - Check in proofreading
- [x] Add genome-wide line to legend of MH plots
  - __DONE__
- [x] Address Saskia's comments
  - __DONE__
- [x] Address Goknur's comments
  - __DONE__
- [x] Address Terry's comments
  - __DONE__
- [x] $\mathbf{\psi}$ doesn't print how I'd hoped.
  - __FIXED__ Use $\bm{\psi}$
- [x] Rename `methsim` plots without the `_1` and `_2` suffixes
  - __DONE__
- [x] Write abstract
  - __DONE__
- [x] Complete `methsim` chapter.
  - __DONE__
- [x] Write concluding remarks
  - __DONE__
- [x] Email communications about printing
  - __DONE__
- [x] Add blank page between title and abstract
  - __DONE__
- [x] Give every figure and table a brief caption.
  - __DONE__
- [x] No full stop in title
  __DONE__
- [x] Does the WGBS bioinformatics analysis chapter need a summary?
  - __DONE: added__
- [x] Re-visit summary of WGBS statistical framework chapter
  - __DONE__
- [x] What is the general restriction on $c$ when $p_{h, i} \neq p_{h', i'}$?
  __DONE: All of $\Pi$ must be $\in [0, 1]$ and alll rows of $\Pi$ must sum to 1.__
- [x] Should gene names and cytosine modifications be _italicised_.
   - Yes to gene names, no to cytosine modifications.
   - __DONE__
- [x] "The conditional estimator, works better than the unconditional estimator when the sample size, n, is small [Agresti 2007]." re conditional odds ratio estimators
  - What does "better" mean.
  - __DONE: Cited page in Agresti's intro to CDA textbook (2007).__
- [x] Figure out why E18BUF is missing chrM for within-fragment $NIL \geq 0$ plots.
  - __DONE: I've noted the E18BUF omission in the caption. E18BUF has no chrM data in the 2ac-tuples `MethPat` object because methtuple was never run on that sample-chromosome-2ac combination (nor was chrY). This was due to some sort of omission in the `run_methtuple.sh` script__
- [x] Figure out what's going on with some of the EPISCOPE and Ziller samples that have a high autosome (chr21, I think).
  - __DONE: Yes, it's chr21. But I have no idea what's causing this. Truly weird. I've commented on this in my thesis__.
- [x] \texttt{--all-combinations} prints as '-all-combinations' and other examples
  - __DONE: Replaced with \texttt{{-}{-}all-combinations}__
- [x] Fix formatting of $2 \times 2$ contingency tables to be consistent.
  - __DONE__
- [x] Find submission documents
  - __DONE__
- [x] Tables in datasets chapter run outside margins and are generally a bit of a mess.
  - __DONE: Fixed these tables__
- [x] Check typesetting of C++
  - __DONE: Used a macro__
- [x] Check Sherman `-CG` and `-CH` parameters (e.g., `--CG` or `-CG`)
  - __DONE__
- [x] "dataset" vs. "data set".
  - __DONE: Using dataset__
- [x] Include timings from Lacey et al. in `methsim` chapter.
  - __DONE__
- [x] Check the size of the simulated data as a `MethPat` object.
  - __DONE: updated__
- [x] Plots of beta-correlations separated by strand in a clearer way.
  - __DONE: Added loess fit__
- [x] Include `coMET` paper in Chapter 6.
  - __DONE__
- [x] Fix `BiocParallel` citation.
  - __DONE__
- [x] Look at some scatterplots of beta-value pairs from simulated data. Include in `methsim` chapter if revealing.
  - __DONE: They are in supplementary code. Nothing revealing enough to justify inclusion in main text; basically, the simulated data lie along the diagonal far more than the real data, hence they have higher correlations.__
- [x] Move estimation of RRBS CpGs from the appendix to online.
  - __DONE__
- [x] Add a representative plot of the within-fragment data for a single sample
  - __DONE: Decided against this; there are enough plots already__
- [x] Increase axis labels in simulation study figures.
  - __DONE__
- [x] Update title of methtuple memory usage plot to refer to 40 samples rather than 49.
  - __DONE: Don't do this. There were 49 runs, so use that number__
- [x] Check that all Ziller plots are titled "Ziller" and not "Ziller merged"
  - __DONE__
- [x] Add $\beta$ to the y-axis of the $\beta$-by-PM boxplots.
  - __DONE__- [ ] Fix MethylationTuples::mantelhaen(). Works if copy-pasted into session but not as part of MethylationTuples package. Then update Software Details in Appendix.
- [x] Explain why Seisenberger only has 2 partitioned methylomes
  - __DONE__
- [x] Address Dineika's suggestions and add Dineika to acknowledgements
  - __DONE__
- [x] Check fonts on `methtuple` performance plots are legible when printed.
  - __DONE: Made figure larger but didn't alter x-axis labels (can't, not enough space).__
- [x] Terry's point re convergence argument (p196)
  - __DONE: Re-worded to make this clearer__
- [x] Add `makeCGI` to list of software.
  - __DONE__
- [x] Read unimelb docs on submission
  - __DONE__
- [x] Check Cokus equations about joint probs
  - __DONE: Correct in text__
- [x] Check Breslow:1981ta for definitions of odds ratio estimates in 3-way tables
  - __DONE__
- [x] Check latex logs for warnings/errors
  - __DONE__


## Thursday tasks


## Friday tasks

- [ ] Add analysis folder to github
  - First upload `.Rmd` analysis scripts.
  - Upload <sample>.CG.<m>.tsv.gz Lister, Seisenberger, Ziller samples (m = 1, 2, 2ac) to figshare. All 1-tuples are under Figshare's 250 MB limit but the same is not true of 2-tuples or 2ac-tuples.
- [ ] Compress PDF

# Long term

- MethylomeParam and simulate,MethylomeParam-method are going to need _another_ re-write :(
  - Slots will be something like MethLevelDT, ComethDT (using Mantel-Haenszel estimates), MethLevelCorDT and all will have to be with respect to the levels of the PartitionedMethylome.
- Once switch from "2" to "1" classes/methods made in methsim, fix other TODOs
- Suggest BPPARAM argument for bsapply(). Suggest bpreplicate() added to BiocParallel.
- Could simulate $B$ from Beta distribution rather than Normal, i.e. beta-binomial model! This is what Lacey et al. or `WGBSsuite` do.
- Figure out why an assay in a SummarizedExperiment can be a DataFrame but cannot be inserted using `assay<-`.
- Read Hadley's recommendations on capitalising parameter descriptions and implement them for all packages.
- Finish ASD2 report.
- Update excludes for coveralls packages
