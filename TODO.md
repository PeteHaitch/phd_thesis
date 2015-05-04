
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
- [ ] Figure out what's going on with some of the EPISCOPE and Ziller samples that have a high autosome (chr21, I think); figure out why E18BUF is missing chrM in $NIL \geq 0$ plots.
- [ ] Should gene names and cytosine modifications be _italicised_.
- [ ] Plots of beta-correlations separated by strand in a clearer way.
  - Use loess to summarise trend
  - In general, plots of beta-correlations are noisy. Just show trend?
  - Plot the difference of correlations between strands?
- [ ] "The conditional estimator, works better than the unconditional estimator when the sample size, n, is small [Agresti 2007]." re conditional odds ratio estimators
  - What does "better" mean.
- [ ] Add a representative plot of the within-fragment data for a single sample
- [ ] Fix formatting of $2 \times 2$ contingency tables to be consistent.
- [ ] Increase axis labels in simulation study figures?
- [ ] Find submission documents
- [ ] Add analysis folder to github
  - First upload `.Rmd` analysis scripts.
  - Upload <sample>.CG.<m>.tsv.gz Lister, Seisenberger, Ziller samples (m = 1, 2, 2ac) to figshare. All 1-tuples are under Figshare's 250 MB limit but the same is not true of 2-tuples or 2ac-tuples.
- [ ] Fix MethylationTuples::mantelhaen(). Works if copy-pasted into session but not as part of MethylationTuples package.
- [ ] Look at some scatterplots of beta-value pairs from simulated data
- "dataset" vs. "data set". Include in `methsim` chapter if revealing.
- [ ] \texttt{--all-combinations} prints as '-all-combinations' and other examples
  - Works when typed in markdown as `--all-combinations`, but this isn't always possible.
- [ ] Update title of methtuple memory usage plot to refer to 40 samples rather than 49.
- [ ] Tables in datasets chapter run outside margins and are generally a bit of a mess.
- [ ] Check that all Ziller plots are titled "Ziller" and not "Ziller merged"
- [ ] Check typesetting of C++
- [ ] Include `coMET` paper in Chapter 6.
- [ ] Clarify the definition of $\psi_{U}$ for $2 \times 2 \times K$; is it the same as the "collapsed" estimator?
- [ ] Check Sherman `-CG` and `-CH` parameters (e.g., `--CG` or `-CG`)
- [ ] Check the size of the simulated data as a `MethPat` object.
- [ ] Add $\beta$ to the y-axis of the $\beta$-by-PM boxplots.
- [ ] Fix `BiocParallel` citation.
- [ ] Include timings from Lacey et al. in `methsim` chapter.
- [ ] Be more specific with what was discovered in the analyses of co-methylation.
- [ ] Move estimation of RRBS CpGs from the appendix to online.
- [ ] Add `makeCGI` to list of software.
- [ ] Multiple $\LaTeX$ typesetting bugs.

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
