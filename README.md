This repository holds the source code and data for a manuscript in progress.

Usage tips:

+ [data](/data/) contains the structures (SMILES) and hit status of compounds screened, enrichment data for macrocycles, and melting temperature data for thermal shift screens.
+ [src](/src/) contains the R source code to analyze data and create figures.
    + This script creates Table 1 and Figures 1A-D, 3A-D, 4A, 5A, S1, and S4B. The remainder of the display items were laid out in Microsoft Word or Adobe Illustrator.
    + Dependencies are [`sqldf`](https://cran.r-project.org/web/packages/sqldf/) and, optionally, [`rcdk`](https://cran.r-project.org/web/packages/rcdk/index.html) (by default the `recompute_fragment_properties` flag, line 61, is set to false, so that the time-consuming step of recalcuating the fragment properties is skipped).
    + To run it and re-generate the figures, just run `Rscript src/binder_figures.R`. It should complete in about 10 seconds.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

