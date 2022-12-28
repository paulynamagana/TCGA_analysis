{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "24c5df42",
   "metadata": {
    "_execution_state": "idle",
    "_uuid": "051d70d956493feee0c6d64651c6a088724dca2a",
    "execution": {
     "iopub.execute_input": "2022-12-28T12:53:42.796380Z",
     "iopub.status.busy": "2022-12-28T12:53:42.794398Z",
     "iopub.status.idle": "2022-12-28T13:53:03.584492Z",
     "shell.execute_reply": "2022-12-28T13:53:03.582063Z"
    },
    "papermill": {
     "duration": 3560.803922,
     "end_time": "2022-12-28T13:53:03.588288",
     "exception": false,
     "start_time": "2022-12-28T12:53:42.784366",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Bioconductor version '3.12' is out-of-date; the current release version '3.16'\n",
      "  is available with R version '4.2'; see https://bioconductor.org/install\n",
      "\n",
      "'getOption(\"repos\")' replaces Bioconductor standard repositories, see\n",
      "'?repositories' for details\n",
      "\n",
      "replacement repositories:\n",
      "    CRAN: http://cran.rstudio.com/\n",
      "\n",
      "\n",
      "Bioconductor version 3.12 (BiocManager 1.30.18), R 4.0.5 (2021-03-31)\n",
      "\n",
      "Installing package(s) 'limma', 'minfi', 'RColorBrewer', 'missMethyl',\n",
      "  'minfiData', 'Gviz', 'sesame', 'DMRcate', 'DMRcatedata', 'stringr', 'readr',\n",
      "  'tidyverse', 'edgeR', 'ggpubr',\n",
      "  'IlluminaHumanMethylation450kanno.ilmn12.hg19',\n",
      "  'IlluminaHumanMethylation450kmanifest', 'stringr', 'TCGAbiolinks'\n",
      "\n",
      "also installing the dependencies ‘TxDb.Hsapiens.UCSC.hg19.knownGene’, ‘interactiveDisplayBase’, ‘MatrixGenerics’, ‘GenomeInfoDbData’, ‘multtest’, ‘scrime’, ‘sparseMatrixStats’, ‘annotate’, ‘FDb.InfiniumMethylation.hg19’, ‘zlibbioc’, ‘BiocFileCache’, ‘AnnotationFilter’, ‘ProtGenerics’, ‘VariantAnnotation’, ‘Rhtslib’, ‘AnnotationHub’, ‘beachmat’, ‘vctrs’, ‘GenomicRanges’, ‘SummarizedExperiment’, ‘Biostrings’, ‘bumphunter’, ‘S4Vectors’, ‘GenomeInfoDb’, ‘Biobase’, ‘IRanges’, ‘nor1mix’, ‘siggenes’, ‘preprocessCore’, ‘illuminaio’, ‘DelayedMatrixStats’, ‘genefilter’, ‘GEOquery’, ‘DelayedArray’, ‘HDF5Array’, ‘BiocParallel’, ‘IlluminaHumanMethylationEPICanno.ilm10b4.hg19’, ‘AnnotationDbi’, ‘GO.db’, ‘IlluminaHumanMethylationEPICmanifest’, ‘methylumi’, ‘org.Hs.eg.db’, ‘ruv’, ‘XVector’, ‘rtracklayer’, ‘biomaRt’, ‘GenomicFeatures’, ‘ensembldb’, ‘BSgenome’, ‘biovizBase’, ‘Rsamtools’, ‘GenomicAlignments’, ‘sesameData’, ‘wheatmap’, ‘DNAcopy’, ‘ExperimentHub’, ‘bsseq’, ‘DSS’, ‘ggplot2’, ‘rstatix’, ‘TCGAbiolinksGUI.data’\n",
      "\n",
      "\n",
      "Warning message in install.packages(...):\n",
      "“installation of package ‘vctrs’ had non-zero exit status”\n",
      "Warning message in install.packages(...):\n",
      "“installation of package ‘ggplot2’ had non-zero exit status”\n",
      "Warning message in install.packages(...):\n",
      "“installation of package ‘ggpubr’ had non-zero exit status”\n",
      "Old packages: 'ade4', 'afex', 'amap', 'Amelia', 'arrow', 'arules',\n",
      "  'AsioHeaders', 'assertr', 'bayesm', 'bayesplot', 'bbotk', 'BiasedUrn',\n",
      "  'bigmemory.sri', 'bigrquery', 'BiocManager', 'BioStatR', 'bipartite', 'bit',\n",
      "  'bonsai', 'bookdown', 'Boom', 'boot', 'Boruta', 'Brobdingnag', 'broom',\n",
      "  'broom.helpers', 'bslib', 'bsts', 'butcher', 'C50', 'callr', 'car',\n",
      "  'CausalImpact', 'cba', 'ChemoSpec', 'ChemoSpecUtils', 'cli', 'clue',\n",
      "  'collections', 'colourpicker', 'compareGroups', 'copula', 'CORElearn',\n",
      "  'credentials', 'cropcircles', 'crypto2', 'Cubist', 'data.table', 'datamods',\n",
      "  'datawizard', 'date', 'dbscan', 'deaR', 'DescTools', 'deSolve', 'devEMF',\n",
      "  'devtools', 'dials', 'digest', 'distr', 'distrEx', 'dlm', 'dlookr', 'doBy',\n",
      "  'doRNG', 'DoubleML', 'DT', 'duckdb', 'e1071', 'effectsize', 'eha', 'energy',\n",
      "  'epiR', 'ergm', 'evaluate', 'exactextractr', 'explore', 'fasterize',\n",
      "  'fBasics', 'fGarch', 'fixest', 'flexsurv', 'flextable', 'fmsb',\n",
      "  'fontawesome', 'forecast', 'foreign', 'formatR', 'fracdiff', 'fst', 'future',\n",
      "  'future.apply', 'gam', 'gamlss', 'gdalUtilities', 'gee', 'geodist',\n",
      "  'geohashTools', 'geosphere', 'geostan', 'gert', 'ggbeeswarm', 'ggborderline',\n",
      "  'ggdag', 'ggeffects', 'ggfortify', 'ggfun', 'gghalves', 'ggiraph', 'ggmap',\n",
      "  'ggpie', 'ggplot2', 'ggpmisc', 'ggpp', 'ggprism', 'ggpubr', 'ggrepel',\n",
      "  'ggshadow', 'ggside', 'ggspatial', 'ggstance', 'ggstatsplot', 'gh',\n",
      "  'gitcreds', 'gld', 'glmmML', 'glmnet', 'globals', 'gower', 'GPArotation',\n",
      "  'graphlayouts', 'greybox', 'grf', 'groupdata2', 'gstat', 'gt', 'gtools',\n",
      "  'gtsummary', 'hablar', 'HDInterval', 'heplots', 'highr', 'Hmisc',\n",
      "  'htmltools', 'htmlwidgets', 'httpuv', 'imager', 'infer', 'insight',\n",
      "  'installr', 'IRkernel', 'ISLR2', 'isoband', 'jpeg', 'jsonify', 'jsonlite',\n",
      "  'jtools', 'keras', 'kit', 'knitr', 'ks', 'L1pack', 'latex2exp', 'lava',\n",
      "  'lemon', 'lessR', 'lhs', 'LiblineaR', 'lidR', 'lightgbm', 'lintr', 'listenv',\n",
      "  'littler', 'lme4', 'locpol', 'LogicReg', 'logspline', 'lpSolveAPI',\n",
      "  'lubridate', 'lwgeom', 'magic', 'mapdata', 'mapproj', 'maps', 'maptools',\n",
      "  'markdown', 'Matching', 'matlib', 'Matrix', 'matrixStats', 'mc2d', 'mclust',\n",
      "  'meta', 'mgcv', 'mice', 'minqa', 'mixtools', 'mlr3', 'mlr3cluster',\n",
      "  'mlr3fselect', 'mlr3hyperband', 'mlr3learners', 'mlr3tuning',\n",
      "  'mlr3tuningspaces', 'mlr3verse', 'modelr', 'modelsummary', 'modeltime',\n",
      "  'moderndive', 'msm', 'ncdf4', 'NCmisc', 'ndjson', 'ndtv', 'nflfastR',\n",
      "  'ngram', 'nlme', 'np', 'officer', 'openair', 'OpenImageR', 'OpenMx',\n",
      "  'openssl', 'openxlsx', 'ordinal', 'ordinalForest', 'overlapping',\n",
      "  'packcircles', 'padr', 'pagedown', 'paletteer', 'paradox', 'parallelly',\n",
      "  'parameters', 'ParBayesianOptimization', 'parsnip', 'pbapply', 'pbdZMQ',\n",
      "  'pcaPP', 'PDtoolkit', 'performance', 'phenofit', 'philentropy', 'phyclust',\n",
      "  'pkgbuild', 'pkgload', 'plgp', 'plotly', 'plsVarSel', 'plyr', 'png',\n",
      "  'polspline', 'polyclip', 'prevR', 'processx', 'progressr', 'proj4',\n",
      "  'protolite', 'ps', 'purrr', 'qgraph', 'qtl', 'qualtRics', 'quanteda',\n",
      "  'R.utils', 'ragg', 'randtoolbox', 'raster', 'rasterVis', 'rayimage',\n",
      "  'rbibutils', 'rcmdcheck', 'RcppAnnoy', 'RcppArmadillo', 'RcppDE',\n",
      "  'RcppEigen', 'reactable', 'rearrr', 'rebus.datetimes', 'recipes', 'reproj',\n",
      "  'rgdal', 'rgeos', 'riskRegression', 'rlas', 'RMariaDB', 'rmarkdown',\n",
      "  'rmutil', 'RMySQL', 'rngWELL', 'RNifti', 'robustlmm', 'RODBC', 'roxygen2',\n",
      "  'rpart', 'rsample', 'rsconnect', 'RSocrata', 'RSQLite', 'rstpm2', 'rsvg',\n",
      "  'rTLS', 'rugarch', 'rversions', 'rvg', 's2', 's20x', 'sasLM', 'sass', 'see',\n",
      "  'segmented', 'seqinr', 'seriation', 'servr', 'sessioninfo', 'Seurat',\n",
      "  'SeuratObject', 'sf', 'sfsmisc', 'shiny', 'shinyWidgets', 'sits', 'skimr',\n",
      "  'slider', 'SmartEDA', 'sp', 'sparklyr', 'spatstat', 'spatstat.data',\n",
      "  'spatstat.geom', 'spatstat.linnet', 'spatstat.random', 'spatstat.sparse',\n",
      "  'spatstat.utils', 'spData', 'spotifyr', 'stacks', 'stars', 'stats19',\n",
      "  'statsExpressions', 'stringdist', 'styler', 'survivalROC', 'sys', 'tables',\n",
      "  'targets', 'taylor', 'tensorflow', 'terra', 'testthat', 'text2vec',\n",
      "  'textrecipes', 'tgp', 'tibbletime', 'tidycensus', 'tidyquant', 'tidytable',\n",
      "  'tidytext', 'tigris', 'timeDate', 'timereg', 'timetk', 'tinytex', 'tm',\n",
      "  'tokenizers', 'topicmodels', 'torch', 'triangle', 'tuneR', 'tvthemes',\n",
      "  'udpipe', 'units', 'usmap', 'V8', 'vctrs', 'vioplot', 'vtable', 'vtree',\n",
      "  'warbleR', 'waterfalls', 'whisker', 'wk', 'workflows', 'worldfootballR',\n",
      "  'writexl', 'xfun', 'XML', 'yaml', 'yfR', 'yulab.utils', 'zip', 'class',\n",
      "  'cluster', 'KernSmooth', 'lattice', 'MASS', 'nnet', 'spatial', 'survival'\n",
      "\n",
      "Loading required package: minfi\n",
      "\n",
      "Loading required package: BiocGenerics\n",
      "\n",
      "Loading required package: parallel\n",
      "\n",
      "\n",
      "Attaching package: ‘BiocGenerics’\n",
      "\n",
      "\n",
      "The following objects are masked from ‘package:parallel’:\n",
      "\n",
      "    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,\n",
      "    clusterExport, clusterMap, parApply, parCapply, parLapply,\n",
      "    parLapplyLB, parRapply, parSapply, parSapplyLB\n",
      "\n",
      "\n",
      "The following objects are masked from ‘package:stats’:\n",
      "\n",
      "    IQR, mad, sd, var, xtabs\n",
      "\n",
      "\n",
      "The following objects are masked from ‘package:base’:\n",
      "\n",
      "    anyDuplicated, append, as.data.frame, basename, cbind, colnames,\n",
      "    dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,\n",
      "    grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,\n",
      "    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,\n",
      "    rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,\n",
      "    union, unique, unsplit, which.max, which.min\n",
      "\n",
      "\n",
      "Loading required package: GenomicRanges\n",
      "\n",
      "Loading required package: stats4\n",
      "\n",
      "Loading required package: S4Vectors\n",
      "\n",
      "\n",
      "Attaching package: ‘S4Vectors’\n",
      "\n",
      "\n",
      "The following object is masked from ‘package:base’:\n",
      "\n",
      "    expand.grid\n",
      "\n",
      "\n",
      "Loading required package: IRanges\n",
      "\n",
      "Loading required package: GenomeInfoDb\n",
      "\n",
      "Loading required package: SummarizedExperiment\n",
      "\n",
      "Loading required package: MatrixGenerics\n",
      "\n",
      "Loading required package: matrixStats\n",
      "\n",
      "\n",
      "Attaching package: ‘MatrixGenerics’\n",
      "\n",
      "\n",
      "The following objects are masked from ‘package:matrixStats’:\n",
      "\n",
      "    colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,\n",
      "    colCounts, colCummaxs, colCummins, colCumprods, colCumsums,\n",
      "    colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,\n",
      "    colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,\n",
      "    colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,\n",
      "    colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,\n",
      "    colWeightedMeans, colWeightedMedians, colWeightedSds,\n",
      "    colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,\n",
      "    rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,\n",
      "    rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,\n",
      "    rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,\n",
      "    rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,\n",
      "    rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,\n",
      "    rowWeightedMads, rowWeightedMeans, rowWeightedMedians,\n",
      "    rowWeightedSds, rowWeightedVars\n",
      "\n",
      "\n",
      "Loading required package: Biobase\n",
      "\n",
      "Welcome to Bioconductor\n",
      "\n",
      "    Vignettes contain introductory material; view with\n",
      "    'browseVignettes()'. To cite Bioconductor, see\n",
      "    'citation(\"Biobase\")', and for packages 'citation(\"pkgname\")'.\n",
      "\n",
      "\n",
      "\n",
      "Attaching package: ‘Biobase’\n",
      "\n",
      "\n",
      "The following object is masked from ‘package:MatrixGenerics’:\n",
      "\n",
      "    rowMedians\n",
      "\n",
      "\n",
      "The following objects are masked from ‘package:matrixStats’:\n",
      "\n",
      "    anyMissing, rowMedians\n",
      "\n",
      "\n",
      "The following object is masked from ‘package:httr’:\n",
      "\n",
      "    content\n",
      "\n",
      "\n",
      "Loading required package: Biostrings\n",
      "\n",
      "Loading required package: XVector\n",
      "\n",
      "\n",
      "Attaching package: ‘Biostrings’\n",
      "\n",
      "\n",
      "The following object is masked from ‘package:base’:\n",
      "\n",
      "    strsplit\n",
      "\n",
      "\n",
      "Loading required package: bumphunter\n",
      "\n",
      "Loading required package: foreach\n",
      "\n",
      "Loading required package: iterators\n",
      "\n",
      "Loading required package: locfit\n",
      "\n",
      "locfit 1.5-9.4 \t 2020-03-24\n",
      "\n",
      "Setting options('download.file.method.GEOquery'='auto')\n",
      "\n",
      "Setting options('GEOquery.inmemory.gpl'=FALSE)\n",
      "\n",
      "\n",
      "Attaching package: ‘minfi’\n",
      "\n",
      "\n",
      "The following object is masked from ‘package:TCGAbiolinks’:\n",
      "\n",
      "    getManifest\n",
      "\n",
      "\n",
      "\n",
      "Attaching package: ‘limma’\n",
      "\n",
      "\n",
      "The following object is masked from ‘package:BiocGenerics’:\n",
      "\n",
      "    plotMA\n",
      "\n",
      "\n",
      "Loading required package: IlluminaHumanMethylationEPICanno.ilm10b4.hg19\n",
      "\n",
      "\n",
      "Attaching package: ‘IlluminaHumanMethylationEPICanno.ilm10b4.hg19’\n",
      "\n",
      "\n",
      "The following objects are masked from ‘package:IlluminaHumanMethylation450kanno.ilmn12.hg19’:\n",
      "\n",
      "    Islands.UCSC, Locations, Manifest, Other, SNPs.132CommonSingle,\n",
      "    SNPs.135CommonSingle, SNPs.137CommonSingle, SNPs.138CommonSingle,\n",
      "    SNPs.141CommonSingle, SNPs.142CommonSingle, SNPs.144CommonSingle,\n",
      "    SNPs.146CommonSingle, SNPs.147CommonSingle, SNPs.Illumina\n",
      "\n",
      "\n",
      "\n",
      "\n",
      "Warning message:\n",
      "\"replacing previous import 'minfi::getMeth' by 'bsseq::getMeth' when loading 'DMRcate'\"\n",
      "Loading required package: grid\n",
      "\n",
      "── \u001b[1mAttaching packages\u001b[22m ─────────────────────────────────────── tidyverse 1.3.2 ──\n",
      "\u001b[32m✔\u001b[39m \u001b[34mtibble \u001b[39m 3.1.8      \u001b[32m✔\u001b[39m \u001b[34mdplyr  \u001b[39m 1.0.10\n",
      "\u001b[32m✔\u001b[39m \u001b[34mtidyr  \u001b[39m 1.2.1      \u001b[32m✔\u001b[39m \u001b[34mforcats\u001b[39m 0.5.2 \n",
      "\u001b[32m✔\u001b[39m \u001b[34mpurrr  \u001b[39m 0.3.5      \n",
      "── \u001b[1mConflicts\u001b[22m ────────────────────────────────────────── tidyverse_conflicts() ──\n",
      "\u001b[31m✖\u001b[39m \u001b[34mpurrr\u001b[39m::\u001b[32maccumulate()\u001b[39m masks \u001b[34mforeach\u001b[39m::accumulate()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mcollapse()\u001b[39m   masks \u001b[34mBiostrings\u001b[39m::collapse(), \u001b[34mIRanges\u001b[39m::collapse()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mcombine()\u001b[39m    masks \u001b[34mminfi\u001b[39m::combine(), \u001b[34mBiobase\u001b[39m::combine(), \u001b[34mBiocGenerics\u001b[39m::combine()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mpurrr\u001b[39m::\u001b[32mcompact()\u001b[39m    masks \u001b[34mXVector\u001b[39m::compact()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mBiobase\u001b[39m::\u001b[32mcontent()\u001b[39m  masks \u001b[34mhttr\u001b[39m::content()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mcount()\u001b[39m      masks \u001b[34mmatrixStats\u001b[39m::count()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mdesc()\u001b[39m       masks \u001b[34mIRanges\u001b[39m::desc()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mtidyr\u001b[39m::\u001b[32mexpand()\u001b[39m     masks \u001b[34mS4Vectors\u001b[39m::expand()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mfilter()\u001b[39m     masks \u001b[34mstats\u001b[39m::filter()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mfirst()\u001b[39m      masks \u001b[34mS4Vectors\u001b[39m::first()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mlag()\u001b[39m        masks \u001b[34mstats\u001b[39m::lag()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mpurrr\u001b[39m::\u001b[32mnone()\u001b[39m       masks \u001b[34mlocfit\u001b[39m::none()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mggplot2\u001b[39m::\u001b[32mPosition()\u001b[39m masks \u001b[34mBiocGenerics\u001b[39m::Position(), \u001b[34mbase\u001b[39m::Position()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mpurrr\u001b[39m::\u001b[32mreduce()\u001b[39m     masks \u001b[34mGenomicRanges\u001b[39m::reduce(), \u001b[34mIRanges\u001b[39m::reduce()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mrename()\u001b[39m     masks \u001b[34mS4Vectors\u001b[39m::rename()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mslice()\u001b[39m      masks \u001b[34mXVector\u001b[39m::slice(), \u001b[34mIRanges\u001b[39m::slice()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mpurrr\u001b[39m::\u001b[32mwhen()\u001b[39m       masks \u001b[34mforeach\u001b[39m::when()\n"
     ]
    }
   ],
   "source": [
    "PACKAGES <- c(\"limma\", \"minfi\", \"RColorBrewer\", \"missMethyl\",\n",
    "              \"minfiData\", \"Gviz\", \"sesame\", \"DMRcate\", \"DMRcatedata\", \"stringr\",\n",
    "              \"readr\", \"tidyverse\", \"edgeR\", \"ggpubr\",\n",
    "              \"IlluminaHumanMethylation450kanno.ilmn12.hg19\",\n",
    "              \"IlluminaHumanMethylation450kmanifest\", \"stringr\", \"TCGAbiolinks\")\n",
    "\n",
    "if (!require(\"BiocManager\", quietly = TRUE))\n",
    "    install.packages(\"BiocManager\")\n",
    "BiocManager::install(PACKAGES, force=TRUE)\n",
    "\n",
    "\n",
    "### Analyse data\n",
    "\n",
    "library(\"TCGAbiolinks\")\n",
    "library(\"IlluminaHumanMethylation450kanno.ilmn12.hg19\")\n",
    "library(\"IlluminaHumanMethylation450kmanifest\")\n",
    "library(\"minfi\")\n",
    "library(\"limma\")\n",
    "library(\"missMethyl\") # Can take a short time...\n",
    "library(\"DMRcate\")\n",
    "library(\"Gviz\")\n",
    "library(\"ggplot2\")\n",
    "library(\"RColorBrewer\")\n",
    "library(\"edgeR\")\n",
    "library(\"stringr\")\n",
    "library(\"readr\")\n",
    "library(\"tidyverse\")\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "0bef706d",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:53:04.048132Z",
     "iopub.status.busy": "2022-12-28T13:53:04.010088Z",
     "iopub.status.idle": "2022-12-28T13:53:04.060792Z",
     "shell.execute_reply": "2022-12-28T13:53:04.058767Z"
    },
    "papermill": {
     "duration": 0.464402,
     "end_time": "2022-12-28T13:53:04.063463",
     "exception": false,
     "start_time": "2022-12-28T13:53:03.599061",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "### NOT RUNNING CODE AS DATA HAS BEEN DOWNLOADED"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "5c32db25",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:53:04.089153Z",
     "iopub.status.busy": "2022-12-28T13:53:04.087213Z",
     "iopub.status.idle": "2022-12-28T13:53:04.101983Z",
     "shell.execute_reply": "2022-12-28T13:53:04.099621Z"
    },
    "papermill": {
     "duration": 0.031574,
     "end_time": "2022-12-28T13:53:04.104571",
     "exception": false,
     "start_time": "2022-12-28T13:53:04.072997",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "#query_met <- GDCquery(project= \"TCGA-LUSC\", \n",
    "                           #data.category = \"DNA methylation\", \n",
    "                           #platform = \"Illumina Human Methylation 450\", \n",
    "                           #legacy = TRUE)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "2fe58ed6",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:53:04.127740Z",
     "iopub.status.busy": "2022-12-28T13:53:04.125871Z",
     "iopub.status.idle": "2022-12-28T13:53:04.140865Z",
     "shell.execute_reply": "2022-12-28T13:53:04.138772Z"
    },
    "papermill": {
     "duration": 0.029437,
     "end_time": "2022-12-28T13:53:04.143573",
     "exception": false,
     "start_time": "2022-12-28T13:53:04.114136",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "#GDCdownload(query_met)\n",
    "#putting files together\n",
    "#data.met <- GDCprepare(query_met)\n",
    "#saving the met object\n",
    "#saveRDS(object = data.met,\n",
    "       # file = \"data.met.RDS\",\n",
    "      #  compress = FALSE)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "dc9d8fb7",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:53:04.167881Z",
     "iopub.status.busy": "2022-12-28T13:53:04.165926Z",
     "iopub.status.idle": "2022-12-28T13:54:03.014292Z",
     "shell.execute_reply": "2022-12-28T13:54:03.012501Z"
    },
    "papermill": {
     "duration": 58.863186,
     "end_time": "2022-12-28T13:54:03.017291",
     "exception": false,
     "start_time": "2022-12-28T13:53:04.154105",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "save_dir = \"/kaggle/working/\"\n",
    "\n",
    "# loading saved session: Once you saved your data, you can load it into your environment: \n",
    "data.met = readRDS(file = \"/kaggle/input/data-met-rds/data.met.RDS\")\n",
    "# met matrix\n",
    "met <- as.data.frame(SummarizedExperiment::assay(data.met))\n",
    "# clinical data\n",
    "clinical <- data.frame(data.met@colData)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "fbd73797",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:54:03.042638Z",
     "iopub.status.busy": "2022-12-28T13:54:03.041048Z",
     "iopub.status.idle": "2022-12-28T13:54:22.528391Z",
     "shell.execute_reply": "2022-12-28T13:54:22.525682Z"
    },
    "papermill": {
     "duration": 19.502431,
     "end_time": "2022-12-28T13:54:22.531867",
     "exception": false,
     "start_time": "2022-12-28T13:54:03.029436",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<strong>png:</strong> 2"
      ],
      "text/latex": [
       "\\textbf{png:} 2"
      ],
      "text/markdown": [
       "**png:** 2"
      ],
      "text/plain": [
       "png \n",
       "  2 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "### MDS ###\n",
    "pal <- brewer.pal(8,\"Dark2\")\n",
    "png(paste0(save_dir, \"MDS plot after GDCprepare.png\"), width=1000, height=750)\n",
    "plotMDS(met, top=1000, gene.selection=\"common\", \n",
    "        col=pal[factor(clinical$sample_type)])\n",
    "legend(\"topright\", legend=levels(factor(clinical$sample_type)), text.col=pal,\n",
    "       bg=\"white\", cex=0.7)\n",
    "dev.off()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "bda2f22c",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:54:22.555108Z",
     "iopub.status.busy": "2022-12-28T13:54:22.553561Z",
     "iopub.status.idle": "2022-12-28T13:54:53.371110Z",
     "shell.execute_reply": "2022-12-28T13:54:53.369254Z"
    },
    "papermill": {
     "duration": 30.840418,
     "end_time": "2022-12-28T13:54:53.382255",
     "exception": false,
     "start_time": "2022-12-28T13:54:22.541837",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<strong>png:</strong> 2"
      ],
      "text/latex": [
       "\\textbf{png:} 2"
      ],
      "text/markdown": [
       "**png:** 2"
      ],
      "text/plain": [
       "png \n",
       "  2 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "\n",
    "#density\n",
    "png(paste0(save_dir, \"densities_after_GDCprepare.png\"), width=1200, height=850)\n",
    "plotDensities(met, legend=FALSE, main= \"Densities sesame data\")\n",
    "dev.off() \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "id": "9ebf8146",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:54:53.405898Z",
     "iopub.status.busy": "2022-12-28T13:54:53.404310Z",
     "iopub.status.idle": "2022-12-28T13:55:01.146524Z",
     "shell.execute_reply": "2022-12-28T13:55:01.144775Z"
    },
    "papermill": {
     "duration": 7.756652,
     "end_time": "2022-12-28T13:55:01.148948",
     "exception": false,
     "start_time": "2022-12-28T13:54:53.392296",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<strong>png:</strong> 2"
      ],
      "text/latex": [
       "\\textbf{png:} 2"
      ],
      "text/markdown": [
       "**png:** 2"
      ],
      "text/plain": [
       "png \n",
       "  2 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "\n",
    "#############plot mean b values\n",
    "# remove probes with NA (similar to na.omit)\n",
    "tcga_na_filtered <- data.met[rowSums(is.na(assay(data.met))) == 0,]\n",
    "\n",
    "\n",
    "df <- data.frame(\n",
    "  \"Sample.mean\" = colMeans(assay(tcga_na_filtered), na.rm = TRUE),\n",
    "  \"groups\" = tcga_na_filtered$sample_type\n",
    ")\n",
    "\n",
    "\n",
    "library(\"ggpubr\")\n",
    "png(paste0(save_dir, \"mean_methylation_groups.png\"), width=1200, height=850)\n",
    "ggpubr::ggboxplot(\n",
    "  data = df,\n",
    "  y = \"Sample.mean\",\n",
    "  x = \"groups\",\n",
    "  color = \"groups\",\n",
    "  add = \"jitter\",\n",
    "  ylab = expression(paste(\"Mean DNA methylation (\", beta, \"-values)\")),\n",
    "  xlab = \"\"\n",
    ") + stat_compare_means() \n",
    "dev.off()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "d54de28f",
   "metadata": {
    "papermill": {
     "duration": 0.010335,
     "end_time": "2022-12-28T13:55:01.169664",
     "exception": false,
     "start_time": "2022-12-28T13:55:01.159329",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "----- **Inspecting Methylation data ** -------------"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "9c1ed5d0",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:55:01.193627Z",
     "iopub.status.busy": "2022-12-28T13:55:01.192028Z",
     "iopub.status.idle": "2022-12-28T13:55:03.524166Z",
     "shell.execute_reply": "2022-12-28T13:55:03.522464Z"
    },
    "papermill": {
     "duration": 2.346919,
     "end_time": "2022-12-28T13:55:03.526588",
     "exception": false,
     "start_time": "2022-12-28T13:55:01.179669",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "\n",
       " FALSE   TRUE \n",
       "104354 381223 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "#___________inspecting methylation data_______________#\n",
    "\n",
    "# get the 450k annotation data\n",
    "ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)\n",
    "\n",
    "## remove probes with NA\n",
    "probe.na <- rowSums(is.na(met))\n",
    "\n",
    "table(probe.na == 0)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "id": "07ae6444",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:55:03.551464Z",
     "iopub.status.busy": "2022-12-28T13:55:03.549922Z",
     "iopub.status.idle": "2022-12-28T13:55:09.665380Z",
     "shell.execute_reply": "2022-12-28T13:55:09.663694Z"
    },
    "papermill": {
     "duration": 6.130515,
     "end_time": "2022-12-28T13:55:09.667706",
     "exception": false,
     "start_time": "2022-12-28T13:55:03.537191",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<style>\n",
       ".list-inline {list-style: none; margin:0; padding: 0}\n",
       ".list-inline>li {display: inline-block}\n",
       ".list-inline>li:not(:last-child)::after {content: \"\\00b7\"; padding: 0 .5ex}\n",
       "</style>\n",
       "<ol class=list-inline><li>381223</li><li>412</li></ol>\n"
      ],
      "text/latex": [
       "\\begin{enumerate*}\n",
       "\\item 381223\n",
       "\\item 412\n",
       "\\end{enumerate*}\n"
      ],
      "text/markdown": [
       "1. 381223\n",
       "2. 412\n",
       "\n",
       "\n"
      ],
      "text/plain": [
       "[1] 381223    412"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "# chose those has no NA values in rows\n",
    "probe <- probe.na[probe.na == 0]\n",
    "met <- met[row.names(met) %in% names(probe), ]\n",
    "\n",
    "dim(met) # 314845    412\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "id": "08a9caee",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:55:09.692818Z",
     "iopub.status.busy": "2022-12-28T13:55:09.691291Z",
     "iopub.status.idle": "2022-12-28T13:55:12.122770Z",
     "shell.execute_reply": "2022-12-28T13:55:12.120956Z"
    },
    "papermill": {
     "duration": 2.447498,
     "end_time": "2022-12-28T13:55:12.125968",
     "exception": false,
     "start_time": "2022-12-28T13:55:09.678470",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "keep\n",
       " FALSE   TRUE \n",
       "  9068 372155 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "## remove probes that match chromosomes X and Y \n",
    "keep <- !(row.names(met) %in% ann450k$Name[ann450k$chr %in% c(\"chrX\",\"chrY\")])\n",
    "table(keep)\n",
    "met <- met[keep, ]\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "id": "b957a9f1",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:55:12.150515Z",
     "iopub.status.busy": "2022-12-28T13:55:12.149003Z",
     "iopub.status.idle": "2022-12-28T13:55:12.170215Z",
     "shell.execute_reply": "2022-12-28T13:55:12.168557Z"
    },
    "papermill": {
     "duration": 0.036711,
     "end_time": "2022-12-28T13:55:12.173353",
     "exception": false,
     "start_time": "2022-12-28T13:55:12.136642",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<style>\n",
       ".list-inline {list-style: none; margin:0; padding: 0}\n",
       ".list-inline>li {display: inline-block}\n",
       ".list-inline>li:not(:last-child)::after {content: \"\\00b7\"; padding: 0 .5ex}\n",
       "</style>\n",
       "<ol class=list-inline><li>372155</li><li>412</li></ol>\n"
      ],
      "text/latex": [
       "\\begin{enumerate*}\n",
       "\\item 372155\n",
       "\\item 412\n",
       "\\end{enumerate*}\n"
      ],
      "text/markdown": [
       "1. 372155\n",
       "2. 412\n",
       "\n",
       "\n"
      ],
      "text/plain": [
       "[1] 372155    412"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "\n",
    "dim(met) ##  308731    412\n",
    "rm(keep) # remove no further needed probes.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "id": "374580b5",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:55:12.204887Z",
     "iopub.status.busy": "2022-12-28T13:55:12.203317Z",
     "iopub.status.idle": "2022-12-28T13:55:12.431955Z",
     "shell.execute_reply": "2022-12-28T13:55:12.430107Z"
    },
    "papermill": {
     "duration": 0.250068,
     "end_time": "2022-12-28T13:55:12.434744",
     "exception": false,
     "start_time": "2022-12-28T13:55:12.184676",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "\n",
       " FALSE   TRUE \n",
       " 87018 398494 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "## remove SNPs overlapped probe\n",
    "table (is.na(ann450k$Probe_rs))\n",
    "# probes without snp\n",
    "no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]\n",
    "\n",
    "\n",
    "snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]\n",
    "#snps with maf <= 0.05\n",
    "snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "id": "39946b5d",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:55:12.460278Z",
     "iopub.status.busy": "2022-12-28T13:55:12.458718Z",
     "iopub.status.idle": "2022-12-28T13:55:17.871803Z",
     "shell.execute_reply": "2022-12-28T13:55:17.870076Z"
    },
    "papermill": {
     "duration": 5.428028,
     "end_time": "2022-12-28T13:55:17.874011",
     "exception": false,
     "start_time": "2022-12-28T13:55:12.445983",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<style>\n",
       ".list-inline {list-style: none; margin:0; padding: 0}\n",
       ".list-inline>li {display: inline-block}\n",
       ".list-inline>li:not(:last-child)::after {content: \"\\00b7\"; padding: 0 .5ex}\n",
       "</style>\n",
       "<ol class=list-inline><li>341803</li><li>412</li></ol>\n"
      ],
      "text/latex": [
       "\\begin{enumerate*}\n",
       "\\item 341803\n",
       "\\item 412\n",
       "\\end{enumerate*}\n"
      ],
      "text/markdown": [
       "1. 341803\n",
       "2. 412\n",
       "\n",
       "\n"
      ],
      "text/plain": [
       "[1] 341803    412"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "# filter met\n",
    "met <- met[row.names(met) %in% c(no.snp.probe, snp5.probe), ]\n",
    "\n",
    "dim(met) #283815    412\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "2cc2e1b6",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:55:17.900426Z",
     "iopub.status.busy": "2022-12-28T13:55:17.898867Z",
     "iopub.status.idle": "2022-12-28T13:55:17.988487Z",
     "shell.execute_reply": "2022-12-28T13:55:17.986528Z"
    },
    "papermill": {
     "duration": 0.104913,
     "end_time": "2022-12-28T13:55:17.990755",
     "exception": false,
     "start_time": "2022-12-28T13:55:17.885842",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "\n",
    "#remove no-further needed dataset\n",
    "rm(no.snp.probe, probe, probe.na, snp.probe, snp5.probe)\n",
    "\n",
    "## Removing probes that have been demonstrated to map to multiple places in the genome.\n",
    "# list adapted from https://www.tandfonline.com/doi/full/10.4161/epi.23470\n",
    "\n",
    "crs.reac <- read.csv(\"/kaggle/input/data-met-rds/cross_reactive_probe.chen2013.csv\")\n",
    "crs.reac <- crs.reac$TargetID[-1]\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "id": "567c8982",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:55:18.016504Z",
     "iopub.status.busy": "2022-12-28T13:55:18.014994Z",
     "iopub.status.idle": "2022-12-28T13:55:43.861034Z",
     "shell.execute_reply": "2022-12-28T13:55:43.859317Z"
    },
    "papermill": {
     "duration": 25.861376,
     "end_time": "2022-12-28T13:55:43.863299",
     "exception": false,
     "start_time": "2022-12-28T13:55:18.001923",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<style>\n",
       ".list-inline {list-style: none; margin:0; padding: 0}\n",
       ".list-inline>li {display: inline-block}\n",
       ".list-inline>li:not(:last-child)::after {content: \"\\00b7\"; padding: 0 .5ex}\n",
       "</style>\n",
       "<ol class=list-inline><li>329304</li><li>412</li></ol>\n"
      ],
      "text/latex": [
       "\\begin{enumerate*}\n",
       "\\item 329304\n",
       "\\item 412\n",
       "\\end{enumerate*}\n"
      ],
      "text/markdown": [
       "1. 329304\n",
       "2. 412\n",
       "\n",
       "\n"
      ],
      "text/plain": [
       "[1] 329304    412"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "text/html": [
       "<strong>png:</strong> 2"
      ],
      "text/latex": [
       "\\textbf{png:} 2"
      ],
      "text/markdown": [
       "**png:** 2"
      ],
      "text/plain": [
       "png \n",
       "  2 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "# filter met\n",
    "met <- met[ -which(row.names(met) %in% crs.reac), ]\n",
    "dim(met)\n",
    "\n",
    "\n",
    "#density\n",
    "png(paste0(save_dir, \"densities_after_filtering.png\"), width=1200, height=850)\n",
    "plotDensities(met, legend=FALSE, main= \"Densities after filtering data\")\n",
    "dev.off() \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "38150a19",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:55:43.889941Z",
     "iopub.status.busy": "2022-12-28T13:55:43.888426Z",
     "iopub.status.idle": "2022-12-28T13:56:17.197954Z",
     "shell.execute_reply": "2022-12-28T13:56:17.195991Z"
    },
    "papermill": {
     "duration": 33.325207,
     "end_time": "2022-12-28T13:56:17.200264",
     "exception": false,
     "start_time": "2022-12-28T13:55:43.875057",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "##################save files#################################################\n",
    "#save as bval\n",
    "bval <- met\n",
    "\n",
    "## converting beta values to m_values\n",
    "## m = log2(beta/1-beta)\n",
    "mval <- t(apply(met, 1, function(x) log2(x/(1-x))))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "id": "9e4a7dbe",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:56:17.226941Z",
     "iopub.status.busy": "2022-12-28T13:56:17.225380Z",
     "iopub.status.idle": "2022-12-28T13:56:23.703648Z",
     "shell.execute_reply": "2022-12-28T13:56:23.701794Z"
    },
    "papermill": {
     "duration": 6.495451,
     "end_time": "2022-12-28T13:56:23.707379",
     "exception": false,
     "start_time": "2022-12-28T13:56:17.211928",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "#______________saving/loading_____________________#\n",
    "# save data sets\n",
    "\n",
    "saveRDS(mval, file = \"mval.RDS\", compress = FALSE)\n",
    "saveRDS (bval, file = \"bval.RDS\", compress = FALSE)\n",
    "\n",
    "\n",
    "#mval <- readRDS(\"mval.RDS\")\n",
    "#bval <- readRDS(\"bval.RDS\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 19,
   "id": "df3eec3c",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:56:23.824412Z",
     "iopub.status.busy": "2022-12-28T13:56:23.821053Z",
     "iopub.status.idle": "2022-12-28T13:56:48.609529Z",
     "shell.execute_reply": "2022-12-28T13:56:48.607892Z"
    },
    "papermill": {
     "duration": 24.856571,
     "end_time": "2022-12-28T13:56:48.612256",
     "exception": false,
     "start_time": "2022-12-28T13:56:23.755685",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<strong>png:</strong> 2"
      ],
      "text/latex": [
       "\\textbf{png:} 2"
      ],
      "text/markdown": [
       "**png:** 2"
      ],
      "text/plain": [
       "png \n",
       "  2 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "\n",
    "#densities Bval\n",
    "png(paste0(save_dir, \"densities_bval.png\"), width=1200, height=850)\n",
    "plotDensities(bval, legend=FALSE, main= \"Densities Bval data\")\n",
    "dev.off()\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 20,
   "id": "4a397917",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:56:48.640965Z",
     "iopub.status.busy": "2022-12-28T13:56:48.639461Z",
     "iopub.status.idle": "2022-12-28T13:57:09.383271Z",
     "shell.execute_reply": "2022-12-28T13:57:09.381481Z"
    },
    "papermill": {
     "duration": 20.760365,
     "end_time": "2022-12-28T13:57:09.385795",
     "exception": false,
     "start_time": "2022-12-28T13:56:48.625430",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<strong>png:</strong> 2"
      ],
      "text/latex": [
       "\\textbf{png:} 2"
      ],
      "text/markdown": [
       "**png:** 2"
      ],
      "text/plain": [
       "png \n",
       "  2 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "\n",
    "#densities mval\n",
    "png(paste0(save_dir, \"densities_mval.png\"), width=1200, height=850)\n",
    "plotDensities(mval, legend=FALSE, main= \"Densities Mval data\")\n",
    "dev.off()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "08efcc70",
   "metadata": {
    "papermill": {
     "duration": 0.013173,
     "end_time": "2022-12-28T13:57:09.411788",
     "exception": false,
     "start_time": "2022-12-28T13:57:09.398615",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    " ---------- **DIFFERENTIAL METHYLATION ANALYSIS** -------------------"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "id": "8bed0737",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:57:09.442179Z",
     "iopub.status.busy": "2022-12-28T13:57:09.440676Z",
     "iopub.status.idle": "2022-12-28T13:57:09.460873Z",
     "shell.execute_reply": "2022-12-28T13:57:09.459120Z"
    },
    "papermill": {
     "duration": 0.038677,
     "end_time": "2022-12-28T13:57:09.463693",
     "exception": false,
     "start_time": "2022-12-28T13:57:09.425016",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "                     \n",
       "                      Lower lobe, lung Lung, NOS Main bronchus\n",
       "  Primary Tumor                    128        19             5\n",
       "  Solid Tissue Normal               17         3             0\n",
       "                     \n",
       "                      Middle lobe, lung Overlapping lesion of lung\n",
       "  Primary Tumor                      11                          8\n",
       "  Solid Tissue Normal                 0                          1\n",
       "                     \n",
       "                      Upper lobe, lung\n",
       "  Primary Tumor                    199\n",
       "  Solid Tissue Normal               21"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "\n",
    "##############--------- Differential methylation analysis -------#########\n",
    "table(clinical$sample_type, clinical$tissue_or_organ_of_origin)\n",
    "\n",
    "##                     Lower lobe, lung Lung, NOS Main bronchus Middle lobe, lung Overlapping lesion of lung Upper lobe, lung\n",
    "#Primary Tumor                    128        19             5                11                          8              199\n",
    "#Solid Tissue Normal               17         3             0                 0                          1               21\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "id": "42785e04",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:57:09.491575Z",
     "iopub.status.busy": "2022-12-28T13:57:09.490092Z",
     "iopub.status.idle": "2022-12-28T13:57:09.884176Z",
     "shell.execute_reply": "2022-12-28T13:57:09.882435Z"
    },
    "papermill": {
     "duration": 0.410606,
     "end_time": "2022-12-28T13:57:09.886562",
     "exception": false,
     "start_time": "2022-12-28T13:57:09.475956",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "\n",
       "TRUE \n",
       " 412 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "text/plain": [
       "\n",
       "TRUE \n",
       " 412 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "text/html": [
       "TRUE"
      ],
      "text/latex": [
       "TRUE"
      ],
      "text/markdown": [
       "TRUE"
      ],
      "text/plain": [
       "[1] TRUE"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "text/html": [
       "TRUE"
      ],
      "text/latex": [
       "TRUE"
      ],
      "text/markdown": [
       "TRUE"
      ],
      "text/plain": [
       "[1] TRUE"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "\n",
    "# filtering and grouping\n",
    "clinical <- clinical[, c(\"barcode\", \"sample_type\", \"site_of_resection_or_biopsy\")]\n",
    "clinical <- na.omit(clinical)\n",
    "barcode <- clinical$barcode\n",
    "\n",
    "# removing samples from meth matrixes\n",
    "bval <- bval[, colnames(bval) %in% barcode]\n",
    "mval <- mval[, colnames(mval) %in% barcode]\n",
    "\n",
    "# Making sure about samples in clinical and matrixes and their order\n",
    "table(colnames(mval) %in% row.names(clinical))\n",
    "table(colnames(bval) %in% row.names(clinical))\n",
    "#\n",
    "all(row.names(clinical) == colnames(bval))\n",
    "all(row.names(clinical) == colnames(mval))\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 23,
   "id": "981d1114",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:57:09.916608Z",
     "iopub.status.busy": "2022-12-28T13:57:09.914939Z",
     "iopub.status.idle": "2022-12-28T13:57:09.934048Z",
     "shell.execute_reply": "2022-12-28T13:57:09.932103Z"
    },
    "papermill": {
     "duration": 0.03708,
     "end_time": "2022-12-28T13:57:09.936561",
     "exception": false,
     "start_time": "2022-12-28T13:57:09.899481",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "\n",
    "#Making grouping variable\n",
    "clinical$sample_type <- as.factor(clinical$sample_type)\n",
    "#levels(clinical$sample_type)\n",
    "clinical$sample_type <- relevel(clinical$sample_type, ref = \"Solid Tissue Normal\")\n",
    "\n",
    "\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "5b910e12",
   "metadata": {
    "papermill": {
     "duration": 0.012928,
     "end_time": "2022-12-28T13:57:09.962980",
     "exception": false,
     "start_time": "2022-12-28T13:57:09.950052",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "----------- **DMC ANALYSIS**"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 24,
   "id": "098295ac",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:57:09.992985Z",
     "iopub.status.busy": "2022-12-28T13:57:09.991455Z",
     "iopub.status.idle": "2022-12-28T13:57:21.664514Z",
     "shell.execute_reply": "2022-12-28T13:57:21.662688Z"
    },
    "papermill": {
     "duration": 11.690993,
     "end_time": "2022-12-28T13:57:21.667311",
     "exception": false,
     "start_time": "2022-12-28T13:57:09.976318",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "                             (Intercept) sample_typePrimary Tumor\n",
      "TCGA-33-6737-01A-11D-1818-05           1                        1\n",
      "TCGA-O2-A52S-01A-11D-A25R-05           1                        1\n",
      "TCGA-77-7338-01A-11D-2043-05           1                        1\n",
      "TCGA-98-A539-01A-31D-A25R-05           1                        1\n",
      "TCGA-63-A5MT-01A-21D-A26N-05           1                        1\n",
      "TCGA-LA-A7SW-01A-11D-A408-05           1                        1\n",
      "TCGA-21-5786-01A-01D-1633-05           1                        1\n",
      "TCGA-85-A513-01A-12D-A26N-05           1                        1\n",
      "TCGA-85-A50M-01A-21D-A25R-05           1                        1\n",
      "TCGA-63-A5MR-01A-31D-A27L-05           1                        1\n",
      "TCGA-NK-A5CR-01A-11D-A26N-05           1                        1\n",
      "TCGA-22-4605-01A-21D-2123-05           1                        1\n",
      "TCGA-43-3920-11B-01D-1551-05           1                        0\n",
      "TCGA-52-7811-01A-11D-2123-05           1                        1\n",
      "TCGA-63-A5MV-01A-21D-A26N-05           1                        1\n",
      "TCGA-85-6798-01A-11D-1947-05           1                        1\n",
      "TCGA-85-7843-01A-11D-2123-05           1                        1\n",
      "TCGA-85-A510-01A-11D-A26N-05           1                        1\n",
      "TCGA-56-8625-01A-11D-2398-05           1                        1\n",
      "TCGA-NK-A7XE-01A-12D-A408-05           1                        1\n",
      "TCGA-39-5019-11A-01D-1818-05           1                        0\n",
      "TCGA-NC-A5HD-01A-11D-A26N-05           1                        1\n",
      "TCGA-22-5477-11A-11D-1633-05           1                        0\n",
      "TCGA-O2-A5IB-01A-11D-A27L-05           1                        1\n",
      "TCGA-22-4593-01A-21D-1818-05           1                        1\n",
      "TCGA-92-7341-01A-31D-2043-05           1                        1\n",
      "TCGA-NK-A5CT-01A-31D-A26N-05           1                        1\n",
      "TCGA-22-5472-01A-01D-1633-05           1                        1\n",
      "TCGA-39-5029-11A-01D-1440-05           1                        0\n",
      "TCGA-39-5016-11A-01D-1440-05           1                        0\n",
      "TCGA-56-8629-01A-11D-2398-05           1                        1\n",
      "TCGA-85-A4CN-01A-11D-A24I-05           1                        1\n",
      "TCGA-56-8305-01A-11D-2294-05           1                        1\n",
      "TCGA-22-5474-11A-01D-1633-05           1                        0\n",
      "TCGA-56-7580-01A-11D-2043-05           1                        1\n",
      "TCGA-56-A4BX-01A-11D-A24I-05           1                        1\n",
      "TCGA-85-A4QR-01A-11D-A258-05           1                        1\n",
      "TCGA-NC-A5HH-01A-11D-A26N-05           1                        1\n",
      "TCGA-34-5240-01A-01D-1440-05           1                        1\n",
      "TCGA-63-A5MH-01A-12D-A27L-05           1                        1\n",
      "TCGA-94-8491-01A-11D-2324-05           1                        1\n",
      "TCGA-18-5595-01A-01D-1633-05           1                        1\n",
      "TCGA-63-A5MU-01A-11D-A26N-05           1                        1\n",
      "TCGA-NK-A5CX-01A-11D-A26N-05           1                        1\n",
      "TCGA-58-A46L-01A-11D-A24I-05           1                        1\n",
      "TCGA-NC-A5HM-01A-12D-A26N-05           1                        1\n",
      "TCGA-85-8664-01A-11D-2398-05           1                        1\n",
      "TCGA-90-7769-01A-11D-2123-05           1                        1\n",
      "TCGA-92-8064-01A-11D-2245-05           1                        1\n",
      "TCGA-39-5034-01A-01D-1440-05           1                        1\n",
      "TCGA-J1-A4AH-01A-31D-A24I-05           1                        1\n",
      "TCGA-43-6770-01A-11D-1818-05           1                        1\n",
      "TCGA-94-8490-01A-11D-2324-05           1                        1\n",
      "TCGA-96-A4JL-01A-11D-A258-05           1                        1\n",
      "TCGA-33-AASB-01A-11D-A408-05           1                        1\n",
      "TCGA-77-7335-01A-11D-2043-05           1                        1\n",
      "TCGA-43-5670-01A-21D-2123-05           1                        1\n",
      "TCGA-90-7766-01A-21D-2123-05           1                        1\n",
      "TCGA-43-8116-01A-11D-2245-05           1                        1\n",
      "TCGA-22-A5C4-01A-12D-A27L-05           1                        1\n",
      "TCGA-63-A5MM-01A-11D-A26N-05           1                        1\n",
      "TCGA-56-A62T-01A-11D-A408-05           1                        1\n",
      "TCGA-77-A5FZ-01A-31D-A27L-05           1                        1\n",
      "TCGA-43-A56U-01A-11D-A26N-05           1                        1\n",
      "TCGA-98-A53C-01A-11D-A25R-05           1                        1\n",
      "TCGA-77-8131-01A-11D-2245-05           1                        1\n",
      "TCGA-58-A46J-01A-11D-A24I-05           1                        1\n",
      "TCGA-43-8115-01A-11D-2245-05           1                        1\n",
      "TCGA-22-4613-11A-01D-1440-05           1                        0\n",
      "TCGA-22-5492-01A-01D-1633-05           1                        1\n",
      "TCGA-NC-A5HK-01A-11D-A26N-05           1                        1\n",
      "TCGA-21-5784-01A-01D-1633-05           1                        1\n",
      "TCGA-63-7022-01A-11D-1947-05           1                        1\n",
      "TCGA-56-A4ZJ-01A-11D-A25R-05           1                        1\n",
      "TCGA-58-A46N-01A-11D-A24I-05           1                        1\n",
      "TCGA-56-5898-01A-11D-1633-05           1                        1\n",
      "TCGA-63-6202-01A-11D-1818-05           1                        1\n",
      "TCGA-22-5481-01A-31D-1947-05           1                        1\n",
      "TCGA-22-4599-11A-01D-1440-05           1                        0\n",
      "TCGA-O2-A52Q-01A-11D-A26N-05           1                        1\n",
      "TCGA-51-6867-01A-11D-2043-05           1                        1\n",
      "TCGA-18-3417-01A-01D-1440-05           1                        1\n",
      "TCGA-O2-A52W-01A-11D-A26N-05           1                        1\n",
      "TCGA-77-8009-01A-11D-2185-05           1                        1\n",
      "TCGA-34-5239-01A-21D-1818-05           1                        1\n",
      "TCGA-43-7657-01A-31D-2123-05           1                        1\n",
      "TCGA-85-8352-01A-31D-2324-05           1                        1\n",
      "TCGA-56-8622-01A-11D-2398-05           1                        1\n",
      "TCGA-98-8022-01A-11D-2245-05           1                        1\n",
      "TCGA-68-A59I-01A-11D-A26N-05           1                        1\n",
      "TCGA-85-8580-01A-31D-2398-05           1                        1\n",
      "TCGA-98-A53D-01A-32D-A26N-05           1                        1\n",
      "TCGA-33-4582-11A-01D-1440-05           1                        0\n",
      "TCGA-46-6026-01A-11D-1818-05           1                        1\n",
      "TCGA-63-A5MP-01A-11D-A26N-05           1                        1\n",
      "TCGA-52-7622-01A-11D-2123-05           1                        1\n",
      "TCGA-39-5036-11A-01D-1440-05           1                        0\n",
      "TCGA-39-5027-01A-21D-1818-05           1                        1\n",
      "TCGA-92-8065-01A-11D-2245-05           1                        1\n",
      "TCGA-39-5011-11A-01D-1440-05           1                        0\n",
      "TCGA-33-4582-01A-01D-1440-05           1                        1\n",
      "TCGA-77-A5GA-01A-11D-A27L-05           1                        1\n",
      "TCGA-37-A5EL-01A-11D-A26N-05           1                        1\n",
      "TCGA-34-5927-01A-11D-1818-05           1                        1\n",
      "TCGA-77-8144-01A-11D-2245-05           1                        1\n",
      "TCGA-85-6560-01A-11D-1818-05           1                        1\n",
      "TCGA-77-8154-01A-11D-2245-05           1                        1\n",
      "TCGA-39-5036-01A-01D-1440-05           1                        1\n",
      "TCGA-56-7582-01A-11D-2043-05           1                        1\n",
      "TCGA-43-6143-01A-11D-1818-05           1                        1\n",
      "TCGA-98-A53J-01A-11D-A26N-05           1                        1\n",
      "TCGA-43-8118-01A-11D-2398-05           1                        1\n",
      "TCGA-22-5489-01A-01D-1633-05           1                        1\n",
      "TCGA-21-5787-01A-01D-1633-05           1                        1\n",
      "TCGA-63-5128-01A-01D-1440-05           1                        1\n",
      "TCGA-56-A4BW-01A-11D-A24I-05           1                        1\n",
      "TCGA-O2-A52V-01A-31D-A25R-05           1                        1\n",
      "TCGA-33-4583-01A-01D-1440-05           1                        1\n",
      "TCGA-58-8390-01A-11D-2324-05           1                        1\n",
      "TCGA-39-5035-11A-01D-1440-05           1                        0\n",
      "TCGA-77-A5G8-01B-11D-A27L-05           1                        1\n",
      "TCGA-56-7221-01A-11D-2043-05           1                        1\n",
      "TCGA-39-5040-01A-21D-2123-05           1                        1\n",
      "TCGA-58-8388-01A-11D-2324-05           1                        1\n",
      "TCGA-22-5489-11A-01D-1633-05           1                        0\n",
      "TCGA-77-8138-01A-11D-2245-05           1                        1\n",
      "TCGA-56-8307-01A-11D-2294-05           1                        1\n",
      "TCGA-85-8052-01A-11D-2245-05           1                        1\n",
      "TCGA-39-5022-01A-21D-1818-05           1                        1\n",
      "TCGA-63-A5MN-01A-22D-A27L-05           1                        1\n",
      "TCGA-18-5595-11A-01D-1633-05           1                        0\n",
      "TCGA-77-A5G6-01A-11D-A27L-05           1                        1\n",
      "TCGA-39-5037-01A-01D-1440-05           1                        1\n",
      "TCGA-NC-A5HE-01A-11D-A26N-05           1                        1\n",
      "TCGA-18-4721-11A-01D-1440-05           1                        0\n",
      "TCGA-96-7545-01A-21D-2043-05           1                        1\n",
      "TCGA-34-8454-01A-11D-2324-05           1                        1\n",
      "TCGA-77-7139-01A-11D-2043-05           1                        1\n",
      "TCGA-18-4721-01A-01D-1440-05           1                        1\n",
      "TCGA-98-7454-01A-11D-2043-05           1                        1\n",
      "TCGA-22-5482-01A-01D-1633-05           1                        1\n",
      "TCGA-22-5491-01A-01D-1633-05           1                        1\n",
      "TCGA-34-A5IX-01A-12D-A27L-05           1                        1\n",
      "TCGA-NC-A5HF-01A-11D-A26N-05           1                        1\n",
      "TCGA-92-8063-01A-11D-2245-05           1                        1\n",
      "TCGA-43-6771-01A-11D-1818-05           1                        1\n",
      "TCGA-85-A4QQ-01A-41D-A25R-05           1                        1\n",
      "TCGA-77-A5GH-01A-11D-A27L-05           1                        1\n",
      "TCGA-56-8628-01A-11D-2398-05           1                        1\n",
      "TCGA-85-A4CL-01A-41D-A26N-05           1                        1\n",
      "TCGA-43-7658-01A-11D-2123-05           1                        1\n",
      "TCGA-77-8130-01A-11D-2245-05           1                        1\n",
      "TCGA-NC-A5HP-01A-11D-A26N-05           1                        1\n",
      "TCGA-77-8153-01A-11D-2398-05           1                        1\n",
      "TCGA-22-5491-11A-01D-1633-05           1                        0\n",
      "TCGA-56-7730-01A-11D-2123-05           1                        1\n",
      "TCGA-63-A5MS-01A-11D-A26N-05           1                        1\n",
      "TCGA-22-4609-01A-21D-2123-05           1                        1\n",
      "TCGA-56-A49D-01A-11D-A24I-05           1                        1\n",
      "TCGA-63-A5ML-01A-31D-A27L-05           1                        1\n",
      "TCGA-52-7809-01A-21D-2123-05           1                        1\n",
      "TCGA-33-6738-01A-11D-1947-05           1                        1\n",
      "TCGA-34-5236-01A-21D-1818-05           1                        1\n",
      "TCGA-85-8288-01A-11D-2294-05           1                        1\n",
      "TCGA-43-A474-01A-11D-A24I-05           1                        1\n",
      "TCGA-85-8049-01A-11D-2245-05           1                        1\n",
      "TCGA-70-6723-01A-11D-1818-05           1                        1\n",
      "TCGA-56-8309-01A-11D-2294-05           1                        1\n",
      "TCGA-77-8008-01A-21D-2185-05           1                        1\n",
      "TCGA-56-6546-01A-11D-1818-05           1                        1\n",
      "TCGA-21-A5DI-01A-31D-A26N-05           1                        1\n",
      "TCGA-39-5028-11A-01D-1440-05           1                        0\n",
      "TCGA-56-8083-01A-11D-2245-05           1                        1\n",
      "TCGA-85-8287-01A-11D-2294-05           1                        1\n",
      "TCGA-22-5483-01A-01D-1818-05           1                        1\n",
      "TCGA-56-7822-01A-11D-2123-05           1                        1\n",
      "TCGA-33-4586-11A-01D-1440-05           1                        0\n",
      "TCGA-63-A5MJ-01A-11D-A27L-05           1                        1\n",
      "TCGA-63-A5M9-01A-11D-A26N-05           1                        1\n",
      "TCGA-63-7021-01A-11D-1947-05           1                        1\n",
      "TCGA-22-5485-11A-01D-1633-05           1                        0\n",
      "TCGA-77-6845-01A-11D-1947-05           1                        1\n",
      "TCGA-39-5011-01A-01D-1440-05           1                        1\n",
      "TCGA-34-5929-01A-11D-1818-05           1                        1\n",
      "TCGA-33-AASI-01A-22D-A408-05           1                        1\n",
      "TCGA-85-8277-01A-11D-2294-05           1                        1\n",
      "TCGA-98-A53I-01A-31D-A25R-05           1                        1\n",
      "TCGA-85-8048-01A-11D-2245-05           1                        1\n",
      "TCGA-85-8070-01A-11D-2245-05           1                        1\n",
      "TCGA-21-5782-01A-01D-1633-05           1                        1\n",
      "TCGA-77-8143-01A-11D-2245-05           1                        1\n",
      "TCGA-85-A511-01A-21D-A25R-05           1                        1\n",
      "TCGA-56-8504-01A-11D-2398-05           1                        1\n",
      "TCGA-85-8276-01A-11D-2294-05           1                        1\n",
      "TCGA-39-5024-01A-21D-1818-05           1                        1\n",
      "TCGA-98-8023-01A-11D-2245-05           1                        1\n",
      "TCGA-43-5668-11A-01D-1633-05           1                        0\n",
      "TCGA-98-8020-01A-11D-2245-05           1                        1\n",
      "TCGA-98-8021-01A-11D-2245-05           1                        1\n",
      "TCGA-43-6773-01A-41D-1947-05           1                        1\n",
      "TCGA-98-A53A-01A-11D-A25R-05           1                        1\n",
      "TCGA-77-6844-01A-11D-1947-05           1                        1\n",
      "TCGA-56-7222-01A-11D-2043-05           1                        1\n",
      "TCGA-39-5030-01A-01D-1440-05           1                        1\n",
      "TCGA-56-A5DR-01A-11D-A27L-05           1                        1\n",
      "TCGA-22-5478-01A-01D-1633-05           1                        1\n",
      "TCGA-85-A4PA-01A-11D-A258-05           1                        1\n",
      "TCGA-NC-A5HI-01A-11D-A26N-05           1                        1\n",
      "TCGA-34-5232-01A-21D-1818-05           1                        1\n",
      "TCGA-77-A5G7-01B-11D-A27L-05           1                        1\n",
      "TCGA-33-AAS8-01A-11D-A408-05           1                        1\n",
      "TCGA-58-A46M-01A-11D-A24I-05           1                        1\n",
      "TCGA-85-8071-01A-11D-2245-05           1                        1\n",
      "TCGA-56-8082-01A-11D-2245-05           1                        1\n",
      "TCGA-58-8387-01A-11D-2294-05           1                        1\n",
      "TCGA-NC-A5HN-01A-11D-A26N-05           1                        1\n",
      "TCGA-85-A50Z-01A-21D-A25R-05           1                        1\n",
      "TCGA-34-5234-01A-01D-1633-05           1                        1\n",
      "TCGA-39-5031-11A-01D-1440-05           1                        0\n",
      "TCGA-96-8169-01A-11D-2294-05           1                        1\n",
      "TCGA-22-5471-01A-01D-1633-05           1                        1\n",
      "TCGA-96-A4JK-01A-11D-A258-05           1                        1\n",
      "TCGA-77-8133-01A-12D-2245-05           1                        1\n",
      "TCGA-56-8304-01A-11D-2324-05           1                        1\n",
      "TCGA-63-A5MB-01A-11D-A26N-05           1                        1\n",
      "TCGA-34-5928-01A-11D-1818-05           1                        1\n",
      "TCGA-63-A5MY-01A-11D-A26N-05           1                        1\n",
      "TCGA-56-6545-01A-11D-1818-05           1                        1\n",
      "TCGA-58-A46K-01A-11D-A24I-05           1                        1\n",
      "TCGA-85-7950-01A-11D-2185-05           1                        1\n",
      "TCGA-58-8393-01A-11D-2324-05           1                        1\n",
      "TCGA-63-A5MI-01A-12D-A27L-05           1                        1\n",
      "TCGA-77-8156-01A-11D-2245-05           1                        1\n",
      "TCGA-94-7943-01A-11D-2185-05           1                        1\n",
      "TCGA-77-8136-01A-11D-2245-05           1                        1\n",
      "TCGA-85-A512-01A-11D-A26N-05           1                        1\n",
      "TCGA-85-A53L-01A-21D-A26N-05           1                        1\n",
      "TCGA-85-8353-01A-21D-2294-05           1                        1\n",
      "TCGA-56-7731-01A-11D-2123-05           1                        1\n",
      "TCGA-85-7698-01A-11D-2123-05           1                        1\n",
      "TCGA-NC-A5HO-01A-11D-A26N-05           1                        1\n",
      "TCGA-85-A4JB-01A-51D-A25R-05           1                        1\n",
      "TCGA-98-A53B-01A-11D-A25R-05           1                        1\n",
      "TCGA-85-8481-01A-11D-2324-05           1                        1\n",
      "TCGA-22-5492-11A-01D-1633-05           1                        0\n",
      "TCGA-34-8456-01A-21D-2324-05           1                        1\n",
      "TCGA-90-A59Q-01A-11D-A26N-05           1                        1\n",
      "TCGA-77-8139-01A-11D-2245-05           1                        1\n",
      "TCGA-22-5485-01A-01D-1633-05           1                        1\n",
      "TCGA-77-8007-01A-11D-2185-05           1                        1\n",
      "TCGA-63-A5MG-01A-12D-A27L-05           1                        1\n",
      "TCGA-85-8582-01A-21D-2398-05           1                        1\n",
      "TCGA-98-A53H-01A-12D-A25R-05           1                        1\n",
      "TCGA-39-5039-11A-01D-1440-05           1                        0\n",
      "TCGA-33-A4WN-01A-11D-A25R-05           1                        1\n",
      "TCGA-58-8392-01A-11D-2324-05           1                        1\n",
      "TCGA-NC-A5HG-01A-11D-A26N-05           1                        1\n",
      "TCGA-85-8351-01A-11D-2294-05           1                        1\n",
      "TCGA-77-A5G1-01A-11D-A27L-05           1                        1\n",
      "TCGA-77-A5G3-01A-31D-A27L-05           1                        1\n",
      "TCGA-90-7767-01A-11D-2123-05           1                        1\n",
      "TCGA-92-7340-01A-21D-2043-05           1                        1\n",
      "TCGA-63-A5MW-01A-11D-A26N-05           1                        1\n",
      "TCGA-37-5819-01A-01D-1633-05           1                        1\n",
      "TCGA-34-5241-01A-01D-1440-05           1                        1\n",
      "TCGA-43-6771-11A-01D-1818-05           1                        0\n",
      "TCGA-68-7755-01A-11D-2123-05           1                        1\n",
      "TCGA-33-4586-01A-01D-1440-05           1                        1\n",
      "TCGA-56-8308-01A-11D-2294-05           1                        1\n",
      "TCGA-85-6561-01A-11D-1818-05           1                        1\n",
      "TCGA-18-3417-11A-01D-1440-05           1                        0\n",
      "TCGA-85-7699-01A-11D-2123-05           1                        1\n",
      "TCGA-77-A5GB-01B-11D-A27L-05           1                        1\n",
      "TCGA-NK-A5D1-01A-11D-A26N-05           1                        1\n",
      "TCGA-21-5783-01A-41D-2185-05           1                        1\n",
      "TCGA-33-AASL-01A-11D-A408-05           1                        1\n",
      "TCGA-63-7023-01A-11D-1947-05           1                        1\n",
      "TCGA-68-8250-01A-11D-2294-05           1                        1\n",
      "TCGA-85-7697-01A-11D-2123-05           1                        1\n",
      "TCGA-77-7138-01A-41D-2043-05           1                        1\n",
      "TCGA-NC-A5HL-01A-11D-A26N-05           1                        1\n",
      "TCGA-77-8150-01A-11D-2245-05           1                        1\n",
      "TCGA-96-8170-01A-11D-2294-05           1                        1\n",
      "TCGA-56-A4ZK-01A-11D-A25R-05           1                        1\n",
      "TCGA-85-8584-01A-11D-2398-05           1                        1\n",
      "TCGA-NC-A5HJ-01A-11D-A26N-05           1                        1\n",
      "TCGA-77-8146-01A-11D-2245-05           1                        1\n",
      "TCGA-52-7810-01A-11D-2123-05           1                        1\n",
      "TCGA-22-5479-01A-31D-1947-05           1                        1\n",
      "TCGA-85-A4JC-01A-11D-A258-05           1                        1\n",
      "TCGA-90-A4EE-01A-11D-A258-05           1                        1\n",
      "TCGA-77-8140-01A-11D-2245-05           1                        1\n",
      "TCGA-85-8072-01A-31D-2245-05           1                        1\n",
      "TCGA-6A-AB49-01A-12D-A408-05           1                        1\n",
      "TCGA-L3-A524-01A-11D-A25R-05           1                        1\n",
      "TCGA-56-8201-01A-11D-2245-05           1                        1\n",
      "TCGA-77-7465-01A-11D-2043-05           1                        1\n",
      "TCGA-77-7142-01A-11D-2043-05           1                        1\n",
      "TCGA-68-7756-01A-11D-2123-05           1                        1\n",
      "TCGA-60-2697-01A-11D-2123-05           1                        1\n",
      "TCGA-85-7696-01A-11D-2123-05           1                        1\n",
      "TCGA-22-4613-01A-01D-1440-05           1                        1\n",
      "TCGA-77-7337-01A-21D-2043-05           1                        1\n",
      "TCGA-56-8503-01A-11D-2398-05           1                        1\n",
      "TCGA-33-AASJ-01A-11D-A408-05           1                        1\n",
      "TCGA-39-5039-01A-01D-1440-05           1                        1\n",
      "TCGA-60-2704-01A-11D-2043-05           1                        1\n",
      "TCGA-56-5897-01A-11D-1633-05           1                        1\n",
      "TCGA-77-7140-01A-41D-2043-05           1                        1\n",
      "TCGA-37-A5EN-01A-21D-A26N-05           1                        1\n",
      "TCGA-39-5021-11A-01D-1440-05           1                        0\n",
      "TCGA-77-8145-01A-11D-2245-05           1                        1\n",
      "TCGA-NC-A5HQ-01A-11D-A26N-05           1                        1\n",
      "TCGA-56-7223-01A-11D-2043-05           1                        1\n",
      "TCGA-56-A4BY-01A-11D-A24I-05           1                        1\n",
      "TCGA-85-6175-01A-11D-1818-05           1                        1\n",
      "TCGA-94-7033-01A-11D-1947-05           1                        1\n",
      "TCGA-43-5668-01A-01D-1633-05           1                        1\n",
      "TCGA-39-5030-11A-01D-1440-05           1                        0\n",
      "TCGA-85-8355-01A-11D-2294-05           1                        1\n",
      "TCGA-34-8455-01A-11D-2324-05           1                        1\n",
      "TCGA-56-7823-01B-11D-2245-05           1                        1\n",
      "TCGA-77-8128-01A-11D-2245-05           1                        1\n",
      "TCGA-85-8350-01A-11D-2294-05           1                        1\n",
      "TCGA-39-5034-11A-01D-1440-05           1                        0\n",
      "TCGA-34-5929-11A-01D-1818-05           1                        0\n",
      "TCGA-39-5031-01A-01D-1440-05           1                        1\n",
      "TCGA-60-2709-01A-21D-1818-05           1                        1\n",
      "TCGA-56-8624-01A-11D-2398-05           1                        1\n",
      "TCGA-77-8148-01A-11D-2245-05           1                        1\n",
      "TCGA-70-6722-01A-11D-1818-05           1                        1\n",
      "TCGA-22-5472-11A-11D-1633-05           1                        0\n",
      "TCGA-85-7844-01A-11D-2123-05           1                        1\n",
      "TCGA-22-5480-11A-01D-1633-05           1                        0\n",
      "TCGA-94-7557-01A-11D-2123-05           1                        1\n",
      "TCGA-56-8623-01A-11D-2398-05           1                        1\n",
      "TCGA-43-7656-01A-11D-2123-05           1                        1\n",
      "TCGA-68-8251-01A-11D-2294-05           1                        1\n",
      "TCGA-33-AASD-01A-11D-A408-05           1                        1\n",
      "TCGA-85-A5B5-01A-21D-A26N-05           1                        1\n",
      "TCGA-22-5473-01A-01D-1633-05           1                        1\n",
      "TCGA-94-8035-01A-11D-2245-05           1                        1\n",
      "TCGA-39-5016-01A-01D-1440-05           1                        1\n",
      "TCGA-39-5019-01A-01D-1818-05           1                        1\n",
      "TCGA-39-5029-01A-01D-1440-05           1                        1\n",
      "TCGA-77-A5GF-01A-21D-A27L-05           1                        1\n",
      "TCGA-XC-AA0X-01A-32D-A408-05           1                        1\n",
      "TCGA-22-5482-11A-01D-1633-05           1                        0\n",
      "TCGA-33-4566-11A-01D-1440-05           1                        0\n",
      "TCGA-68-7757-01B-11D-2294-05           1                        1\n",
      "TCGA-33-A5GW-01A-11D-A27L-05           1                        1\n",
      "TCGA-34-5231-01A-21D-1818-05           1                        1\n",
      "TCGA-39-5037-11A-01D-1440-05           1                        0\n",
      "TCGA-58-8391-01A-11D-2324-05           1                        1\n",
      "TCGA-33-4587-01A-11D-2123-05           1                        1\n",
      "TCGA-22-4599-01A-01D-1440-05           1                        1\n",
      "TCGA-79-5596-01A-31D-1947-05           1                        1\n",
      "TCGA-33-4566-01A-01D-1440-05           1                        1\n",
      "TCGA-77-7141-01A-11D-2043-05           1                        1\n",
      "TCGA-85-8666-01A-11D-2398-05           1                        1\n",
      "TCGA-O2-A52N-01A-11D-A26N-05           1                        1\n",
      "TCGA-22-5473-11A-11D-1633-05           1                        0\n",
      "TCGA-33-4589-01A-01D-1440-05           1                        1\n",
      "TCGA-39-5035-01A-01D-1440-05           1                        1\n",
      "TCGA-43-3394-11A-01D-1551-05           1                        0\n",
      "TCGA-68-A59J-01A-21D-A26N-05           1                        1\n",
      "TCGA-52-7812-01A-11D-2123-05           1                        1\n",
      "TCGA-90-6837-01A-11D-1947-05           1                        1\n",
      "TCGA-94-A5I4-01A-11D-A26N-05           1                        1\n",
      "TCGA-98-A538-01A-11D-A25R-05           1                        1\n",
      "TCGA-94-A5I6-01A-21D-A27L-05           1                        1\n",
      "TCGA-56-A5DS-01A-11D-A27L-05           1                        1\n",
      "TCGA-22-5471-11A-01D-1633-05           1                        0\n",
      "TCGA-60-2703-01A-11D-2043-05           1                        1\n",
      "TCGA-43-A475-01A-11D-A24I-05           1                        1\n",
      "TCGA-MF-A522-01A-11D-A25R-05           1                        1\n",
      "TCGA-39-5021-01A-01D-1440-05           1                        1\n",
      "TCGA-22-5478-11A-11D-1633-05           1                        0\n",
      "TCGA-LA-A446-01A-21D-A258-05           1                        1\n",
      "TCGA-94-A4VJ-01A-11D-A258-05           1                        1\n",
      "TCGA-46-6025-01A-11D-1818-05           1                        1\n",
      "TCGA-58-8386-01A-11D-2294-05           1                        1\n",
      "TCGA-85-8479-01A-11D-2324-05           1                        1\n",
      "TCGA-56-7579-01A-11D-2043-05           1                        1\n",
      "TCGA-22-5474-01A-01D-1633-05           1                        1\n",
      "TCGA-77-7463-01A-11D-2043-05           1                        1\n",
      "TCGA-22-5477-01A-01D-1633-05           1                        1\n",
      "TCGA-43-6647-01A-11D-1818-05           1                        1\n",
      "TCGA-22-4601-11A-01D-1440-05           1                        0\n",
      "TCGA-85-7710-01A-11D-2123-05           1                        1\n",
      "TCGA-63-7020-01A-11D-1947-05           1                        1\n",
      "TCGA-NC-A5HR-01A-21D-A26N-05           1                        1\n",
      "TCGA-33-4589-11A-01D-1440-05           1                        0\n",
      "TCGA-39-5028-01A-01D-1440-05           1                        1\n",
      "TCGA-90-A4ED-01A-31D-A258-05           1                        1\n",
      "TCGA-63-5131-01A-01D-1440-05           1                        1\n",
      "TCGA-77-6843-01A-11D-1947-05           1                        1\n",
      "TCGA-85-8354-01A-31D-2324-05           1                        1\n",
      "TCGA-33-4583-11A-01D-1440-05           1                        0\n",
      "TCGA-96-7544-01A-11D-2043-05           1                        1\n",
      "TCGA-56-8626-01A-11D-2398-05           1                        1\n",
      "TCGA-L3-A4E7-01A-11D-A258-05           1                        1\n",
      "TCGA-43-A56V-01A-11D-A26N-05           1                        1\n",
      "TCGA-77-6842-01A-11D-1947-05           1                        1\n",
      "TCGA-37-A5EM-01A-21D-A27L-05           1                        1\n",
      "TCGA-22-4601-01A-01D-1440-05           1                        1\n",
      "TCGA-22-5480-01A-01D-1633-05           1                        1\n",
      "TCGA-18-5592-11A-11D-1633-05           1                        0\n",
      "TCGA-NC-A5HT-01A-11D-A26N-05           1                        1\n",
      "TCGA-18-5592-01A-01D-1633-05           1                        1\n",
      "TCGA-90-7964-01A-21D-2185-05           1                        1\n",
      "TCGA-34-7107-01A-11D-1947-05           1                        1\n",
      "attr(,\"assign\")\n",
      "[1] 0 1\n",
      "attr(,\"contrasts\")\n",
      "attr(,\"contrasts\")$sample_type\n",
      "[1] \"contr.treatment\"\n",
      "\n",
      "An object of class \"MArrayLM\"\n",
      "$coefficients\n",
      "           (Intercept) sample_typePrimary Tumor\n",
      "cg00000029   -1.754809             -0.007678149\n",
      "cg00000165   -2.525597              2.687733713\n",
      "cg00000236    2.974544              0.065941613\n",
      "cg00000289    1.737575             -0.662286385\n",
      "cg00000321   -1.431677              1.046537149\n",
      "329299 more rows ...\n",
      "\n",
      "$rank\n",
      "[1] 2\n",
      "\n",
      "$assign\n",
      "[1] 0 1\n",
      "\n",
      "$qr\n",
      "$qr\n",
      "                              (Intercept) sample_typePrimary Tumor\n",
      "TCGA-33-6737-01A-11D-1818-05 -20.29778313             -18.22859165\n",
      "TCGA-O2-A52S-01A-11D-A25R-05   0.04926646              -6.14153455\n",
      "TCGA-77-7338-01A-11D-2043-05   0.04926646               0.01581938\n",
      "TCGA-98-A539-01A-31D-A25R-05   0.04926646               0.01581938\n",
      "TCGA-63-A5MT-01A-21D-A26N-05   0.04926646               0.01581938\n",
      "407 more rows ...\n",
      "\n",
      "$qraux\n",
      "[1] 1.049266 1.015819\n",
      "\n",
      "$pivot\n",
      "[1] 1 2\n",
      "\n",
      "$tol\n",
      "[1] 1e-07\n",
      "\n",
      "$rank\n",
      "[1] 2\n",
      "\n",
      "\n",
      "$df.residual\n",
      "[1] 410 410 410 410 410\n",
      "329299 more elements ...\n",
      "\n",
      "$sigma\n",
      "cg00000029 cg00000165 cg00000236 cg00000289 cg00000321 \n",
      " 0.9108121  1.3193328  0.4724132  0.5990647  1.0667723 \n",
      "329299 more elements ...\n",
      "\n",
      "$cov.coefficients\n",
      "                         (Intercept) sample_typePrimary Tumor\n",
      "(Intercept)               0.02380952              -0.02380952\n",
      "sample_typePrimary Tumor -0.02380952               0.02651223\n",
      "\n",
      "$stdev.unscaled\n",
      "           (Intercept) sample_typePrimary Tumor\n",
      "cg00000029   0.1543033                0.1628258\n",
      "cg00000165   0.1543033                0.1628258\n",
      "cg00000236   0.1543033                0.1628258\n",
      "cg00000289   0.1543033                0.1628258\n",
      "cg00000321   0.1543033                0.1628258\n",
      "329299 more rows ...\n",
      "\n",
      "$pivot\n",
      "[1] 1 2\n",
      "\n",
      "$Amean\n",
      "cg00000029 cg00000165 cg00000236 cg00000289 cg00000321 \n",
      "-1.7617046 -0.1118555  3.0337632  1.1428028 -0.4918257 \n",
      "329299 more elements ...\n",
      "\n",
      "$method\n",
      "[1] \"ls\"\n",
      "\n",
      "$design\n",
      "                             (Intercept) sample_typePrimary Tumor\n",
      "TCGA-33-6737-01A-11D-1818-05           1                        1\n",
      "TCGA-O2-A52S-01A-11D-A25R-05           1                        1\n",
      "TCGA-77-7338-01A-11D-2043-05           1                        1\n",
      "TCGA-98-A539-01A-31D-A25R-05           1                        1\n",
      "TCGA-63-A5MT-01A-21D-A26N-05           1                        1\n",
      "407 more rows ...\n",
      "\n"
     ]
    }
   ],
   "source": [
    "###########_____________ DMC analysis________________###########\n",
    "design <- model.matrix(~ sample_type, data = clinical)\n",
    "print(design)\n",
    "\n",
    "# fit the linear model \n",
    "fit <- lmFit(mval, design)\n",
    "fit2 <- eBayes(fit)\n",
    "\n",
    "print(fit)\n",
    "\n",
    "#extracting significantly methylated probes\n",
    "deff.meth = topTable(fit2, coef=ncol(design), sort.by=\"p\",number = nrow(mval), adjust.method = \"BY\")\n",
    "\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 25,
   "id": "8c1338e0",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:57:21.698072Z",
     "iopub.status.busy": "2022-12-28T13:57:21.696620Z",
     "iopub.status.idle": "2022-12-28T13:57:39.050915Z",
     "shell.execute_reply": "2022-12-28T13:57:39.049046Z"
    },
    "papermill": {
     "duration": 17.371964,
     "end_time": "2022-12-28T13:57:39.053099",
     "exception": false,
     "start_time": "2022-12-28T13:57:21.681135",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<dl>\n",
       "\t<dt>$cg02627286</dt>\n",
       "\t\t<dd>NULL</dd>\n",
       "\t<dt>$cg16732616</dt>\n",
       "\t\t<dd>NULL</dd>\n",
       "\t<dt>$cg06962177</dt>\n",
       "\t\t<dd>NULL</dd>\n",
       "\t<dt>$cg26521404</dt>\n",
       "\t\t<dd>NULL</dd>\n",
       "\t<dt>$cg25191628</dt>\n",
       "\t\t<dd>NULL</dd>\n",
       "\t<dt>$cg02443967</dt>\n",
       "\t\t<dd>NULL</dd>\n",
       "\t<dt>$cg22167515</dt>\n",
       "\t\t<dd>NULL</dd>\n",
       "\t<dt>$cg07078225</dt>\n",
       "\t\t<dd>NULL</dd>\n",
       "\t<dt>$cg16768018</dt>\n",
       "\t\t<dd>NULL</dd>\n",
       "\t<dt>$cg14823851</dt>\n",
       "\t\t<dd>NULL</dd>\n",
       "</dl>\n"
      ],
      "text/latex": [
       "\\begin{description}\n",
       "\\item[\\$cg02627286] NULL\n",
       "\\item[\\$cg16732616] NULL\n",
       "\\item[\\$cg06962177] NULL\n",
       "\\item[\\$cg26521404] NULL\n",
       "\\item[\\$cg25191628] NULL\n",
       "\\item[\\$cg02443967] NULL\n",
       "\\item[\\$cg22167515] NULL\n",
       "\\item[\\$cg07078225] NULL\n",
       "\\item[\\$cg16768018] NULL\n",
       "\\item[\\$cg14823851] NULL\n",
       "\\end{description}\n"
      ],
      "text/markdown": [
       "$cg02627286\n",
       ":   NULL\n",
       "$cg16732616\n",
       ":   NULL\n",
       "$cg06962177\n",
       ":   NULL\n",
       "$cg26521404\n",
       ":   NULL\n",
       "$cg25191628\n",
       ":   NULL\n",
       "$cg02443967\n",
       ":   NULL\n",
       "$cg22167515\n",
       ":   NULL\n",
       "$cg07078225\n",
       ":   NULL\n",
       "$cg16768018\n",
       ":   NULL\n",
       "$cg14823851\n",
       ":   NULL\n",
       "\n",
       "\n"
      ],
      "text/plain": [
       "$cg02627286\n",
       "NULL\n",
       "\n",
       "$cg16732616\n",
       "NULL\n",
       "\n",
       "$cg06962177\n",
       "NULL\n",
       "\n",
       "$cg26521404\n",
       "NULL\n",
       "\n",
       "$cg25191628\n",
       "NULL\n",
       "\n",
       "$cg02443967\n",
       "NULL\n",
       "\n",
       "$cg22167515\n",
       "NULL\n",
       "\n",
       "$cg07078225\n",
       "NULL\n",
       "\n",
       "$cg16768018\n",
       "NULL\n",
       "\n",
       "$cg14823851\n",
       "NULL\n"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "text/html": [
       "<strong>png:</strong> 2"
      ],
      "text/latex": [
       "\\textbf{png:} 2"
      ],
      "text/markdown": [
       "**png:** 2"
      ],
      "text/plain": [
       "png \n",
       "  2 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "# Visualization\n",
    "# plot the top 10 most significantly differentially methylated CpGs \n",
    "png(paste0(save_dir, \"top 10 most significantly differentially methylated CpGs.png\"), width=1000, height=750)\n",
    "par(mfrow=c(2,5))\n",
    "sapply(rownames(deff.meth)[1:10], function(cpg){\n",
    "  plotCpg(bval, cpg=cpg, pheno=clinical$sample_type, ylab = \"Beta values\")\n",
    "})\n",
    "dev.off()\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 26,
   "id": "5278514c",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:57:39.085359Z",
     "iopub.status.busy": "2022-12-28T13:57:39.083722Z",
     "iopub.status.idle": "2022-12-28T13:57:45.091463Z",
     "shell.execute_reply": "2022-12-28T13:57:45.088779Z"
    },
    "papermill": {
     "duration": 6.027303,
     "end_time": "2022-12-28T13:57:45.095005",
     "exception": false,
     "start_time": "2022-12-28T13:57:39.067702",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<strong>png:</strong> 2"
      ],
      "text/latex": [
       "\\textbf{png:} 2"
      ],
      "text/markdown": [
       "**png:** 2"
      ],
      "text/plain": [
       "png \n",
       "  2 "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "\n",
    "#making a volcano plot\n",
    "#making dataset\n",
    "dat <- data.frame(foldchange = fit[[\"coefficients\"]][,2], logPvalue =  -log10(fit2[[\"p.value\"]][,2]))\n",
    "dat$threshold <- as.factor(abs(dat$foldchange) < 0.4)\n",
    "\n",
    "#Visualization\n",
    "png(paste0(save_dir, \"volcano_plot_meth.png\"), width=1000, height=750)\n",
    "cols <- c(\"TRUE\" = \"grey\", \"FALSE\" = \"blue\")\n",
    "ggplot(data=dat, aes(x=foldchange, y = logPvalue, color=threshold)) +\n",
    "  geom_point(alpha=.6, size=1.2) +\n",
    "  scale_colour_manual(values = cols) +\n",
    "  geom_vline(xintercept = 0.4, colour=\"#990000\", linetype=\"dashed\") + \n",
    "  geom_vline(xintercept = - 0.4, colour=\"#990000\", linetype=\"dashed\") +\n",
    "  theme(legend.position=\"none\") +\n",
    "  xlab(\"Fold Change\") +\n",
    "  ylab(\"-log10 p value\") +\n",
    "  theme_bw() +\n",
    "  theme(legend.position = \"none\")\n",
    "dev.off()\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "38ca4899",
   "metadata": {
    "papermill": {
     "duration": 0.014003,
     "end_time": "2022-12-28T13:57:45.123046",
     "exception": false,
     "start_time": "2022-12-28T13:57:45.109043",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "------------ **DIRRERENTIALLY METHYLATED REGIONS** -----------------"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 27,
   "id": "254e92a8",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:57:45.154534Z",
     "iopub.status.busy": "2022-12-28T13:57:45.152813Z",
     "iopub.status.idle": "2022-12-28T13:58:26.133784Z",
     "shell.execute_reply": "2022-12-28T13:58:26.131808Z"
    },
    "papermill": {
     "duration": 41.012666,
     "end_time": "2022-12-28T13:58:26.149773",
     "exception": false,
     "start_time": "2022-12-28T13:57:45.137107",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Your contrast returned 211164 individually significant probes. We recommend the default setting of pcutoff in dmrcate().\n",
      "\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Formal class 'CpGannotated' [package \"DMRcate\"] with 1 slot\n",
      "  ..@ ranges:Formal class 'GRanges' [package \"GenomicRanges\"] with 7 slots\n",
      "  .. .. ..@ seqnames       :Formal class 'Rle' [package \"S4Vectors\"] with 4 slots\n",
      "  .. .. .. .. ..@ values         : Factor w/ 22 levels \"chr1\",\"chr2\",..: 1 2 3 4 5 6 7 8 9 10 ...\n",
      "  .. .. .. .. ..@ lengths        : int [1:22] 33105 24193 17544 13526 16618 24610 20190 14099 6878 16943 ...\n",
      "  .. .. .. .. ..@ elementMetadata: NULL\n",
      "  .. .. .. .. ..@ metadata       : list()\n",
      "  .. .. ..@ ranges         :Formal class 'IRanges' [package \"IRanges\"] with 6 slots\n",
      "  .. .. .. .. ..@ start          : int [1:329304] 15865 710097 714177 758829 763119 779995 805338 812539 812679 813317 ...\n",
      "  .. .. .. .. ..@ width          : int [1:329304] 1 1 1 1 1 1 1 1 1 1 ...\n",
      "  .. .. .. .. ..@ NAMES          : chr [1:329304] \"cg13869341\" \"cg15560884\" \"cg01014490\" \"cg11954957\" ...\n",
      "  .. .. .. .. ..@ elementType    : chr \"ANY\"\n",
      "  .. .. .. .. ..@ elementMetadata: NULL\n",
      "  .. .. .. .. ..@ metadata       : list()\n",
      "  .. .. ..@ strand         :Formal class 'Rle' [package \"S4Vectors\"] with 4 slots\n",
      "  .. .. .. .. ..@ values         : Factor w/ 3 levels \"+\",\"-\",\"*\": 3\n",
      "  .. .. .. .. ..@ lengths        : int 329304\n",
      "  .. .. .. .. ..@ elementMetadata: NULL\n",
      "  .. .. .. .. ..@ metadata       : list()\n",
      "  .. .. ..@ seqinfo        :Formal class 'Seqinfo' [package \"GenomeInfoDb\"] with 4 slots\n",
      "  .. .. .. .. ..@ seqnames   : chr [1:22] \"chr1\" \"chr2\" \"chr3\" \"chr4\" ...\n",
      "  .. .. .. .. ..@ seqlengths : int [1:22] NA NA NA NA NA NA NA NA NA NA ...\n",
      "  .. .. .. .. ..@ is_circular: logi [1:22] NA NA NA NA NA NA ...\n",
      "  .. .. .. .. ..@ genome     : chr [1:22] NA NA NA NA ...\n",
      "  .. .. ..@ elementMetadata:Formal class 'DFrame' [package \"S4Vectors\"] with 6 slots\n",
      "  .. .. .. .. ..@ rownames       : NULL\n",
      "  .. .. .. .. ..@ nrows          : int 329304\n",
      "  .. .. .. .. ..@ listData       :List of 4\n",
      "  .. .. .. .. .. ..$ stat   : num [1:329304] -3.79 -5.04 3.35 -7.4 7.3 ...\n",
      "  .. .. .. .. .. ..$ diff   : num [1:329304] -0.02958 -0.0574 0.00403 -0.09238 0.00771 ...\n",
      "  .. .. .. .. .. ..$ ind.fdr: num [1:329304] 2.89e-04 1.66e-06 1.34e-03 3.89e-12 7.66e-12 ...\n",
      "  .. .. .. .. .. ..$ is.sig : logi [1:329304] TRUE TRUE FALSE TRUE TRUE TRUE ...\n",
      "  .. .. .. .. ..@ elementType    : chr \"ANY\"\n",
      "  .. .. .. .. ..@ elementMetadata: NULL\n",
      "  .. .. .. .. ..@ metadata       : list()\n",
      "  .. .. ..@ elementType    : chr \"ANY\"\n",
      "  .. .. ..@ metadata       : list()\n"
     ]
    }
   ],
   "source": [
    "\n",
    "\n",
    "#######------------------------- Differentially methylated regions (DMR) analysis\n",
    "# setting some annotation\n",
    "myAnnotation <- cpg.annotate(object = mval,\n",
    "                             datatype = \"array\", \n",
    "                             what = \"M\", \n",
    "                             analysis.type = \"differential\", \n",
    "                             design = design, \n",
    "                             contrasts = FALSE, \n",
    "                             coef = \"sample_typePrimary Tumor\", \n",
    "                             arraytype = \"450K\",\n",
    "                             fdr = 0.001)\n",
    "## Your contrast returned 184682 individually significant probes. We recommend the default setting of pcutoff in dmrcate().\n",
    "\n",
    "\n",
    "str(myAnnotation)\n",
    "\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 28,
   "id": "55ec1e29",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:58:26.185317Z",
     "iopub.status.busy": "2022-12-28T13:58:26.183718Z",
     "iopub.status.idle": "2022-12-28T13:59:06.684819Z",
     "shell.execute_reply": "2022-12-28T13:59:06.683068Z"
    },
    "papermill": {
     "duration": 40.536395,
     "end_time": "2022-12-28T13:59:06.702719",
     "exception": false,
     "start_time": "2022-12-28T13:58:26.166324",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Fitting chr1...\n",
      "\n",
      "Fitting chr2...\n",
      "\n",
      "Fitting chr3...\n",
      "\n",
      "Fitting chr4...\n",
      "\n",
      "Fitting chr5...\n",
      "\n",
      "Fitting chr6...\n",
      "\n",
      "Fitting chr7...\n",
      "\n",
      "Fitting chr8...\n",
      "\n",
      "Fitting chr9...\n",
      "\n",
      "Fitting chr10...\n",
      "\n",
      "Fitting chr11...\n",
      "\n",
      "Fitting chr12...\n",
      "\n",
      "Fitting chr13...\n",
      "\n",
      "Fitting chr14...\n",
      "\n",
      "Fitting chr15...\n",
      "\n",
      "Fitting chr16...\n",
      "\n",
      "Fitting chr17...\n",
      "\n",
      "Fitting chr18...\n",
      "\n",
      "Fitting chr19...\n",
      "\n",
      "Fitting chr20...\n",
      "\n",
      "Fitting chr21...\n",
      "\n",
      "Fitting chr22...\n",
      "\n",
      "Demarcating regions...\n",
      "\n",
      "Done!\n",
      "\n",
      "using temporary cache /tmp/RtmpkkiFci/BiocFileCache\n",
      "\n",
      "snapshotDate(): 2020-10-27\n",
      "\n",
      "see ?DMRcatedata and browseVignettes('DMRcatedata') for documentation\n",
      "\n",
      "downloading 1 resources\n",
      "\n",
      "retrieving 1 resource\n",
      "\n",
      "loading from cache\n",
      "\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "GRanges object with 28555 ranges and 8 metadata columns:\n",
       "          seqnames              ranges strand |   no.cpgs min_smoothed_fdr\n",
       "             <Rle>           <IRanges>  <Rle> | <integer>        <numeric>\n",
       "      [1]     chr6   33156164-33181870      * |       258     1.71822e-306\n",
       "      [2]     chr6   33279563-33291947      * |       205     6.18577e-200\n",
       "      [3]     chr6   33128825-33155135      * |       157      0.00000e+00\n",
       "      [4]     chr6   32115964-32123701      * |       122     5.31102e-283\n",
       "      [5]     chr6   32935236-32942808      * |       118     1.43344e-245\n",
       "      ...      ...                 ...    ... .       ...              ...\n",
       "  [28551]    chr12   41086680-41086879      * |         3      2.46537e-14\n",
       "  [28552]     chr2 240499800-240499816      * |         2      1.87106e-13\n",
       "  [28553]     chr8   33457499-33457546      * |         2      1.33276e-13\n",
       "  [28554]    chr20   30327108-30327157      * |         3      2.20283e-14\n",
       "  [28555]     chr1   35250346-35250360      * |         2      1.45887e-18\n",
       "            Stouffer       HMFDR      Fisher     maxdiff     meandiff\n",
       "           <numeric>   <numeric>   <numeric>   <numeric>    <numeric>\n",
       "      [1]          0 1.32115e-31           0   0.2718494  -0.01179163\n",
       "      [2]          0 2.30206e-23           0  -0.2598073  -0.02406231\n",
       "      [3]          0 3.04743e-35           0  -0.2906048  -0.13607969\n",
       "      [4]          0 2.19483e-38           0   0.5123341   0.07420827\n",
       "      [5]          0 1.13963e-16           0  -0.0656017   0.00509796\n",
       "      ...        ...         ...         ...         ...          ...\n",
       "  [28551] 0.00388941 0.000448119 0.000971226  0.06247325  0.040056812\n",
       "  [28552] 0.00295381 0.007669314 0.003748445 -0.09465605 -0.067122118\n",
       "  [28553] 0.00462257 0.032406861 0.008420450  0.01617101  0.013162856\n",
       "  [28554] 0.34127611 0.056561793 0.090956543  0.00150483  0.000942769\n",
       "  [28555] 0.11147044 0.085197626 0.107892397 -0.01593387 -0.013335469\n",
       "               overlapping.genes\n",
       "                     <character>\n",
       "      [1] RNY4P10, SLC39A7, HS..\n",
       "      [2]    TAPBP, ZBTB22, DAXX\n",
       "      [3]                COL11A2\n",
       "      [4] PPT2, PPT2-EGFL8, PR..\n",
       "      [5] BRD2, BRD2-IT1, HLA-..\n",
       "      ...                    ...\n",
       "  [28551]                  CNTN1\n",
       "  [28552]                   <NA>\n",
       "  [28553]                 DUSP26\n",
       "  [28554]                   TPX2\n",
       "  [28555] GJB3, SMIM12, RP1-34..\n",
       "  -------\n",
       "  seqinfo: 22 sequences from an unspecified genome; no seqlengths"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "# DMR analysis\n",
    "DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)\n",
    "results.ranges <- extractRanges(DMRs)\n",
    "results.ranges\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 29,
   "id": "cd6c865e",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:59:06.737545Z",
     "iopub.status.busy": "2022-12-28T13:59:06.735954Z",
     "iopub.status.idle": "2022-12-28T13:59:06.764696Z",
     "shell.execute_reply": "2022-12-28T13:59:06.762938Z"
    },
    "papermill": {
     "duration": 0.048849,
     "end_time": "2022-12-28T13:59:06.767084",
     "exception": false,
     "start_time": "2022-12-28T13:59:06.718235",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "# visualization\n",
    "dmr.table <- data.frame(results.ranges)\n",
    "\n",
    "# setting up a variable for grouping and color\n",
    "\n",
    "pal <- brewer.pal(8,\"Dark2\")\n",
    "groups <- pal[1:length(unique(clinical$paper_Histologic.grade))]\n",
    "names(groups) <- levels(factor(clinical$paper_Histologic.grade))\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 30,
   "id": "d268728f",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:59:06.802728Z",
     "iopub.status.busy": "2022-12-28T13:59:06.801198Z",
     "iopub.status.idle": "2022-12-28T13:59:06.843493Z",
     "shell.execute_reply": "2022-12-28T13:59:06.841677Z"
    },
    "papermill": {
     "duration": 0.06247,
     "end_time": "2022-12-28T13:59:06.845891",
     "exception": false,
     "start_time": "2022-12-28T13:59:06.783421",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "#setting up the genomic region \n",
    "gen <- \"hg19\"\n",
    "# the index of the DMR that we will plot \n",
    "dmrIndex <- 2\n",
    "# coordinates are stored under results.ranges[dmrIndex]\n",
    "\n",
    "chrom <- as.character(seqnames(results.ranges[dmrIndex]))\n",
    "start <- as.numeric(start(results.ranges[dmrIndex]))\n",
    "end <- as.numeric(end(results.ranges[dmrIndex]))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 31,
   "id": "7e3cbb45",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T13:59:06.881686Z",
     "iopub.status.busy": "2022-12-28T13:59:06.880081Z",
     "iopub.status.idle": "2022-12-28T13:59:06.897885Z",
     "shell.execute_reply": "2022-12-28T13:59:06.896072Z"
    },
    "papermill": {
     "duration": 0.038178,
     "end_time": "2022-12-28T13:59:06.900385",
     "exception": false,
     "start_time": "2022-12-28T13:59:06.862207",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "# add 25% extra space to plot\n",
    "minbase <- start - (0.25*(end-start))\n",
    "maxbase <- end + (0.25*(end-start))\n",
    "\n",
    "\n",
    "# defining CpG islands track\n",
    "# download cpgislands for chromosome number 6 from ucsc"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "R",
   "language": "R",
   "name": "ir"
  },
  "language_info": {
   "codemirror_mode": "r",
   "file_extension": ".r",
   "mimetype": "text/x-r-source",
   "name": "R",
   "pygments_lexer": "r",
   "version": "4.0.5"
  },
  "papermill": {
   "default_parameters": {},
   "duration": 3929.146738,
   "end_time": "2022-12-28T13:59:08.247274",
   "environment_variables": {},
   "exception": null,
   "input_path": "__notebook__.ipynb",
   "output_path": "__notebook__.ipynb",
   "parameters": {},
   "start_time": "2022-12-28T12:53:39.100536",
   "version": "2.4.0"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
