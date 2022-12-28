{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "974fa625",
   "metadata": {
    "_execution_state": "idle",
    "_uuid": "051d70d956493feee0c6d64651c6a088724dca2a",
    "execution": {
     "iopub.execute_input": "2022-12-28T00:21:14.555981Z",
     "iopub.status.busy": "2022-12-28T00:21:14.552498Z",
     "iopub.status.idle": "2022-12-28T01:20:16.658168Z",
     "shell.execute_reply": "2022-12-28T01:20:16.655758Z"
    },
    "papermill": {
     "duration": 3542.121506,
     "end_time": "2022-12-28T01:20:16.663776",
     "exception": false,
     "start_time": "2022-12-28T00:21:14.542270",
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
      "The following object is masked from ‘package:limma’:\n",
      "\n",
      "    plotMA\n",
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
      "Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19\n",
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
      "Loading required package: IlluminaHumanMethylation450kmanifest\n",
      "\n",
      "Loading required package: grid\n",
      "\n",
      "Loading required package: sesameData\n",
      "\n",
      "Loading required package: ExperimentHub\n",
      "\n",
      "Loading required package: AnnotationHub\n",
      "\n",
      "Loading required package: BiocFileCache\n",
      "\n",
      "Loading required package: dbplyr\n",
      "\n",
      "\n",
      "Attaching package: 'AnnotationHub'\n",
      "\n",
      "\n",
      "The following object is masked from 'package:Biobase':\n",
      "\n",
      "    cache\n",
      "\n",
      "\n",
      "Loading required package: rmarkdown\n",
      "\n",
      "Loading sesameData.\n",
      "\n",
      "\n",
      "----------------------------------------------------------\n",
      "| SEnsible Step-wise Analysis of DNA MEthylation (SeSAMe)\n",
      "| --------------------------------------------------------\n",
      "| Please cache the annotation data for your array platform\n",
      "| (e.g. EPIC) by calling \"sesameDataCache(\"EPIC\")\"\n",
      "| or \"sesameDataCacheAll()\". This needs to be done only\n",
      "| once per SeSAMe installation.\n",
      "----------------------------------------------------------\n",
      "\n",
      "\n",
      "Warning message:\n",
      "\"replacing previous import 'minfi::getMeth' by 'bsseq::getMeth' when loading 'DMRcate'\"\n",
      "── \u001b[1mAttaching packages\u001b[22m ─────────────────────────────────────── tidyverse 1.3.2 ──\n",
      "\u001b[32m✔\u001b[39m \u001b[34mggplot2\u001b[39m 3.3.6      \u001b[32m✔\u001b[39m \u001b[34mpurrr  \u001b[39m 0.3.5 \n",
      "\u001b[32m✔\u001b[39m \u001b[34mtibble \u001b[39m 3.1.8      \u001b[32m✔\u001b[39m \u001b[34mdplyr  \u001b[39m 1.0.10\n",
      "\u001b[32m✔\u001b[39m \u001b[34mtidyr  \u001b[39m 1.2.1      \u001b[32m✔\u001b[39m \u001b[34mforcats\u001b[39m 0.5.2 \n",
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
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mident()\u001b[39m      masks \u001b[34mdbplyr\u001b[39m::ident()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mlag()\u001b[39m        masks \u001b[34mstats\u001b[39m::lag()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mpurrr\u001b[39m::\u001b[32mnone()\u001b[39m       masks \u001b[34mlocfit\u001b[39m::none()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mggplot2\u001b[39m::\u001b[32mPosition()\u001b[39m masks \u001b[34mBiocGenerics\u001b[39m::Position(), \u001b[34mbase\u001b[39m::Position()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mpurrr\u001b[39m::\u001b[32mreduce()\u001b[39m     masks \u001b[34mGenomicRanges\u001b[39m::reduce(), \u001b[34mIRanges\u001b[39m::reduce()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mrename()\u001b[39m     masks \u001b[34mS4Vectors\u001b[39m::rename()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32mslice()\u001b[39m      masks \u001b[34mXVector\u001b[39m::slice(), \u001b[34mIRanges\u001b[39m::slice()\n",
      "\u001b[31m✖\u001b[39m \u001b[34mdplyr\u001b[39m::\u001b[32msql()\u001b[39m        masks \u001b[34mdbplyr\u001b[39m::sql()\n",
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
    "library(\"limma\")\n",
    "library(\"minfi\")\n",
    "library(\"RColorBrewer\")\n",
    "library(\"missMethyl\") # Can take a short time...\n",
    "library(\"minfiData\")\n",
    "library(\"Gviz\")\n",
    "library(\"sesame\") #for analysis\n",
    "library(\"DMRcate\")\n",
    "library(\"DMRcatedata\")\n",
    "library(\"stringr\")\n",
    "library(\"readr\")\n",
    "library(\"tidyverse\")\n",
    "library(\"edgeR\")\n",
    "library(\"IlluminaHumanMethylation450kanno.ilmn12.hg19\")\n",
    "library(\"IlluminaHumanMethylation450kmanifest\")\n",
    "library(\"stringr\")\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "592151ce",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:20:17.039066Z",
     "iopub.status.busy": "2022-12-28T01:20:17.001182Z",
     "iopub.status.idle": "2022-12-28T01:20:17.052138Z",
     "shell.execute_reply": "2022-12-28T01:20:17.050040Z"
    },
    "papermill": {
     "duration": 0.382177,
     "end_time": "2022-12-28T01:20:17.055646",
     "exception": false,
     "start_time": "2022-12-28T01:20:16.673469",
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
   "id": "f7fc813b",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:20:17.082628Z",
     "iopub.status.busy": "2022-12-28T01:20:17.080651Z",
     "iopub.status.idle": "2022-12-28T01:20:17.095117Z",
     "shell.execute_reply": "2022-12-28T01:20:17.093059Z"
    },
    "papermill": {
     "duration": 0.03121,
     "end_time": "2022-12-28T01:20:17.097658",
     "exception": false,
     "start_time": "2022-12-28T01:20:17.066448",
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
   "id": "1535f325",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:20:17.119823Z",
     "iopub.status.busy": "2022-12-28T01:20:17.117860Z",
     "iopub.status.idle": "2022-12-28T01:20:17.132287Z",
     "shell.execute_reply": "2022-12-28T01:20:17.130245Z"
    },
    "papermill": {
     "duration": 0.028214,
     "end_time": "2022-12-28T01:20:17.134931",
     "exception": false,
     "start_time": "2022-12-28T01:20:17.106717",
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
   "id": "600d6ba0",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:20:17.157151Z",
     "iopub.status.busy": "2022-12-28T01:20:17.155211Z",
     "iopub.status.idle": "2022-12-28T01:21:17.076590Z",
     "shell.execute_reply": "2022-12-28T01:21:17.074773Z"
    },
    "papermill": {
     "duration": 59.935471,
     "end_time": "2022-12-28T01:21:17.079464",
     "exception": false,
     "start_time": "2022-12-28T01:20:17.143993",
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
   "id": "29dd6aa4",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:21:17.103000Z",
     "iopub.status.busy": "2022-12-28T01:21:17.101438Z",
     "iopub.status.idle": "2022-12-28T01:21:19.183554Z",
     "shell.execute_reply": "2022-12-28T01:21:19.181865Z"
    },
    "papermill": {
     "duration": 2.096169,
     "end_time": "2022-12-28T01:21:19.185950",
     "exception": false,
     "start_time": "2022-12-28T01:21:17.089781",
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
   "execution_count": 7,
   "id": "9153752d",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:21:19.207664Z",
     "iopub.status.busy": "2022-12-28T01:21:19.206245Z",
     "iopub.status.idle": "2022-12-28T01:21:22.460382Z",
     "shell.execute_reply": "2022-12-28T01:21:22.458567Z"
    },
    "papermill": {
     "duration": 3.26743,
     "end_time": "2022-12-28T01:21:22.462621",
     "exception": false,
     "start_time": "2022-12-28T01:21:19.195191",
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
    "# chose those has no NA values in rows\n",
    "probe <- probe.na[probe.na == 0]\n",
    "met <- met[row.names(met) %in% names(probe), ]\n",
    "\n",
    "## remove probes that match chromosomes X and Y \n",
    "keep <- !(row.names(met) %in% ann450k$Name[ann450k$chr %in% c(\"chrX\",\"chrY\")])\n",
    "table(keep)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "id": "3c9b9cdb",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:21:22.485432Z",
     "iopub.status.busy": "2022-12-28T01:21:22.483803Z",
     "iopub.status.idle": "2022-12-28T01:21:29.093335Z",
     "shell.execute_reply": "2022-12-28T01:21:29.091077Z"
    },
    "papermill": {
     "duration": 6.62547,
     "end_time": "2022-12-28T01:21:29.097486",
     "exception": false,
     "start_time": "2022-12-28T01:21:22.472016",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "met <- met[keep, ]\n",
    "rm(keep) # remove no further needed probes."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "bcfc003e",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:21:29.121957Z",
     "iopub.status.busy": "2022-12-28T01:21:29.120002Z",
     "iopub.status.idle": "2022-12-28T01:21:29.272596Z",
     "shell.execute_reply": "2022-12-28T01:21:29.270316Z"
    },
    "papermill": {
     "duration": 0.169448,
     "end_time": "2022-12-28T01:21:29.276566",
     "exception": false,
     "start_time": "2022-12-28T01:21:29.107118",
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
    "no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "id": "a768947d",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:21:29.302597Z",
     "iopub.status.busy": "2022-12-28T01:21:29.300673Z",
     "iopub.status.idle": "2022-12-28T01:21:32.240916Z",
     "shell.execute_reply": "2022-12-28T01:21:32.238726Z"
    },
    "papermill": {
     "duration": 2.956826,
     "end_time": "2022-12-28T01:21:32.243484",
     "exception": false,
     "start_time": "2022-12-28T01:21:29.286658",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]\n",
    "#snps with maf <= 0.05\n",
    "snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]\n",
    "\n",
    "# filter met\n",
    "met <- met[row.names(met) %in% c(no.snp.probe, snp5.probe), ]\n",
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
   "execution_count": 11,
   "id": "74313541",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:21:32.267917Z",
     "iopub.status.busy": "2022-12-28T01:21:32.266180Z",
     "iopub.status.idle": "2022-12-28T01:22:12.654607Z",
     "shell.execute_reply": "2022-12-28T01:22:12.652530Z"
    },
    "papermill": {
     "duration": 40.403379,
     "end_time": "2022-12-28T01:22:12.657217",
     "exception": false,
     "start_time": "2022-12-28T01:21:32.253838",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "# filtre met\n",
    "met <- met[ -which(row.names(met) %in% crs.reac), ]\n",
    "bval <- met\n",
    "\n",
    "## converting beta values to m_values\n",
    "## m = log2(beta/1-beta)\n",
    "mval <- t(apply(met, 1, function(x) log2(x/(1-x))))\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "id": "3e96a124",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:22:12.681306Z",
     "iopub.status.busy": "2022-12-28T01:22:12.679651Z",
     "iopub.status.idle": "2022-12-28T01:22:18.606313Z",
     "shell.execute_reply": "2022-12-28T01:22:18.604372Z"
    },
    "papermill": {
     "duration": 5.941187,
     "end_time": "2022-12-28T01:22:18.608938",
     "exception": false,
     "start_time": "2022-12-28T01:22:12.667751",
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
   "execution_count": null,
   "id": "55078f91",
   "metadata": {
    "papermill": {
     "duration": 0.057555,
     "end_time": "2022-12-28T01:22:20.334127",
     "exception": false,
     "start_time": "2022-12-28T01:22:20.276572",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "id": "7d3922a9",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:22:20.454257Z",
     "iopub.status.busy": "2022-12-28T01:22:20.452077Z",
     "iopub.status.idle": "2022-12-28T01:22:32.845010Z",
     "shell.execute_reply": "2022-12-28T01:22:32.841728Z"
    },
    "papermill": {
     "duration": 12.458254,
     "end_time": "2022-12-28T01:22:32.849418",
     "exception": false,
     "start_time": "2022-12-28T01:22:20.391164",
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
    "dev.off()\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "id": "4f7033c1",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:22:32.879101Z",
     "iopub.status.busy": "2022-12-28T01:22:32.875259Z",
     "iopub.status.idle": "2022-12-28T01:23:05.343784Z",
     "shell.execute_reply": "2022-12-28T01:23:05.341629Z"
    },
    "papermill": {
     "duration": 32.495989,
     "end_time": "2022-12-28T01:23:05.357543",
     "exception": false,
     "start_time": "2022-12-28T01:22:32.861554",
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
    "#density\n",
    "png(paste0(save_dir, \"Densities_after_GDCprepare.png\"), width=1200, height=850)\n",
    "plotDensities(met, legend=FALSE, main= \"Densities sesame data\")\n",
    "dev.off() \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "dbd9e76d",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:23:05.383631Z",
     "iopub.status.busy": "2022-12-28T01:23:05.381458Z",
     "iopub.status.idle": "2022-12-28T01:23:09.703420Z",
     "shell.execute_reply": "2022-12-28T01:23:09.701430Z"
    },
    "papermill": {
     "duration": 4.338038,
     "end_time": "2022-12-28T01:23:09.705926",
     "exception": false,
     "start_time": "2022-12-28T01:23:05.367888",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "#############plot mean b values\n",
    "# remove probes with NA (similar to na.omit)\n",
    "tcga_na_filtered <- data.met[rowSums(is.na(assay(data.met))) == 0,]\n",
    "\n",
    "\n",
    "df <- data.frame(\n",
    "  \"Sample.mean\" = colMeans(assay(tcga_na_filtered), na.rm = TRUE),\n",
    "  \"groups\" = tcga_na_filtered$sample_type\n",
    ")\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "id": "602f5fcb",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:23:09.729648Z",
     "iopub.status.busy": "2022-12-28T01:23:09.728142Z",
     "iopub.status.idle": "2022-12-28T01:23:13.957721Z",
     "shell.execute_reply": "2022-12-28T01:23:13.955330Z"
    },
    "papermill": {
     "duration": 4.24496,
     "end_time": "2022-12-28T01:23:13.961050",
     "exception": false,
     "start_time": "2022-12-28T01:23:09.716090",
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
   "cell_type": "code",
   "execution_count": 17,
   "id": "96ff7112",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:23:13.988219Z",
     "iopub.status.busy": "2022-12-28T01:23:13.986075Z",
     "iopub.status.idle": "2022-12-28T01:23:37.060264Z",
     "shell.execute_reply": "2022-12-28T01:23:37.058410Z"
    },
    "papermill": {
     "duration": 23.09052,
     "end_time": "2022-12-28T01:23:37.062818",
     "exception": false,
     "start_time": "2022-12-28T01:23:13.972298",
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
   "execution_count": 18,
   "id": "06776a26",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:23:37.098060Z",
     "iopub.status.busy": "2022-12-28T01:23:37.096509Z",
     "iopub.status.idle": "2022-12-28T01:24:00.203243Z",
     "shell.execute_reply": "2022-12-28T01:24:00.201417Z"
    },
    "papermill": {
     "duration": 23.132098,
     "end_time": "2022-12-28T01:24:00.206302",
     "exception": false,
     "start_time": "2022-12-28T01:23:37.074204",
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
    "###########################################################\n",
    "\n",
    "\n",
    "#densities mval\n",
    "png(paste0(save_dir, \"densities_mval.png\"), width=1200, height=850)\n",
    "plotDensities(mval, legend=FALSE, main= \"Densities Mval data\")\n",
    "dev.off()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "dc29d071",
   "metadata": {
    "papermill": {
     "duration": 0.016177,
     "end_time": "2022-12-28T01:24:00.238563",
     "exception": false,
     "start_time": "2022-12-28T01:24:00.222386",
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
   "execution_count": 19,
   "id": "761e49cc",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:24:00.270879Z",
     "iopub.status.busy": "2022-12-28T01:24:00.269237Z",
     "iopub.status.idle": "2022-12-28T01:24:00.290822Z",
     "shell.execute_reply": "2022-12-28T01:24:00.288932Z"
    },
    "papermill": {
     "duration": 0.038426,
     "end_time": "2022-12-28T01:24:00.293424",
     "exception": false,
     "start_time": "2022-12-28T01:24:00.254998",
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
   "execution_count": 20,
   "id": "62419a3a",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:24:00.324723Z",
     "iopub.status.busy": "2022-12-28T01:24:00.322779Z",
     "iopub.status.idle": "2022-12-28T01:24:00.744138Z",
     "shell.execute_reply": "2022-12-28T01:24:00.742305Z"
    },
    "papermill": {
     "duration": 0.438506,
     "end_time": "2022-12-28T01:24:00.747208",
     "exception": false,
     "start_time": "2022-12-28T01:24:00.308702",
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
   "execution_count": 21,
   "id": "b48b51f5",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:24:00.774209Z",
     "iopub.status.busy": "2022-12-28T01:24:00.772536Z",
     "iopub.status.idle": "2022-12-28T01:24:18.278952Z",
     "shell.execute_reply": "2022-12-28T01:24:18.276426Z"
    },
    "papermill": {
     "duration": 17.522707,
     "end_time": "2022-12-28T01:24:18.282023",
     "exception": false,
     "start_time": "2022-12-28T01:24:00.759316",
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
    "\n",
    "###########_____________ DMC analysis________________###########\n",
    "design <- model.matrix(~ sample_type, data = clinical)\n",
    "# fit the linear model \n",
    "fit <- lmFit(mval, design)\n",
    "fit2 <- eBayes(fit)\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "id": "8c1de36e",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:24:18.313441Z",
     "iopub.status.busy": "2022-12-28T01:24:18.310803Z",
     "iopub.status.idle": "2022-12-28T01:24:35.551093Z",
     "shell.execute_reply": "2022-12-28T01:24:35.549367Z"
    },
    "papermill": {
     "duration": 17.259083,
     "end_time": "2022-12-28T01:24:35.553845",
     "exception": false,
     "start_time": "2022-12-28T01:24:18.294762",
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
    "\n",
    "#extracting significantly methylated probes\n",
    "deff.meth = topTable(fit2, coef=ncol(design), sort.by=\"p\",number = nrow(mval), adjust.method = \"BY\")\n",
    "\n",
    "\n",
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
   "execution_count": 23,
   "id": "e7e16f20",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:24:35.583202Z",
     "iopub.status.busy": "2022-12-28T01:24:35.581694Z",
     "iopub.status.idle": "2022-12-28T01:24:41.604910Z",
     "shell.execute_reply": "2022-12-28T01:24:41.602648Z"
    },
    "papermill": {
     "duration": 6.041034,
     "end_time": "2022-12-28T01:24:41.607780",
     "exception": false,
     "start_time": "2022-12-28T01:24:35.566746",
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
   "cell_type": "code",
   "execution_count": 24,
   "id": "e6d9e4af",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2022-12-28T01:24:41.637597Z",
     "iopub.status.busy": "2022-12-28T01:24:41.635588Z",
     "iopub.status.idle": "2022-12-28T01:25:25.378766Z",
     "shell.execute_reply": "2022-12-28T01:25:25.377043Z"
    },
    "papermill": {
     "duration": 43.772984,
     "end_time": "2022-12-28T01:25:25.393392",
     "exception": false,
     "start_time": "2022-12-28T01:24:41.620408",
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
   "duration": 3856.674485,
   "end_time": "2022-12-28T01:25:27.134091",
   "environment_variables": {},
   "exception": null,
   "input_path": "__notebook__.ipynb",
   "output_path": "__notebook__.ipynb",
   "parameters": {},
   "start_time": "2022-12-28T00:21:10.459606",
   "version": "2.4.0"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
