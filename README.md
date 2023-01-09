# SE_analysis
Combined analysis of super-enhancer data (ROSE), enhancer quantification, gene body expression and TF motifs

## What to integrate:

- super-enhancer (SE) output from ROSE
- quantification of regular enhancers
- gene body expression (quantified by HOMER)
- TF motif analysis (HOMER)


## How-to/Steps:

- make constituent enhancer dataset from all AL timepoints
- run ROSE on all timepoints
- quantify regular enhancers for all timepoints
- per SE, back-reference regular enhancers and quantification data
- per SE, find nearest genes (e.g. within 50 kbp)
- per SE, do TF motif analysis with HOMER

