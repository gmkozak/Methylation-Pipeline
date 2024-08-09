# Methylation-Pipeline
Code for pipelines used to analyze enzymatic methylation sequencing (EM-seq) data from Lepidopteran insects.

Pipeline A uses Trimmomatic, bwa-meth, and MethylDackel.

Pipeline B uses Trim Galore and Biscuit.

Differential methylation calling from bedGraph files (Pipeline A) or cov files (Pipeline B) is performed in R using code in the MethylSig.All.Final.R
