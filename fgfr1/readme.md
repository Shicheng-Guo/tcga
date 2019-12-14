sx We downloaded 11,093 gene expression quantification data derived from RNA-seq data in TCGA database on February 24th, 2019. The RNA-seq data covered 32 cancer types.  However 9 cancer types were excluded since the low samples size for control samples (N<=1). Log2-transformed fragments per kilobase of transcript per million mapped reads upper quartile (FPKM-UQ) derived from HTSeq49, is applied for differential gene expression analysis. Bayesian generalized linear model (bayesglm) from ARM package (v1.10-1) was applied for differential gene expression analysis. metafor package (v2.1-0) was applied for meta-analysis cross 23 cancer types. Cox proportional hazards regression model was applied for survival analysis to TCGA overall survival times (R survival package v0.9). 
x
x
x
x
s
x
x
x
x
x
x
