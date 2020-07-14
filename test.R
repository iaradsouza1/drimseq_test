library(DRIMSeq)

# Load data
load("drimseq_objects.rda")

# ann = samples df
# genotypes = genotypes df for chr 1 (only SNPs in biallelic sites)
# g_rng = GRange object for gene coordinates in chr 1
# rng = GRange object for SNPs coordinates
# cts = counts df, in scaledTPM 

# Build dmSQTL object 
w <- 5000
d <- dmSQTLdata(counts = cts, gene_ranges = g_rng,
                genotypes = genotypes, snp_ranges = rng,
                samples = ann, window = w)

# Filter transcritps
# Get the minimal number of samples where genes should be expressed. I'm using 70% of group samples.
min_samples <- floor(length(ann$sample_id) * 0.7)
d <- dmFilter(d, min_samps_gene_expr = min_samples, min_samps_feature_expr = 5,
              minor_allele_freq = 5, min_gene_expr = 10, min_feature_expr = 10)

d <- d[1:10,]
# Estimate precision
d <- dmPrecision(d)

# Fit model
d <- dmFit(d)

# Test
d <- dmTest(d)

sessionInfo() 


  