# BD TWAS Workflow
This directory contains the steps taken to generate TWAS results for both the amygdala and sACC samples using the bipolar disorder GWAS made available by the Psychiatric Genomics Consortium (PGC) and published by Stahl, E.A. et al [[1](#references)].

## 1) Filter SNPs

```
qsub filter_snps/filter_snps_amygdala_gene.sh
# and/or
qsub filter_snps/filter_snps_sacc_gene.sh
```

The first step of this workflow outputs PLINK `.bed`, `.fam`, `.bim` and `.bed` files with Lieber Institute of Brain Development (LIBD) DNA genotype data for the samples used in this RNA-seq project.

## 2) Build gene-level PLINK .bed files
```
sh run_build_bims_amygdala_gene.sh
# and/or
sh run_build_bims_sacc_gene.sh
```

This script builds the `.bim` and `.bed` files that contain gene windows with 500kb of coverage on each side of the gene. The outputs are in separate gene window subdirectories in `{subregion}_gene/bim_files/{subregion}_gene_{1:n gene windows}`.

## 3) Compute TWAS Weights individually
```
sh run_compute_weights_indv_amygdala_full_gene.sh
# and/or
sh run_compute_weights_indv_sacc_full_gene.sh
```

Here, we make use of Gusev et al's `FUSION.compute_weights.R` script [[2](#references)], which produces expression weights one gene at a time, taking in the above `build_bims.R` output. The flags are customizable and allow the user to specify multiple different models for testing. We set the heritability score P-value threshold to above `1` so heritability is not filtered. The output for this script can be found in `{subregion}_gene/out_files/gene_{1:n gene windows}`.

## 4) Merge Individual TWAS gene weights
```
qsub compute_weights_amygdala_gene.sh
# and/or
qsub compute_weights_sacc_gene.sh
```

The results from step #3 are merged together by generating an index `pos_info.Rdata`. `FUSION.profile_wgt.R` from Gusev et al.'s TWAS workflow is also called to produce a per-gene profile with summary of the data and model performance in `{subregion}_gene.profile.err`.

## 5) Convert hg19 GWAS to hg38
```
qsub run_process-hg19-gwas.sh
```

The GWAS file sourced from PGC is converted from an hg19 to an hg38 coordinate map and saved in the working directory as `PGC_BIP_hg38_clean.txt`.

## 6) TWAS results into a single Rdata file
```
qsub run_read_twas_amygdala.sh
# and/or
qsub run_read_twas_sacc.sh
```

All of the TWAS results from both the sACC and amygdala subregions are combined into a single Rdata file that is used for downstream analysis.

## 7) TWAS plot generation
```
qsub run_generate_twas_plots.sh
```

This script generates a number of different plots and outputs their corresponding tables into `analysis/plots/` and `analysis/tables/`, respectively.

## 8) Enrichment tests on TWAS FDR < 5% genes
```
qsub run_enrichment_test.sh
```

Executes the enrichment test that ran with `{clusterProfiler}` [[3](#references)]. Results found in `analysis/tables/BD_EnrichmentTest.xlsx`.

## References
1. Stahl, E.A., Breen, G., Forstner, A.J. et al. Genome-wide association study identifies 30 loci associated with bipolar disorder. Nat Genet 51, 793-803 (2019). https://doi.org/10.1038/s41588-019-0397-8
2. Gusev, A., Ko, A., Shi, H. et al. Integrative approaches for large-scale transcriptome-wide association studies. Nat Genet 48, 245–252 (2016). https://doi.org/10.1038/ng.3506
3. Yu G, Wang L, Han Y, He Q (2012). “clusterProfiler: an R package for comparing biological themes among gene clusters.” OMICS: A Journal of Integrative Biology, 16(5), 284-287. doi: 10.1089/omi.2011.0118.
