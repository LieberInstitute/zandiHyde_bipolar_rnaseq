
for i in $(seq 1 22); do
	fastQTL.static \
	--vcf /users/schadinh/lieber/vcf_data/sacc_raggr-only_20190308.vcf.gz \
	--bed /users/schadinh/lieber/qqnorm/sacc_qqnorm_20190220.txt.gz \
	--window 5e5 \
	--region ${i}:0-258956422 \
	--out /users/schadinh/lieber/fastqtl_output/bipseq_20190308/sacc_chr${i}_fastqtl_20190308.out \
	--cov /users/schadinh/lieber/covariates/sacc_sex-bp-5geno-pc_covs_20190306.txt;
done
