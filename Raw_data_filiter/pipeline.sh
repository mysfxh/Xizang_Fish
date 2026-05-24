###鱼类数据返回后的第一次分析

###合并外群TOgz
bcftools merge   --threads 200   -m none   -Oz   -o ./All135_merged.vcf.gz   /home/xiongh/2026/Fish/fitered/new_1/tos/vcf/TO.snp.vcf.gz  ../All.filtered.SNP.vcf.gz


##数据过滤step1
QUAL_MIN=30
QD_MIN=2
FS_MAX=60
MQ_MIN=40
SOR_MAX=4
THREADS=16

SITE_KEEP_EXPR="QUAL>=${QUAL_MIN} && INFO/QD>=${QD_MIN} && INFO/FS<=${FS_MAX} && INFO/MQ>=${MQ_MIN} && INFO/SOR<=${SOR_MAX}"

nohup bash -c "
bcftools view --threads ${THREADS} \
  -i '${SITE_KEEP_EXPR}' \
  -Oz -o fish.01.sitefiltered.vcf.gz \
  ./All135_merged.vcf.gz && \
tabix -p vcf fish.01.sitefiltered.vcf.gz
" > fish.01.sitefiltered.log 2>&1 &

##数据过滤step2
nohup bash -c "vcftools \
--gzvcf fish.01.sitefiltered.vcf.gz \
--minDP 20 \
--max-missing 0.75 \
--maf 0.05 \
--min-alleles 2 \
--max-alleles 2 \
--recode \
--recode-INFO-all \
--stdout | bgzip -c > fish.02_0.75.10filtered.vcf.gz" > fish.02_0.75.10filtered.log 2>&1 & 

