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

mkdir plink_result

plink --vcf ../fish.02_0.75.10filtered.vcf.gz \
  --allow-extra-chr \
  --double-id \
 --make-bed \
  --out fish.02_0.75.10filtered
# 1. 重新设置 SNP ID，避免重复 ID
plink2 \
  --bfile fish.02_0.75.10filtered\
  --set-all-var-ids @:#\$r,\$a \
  --new-id-max-allele-len 1000 truncate \
  --make-bed \
  --out fish.02_0.75.10filtered.unique


# 2. 进行 LD pruning，生成 .prune.in 和 .prune.out
plink \
  --bfile fish.02_0.75.10filtered.unique \
  --allow-extra-chr \
  --indep-pairwise 100 1 0.8 \
  --out fish.02_0.75.10filtered.unique.unique.100_1_0.8


# 3. 提取 LD 过滤后保留下来的 SNP
plink \
  --bfile fish.02_0.75.10filtered.unique \
  --allow-extra-chr \
  --extract fish.02_0.75.10filtered.unique.unique.100_1_0.8.prune.in \
  --make-bed \
  --out fish.02_0.75.10filtered.snapp.filtered.unique.100_1_0.8


# 4. 把 LD 过滤后的数据导出为 VCF
plink \
  --bfile fish.02_0.75.10filtered.snapp.filtered.unique.100_1_0.8 \
  --allow-extra-chr \
  --recode vcf-iid bgz \
