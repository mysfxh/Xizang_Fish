#!/bin/bash
set -euo pipefail

############################################################
# 计算 π、Ho、ROH、FROH
# 注意：输入 VCF 必须是基础过滤后的高质量 SNP 数据
# 不要使用 LD pruning 后的数据
############################################################

########################
# 1. 用户需要修改的参数
########################

# 高质量过滤后的 VCF，不要用 LD pruning 后的 VCF
VCF="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/03.variants/Fish_unout/fish.without_removed_samples.vcf.gz"

# 群体分组文件：第一列样本名，第二列群体名
POPMAP="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/03.variants/Fish_unout/popname.txt"

# 参考基因组 fai 文件，用于计算 FROH 分母
# 如果没有 fai，先运行：samtools faidx reference.fa
REF_FAI="/media/server/Elements/BN250801NJ01S02N1-安徽大学1例高原鳅全基因组图谱开发/03.Hi-C/HiC.review.assembly.chr.fa.fai"

# 输出目录
OUTDIR="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/03.variants/Fish_unout//diversity_ROH_results"

# π 滑窗大小和步长
WINDOW_SIZE=100000
WINDOW_STEP=10000

# PLINK ROH 参数
ROH_KB=100
ROH_SNP=50
ROH_DENSITY=50
ROH_GAP=1000
ROH_WINDOW_SNP=50
ROH_WINDOW_HET=1
ROH_WINDOW_MISSING=5
ROH_WINDOW_THRESHOLD=0.05

########################
# 2. 创建输出目录
########################

mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/pop_list
mkdir -p ${OUTDIR}/pi_site
mkdir -p ${OUTDIR}/pi_window
mkdir -p ${OUTDIR}/heterozygosity
mkdir -p ${OUTDIR}/plink_roh
mkdir -p ${OUTDIR}/summary

########################
# 3. 检查输入文件
########################

echo "检查输入文件..."

if [ ! -f "$VCF" ]; then
    echo "错误：找不到 VCF 文件：$VCF"
    exit 1
fi

if [ ! -f "${VCF}.tbi" ] && [ ! -f "${VCF}.csi" ]; then
    echo "VCF 没有索引，正在建立索引..."
    bcftools index -t "$VCF"
fi

if [ ! -f "$POPMAP" ]; then
    echo "错误：找不到 popmap 文件：$POPMAP"
    exit 1
fi

if [ ! -f "$REF_FAI" ]; then
    echo "错误：找不到参考基因组 fai 文件：$REF_FAI"
    echo "请先运行：samtools faidx reference.fa"
    exit 1
fi

########################
# 4. 检查样本名是否匹配
########################

echo "检查 VCF 和 popmap 样本名是否一致..."

bcftools query -l "$VCF" | sort > ${OUTDIR}/summary/vcf_samples.txt
awk '{print $1}' "$POPMAP" | sort > ${OUTDIR}/summary/popmap_samples.txt

comm -23 ${OUTDIR}/summary/popmap_samples.txt ${OUTDIR}/summary/vcf_samples.txt > ${OUTDIR}/summary/popmap_not_in_vcf.txt
comm -13 ${OUTDIR}/summary/popmap_samples.txt ${OUTDIR}/summary/vcf_samples.txt > ${OUTDIR}/summary/vcf_not_in_popmap.txt

if [ -s ${OUTDIR}/summary/popmap_not_in_vcf.txt ]; then
    echo "警告：以下 popmap 样本不在 VCF 中："
    cat ${OUTDIR}/summary/popmap_not_in_vcf.txt
    echo "请检查样本名。"
fi

if [ -s ${OUTDIR}/summary/vcf_not_in_popmap.txt ]; then
    echo "警告：以下 VCF 样本不在 popmap 中："
    cat ${OUTDIR}/summary/vcf_not_in_popmap.txt
fi

########################
# 5. 按群体生成样本列表
########################

echo "生成每个群体的样本列表..."

awk '{print $2}' "$POPMAP" | sort | uniq | while read pop
do
    awk -v p="$pop" '$2==p {print $1}' "$POPMAP" > ${OUTDIR}/pop_list/${pop}.txt
done

echo "群体列表："
ls ${OUTDIR}/pop_list/*.txt

########################
# 6. 计算每个位点的 π
########################

echo "开始计算 site π..."

for popfile in ${OUTDIR}/pop_list/*.txt
do
    pop=$(basename "$popfile" .txt)

    echo "计算群体 ${pop} 的 site π..."

    vcftools \
        --gzvcf "$VCF" \
        --keep "$popfile" \
        --site-pi \
        --out ${OUTDIR}/pi_site/${pop}
done

########################
# 7. 汇总每个群体平均 π
########################

echo "汇总 mean π..."

echo -e "Pop\tMean_pi\tN_sites" > ${OUTDIR}/summary/pop_mean_pi.txt

for file in ${OUTDIR}/pi_site/*.sites.pi
do
    pop=$(basename "$file" .sites.pi)

    awk -v p="$pop" '
    BEGIN{OFS="\t"}
    NR>1 && $3!="nan" && $3!="NA" {
        sum += $3
        n++
    }
    END{
        if(n>0) print p, sum/n, n
        else print p, "NA", 0
    }' "$file" >> ${OUTDIR}/summary/pop_mean_pi.txt
done

########################
# 8. 计算滑窗 π
########################

echo "开始计算 window π..."

for popfile in ${OUTDIR}/pop_list/*.txt
do
    pop=$(basename "$popfile" .txt)

    echo "计算群体 ${pop} 的 window π..."

    vcftools \
        --gzvcf "$VCF" \
        --keep "$popfile" \
        --window-pi ${WINDOW_SIZE} \
        --window-pi-step ${WINDOW_STEP} \
        --out ${OUTDIR}/pi_window/${pop}
done

########################
# 9. 计算 Ho
########################

echo "开始计算 Ho..."

vcftools \
    --gzvcf "$VCF" \
    --het \
    --out ${OUTDIR}/heterozygosity/all_samples

# vcftools --het 输出通常为：
# INDV O(HOM) E(HOM) N_SITES F
# Ho = 1 - O(HOM) / N_SITES

awk '
BEGIN{
    OFS="\t";
    print "Sample","O_HOM","E_HOM","N_SITES","F","Ho"
}
NR>1{
    Ho = 1 - $2/$4
    print $1,$2,$3,$4,$5,Ho
}' ${OUTDIR}/heterozygosity/all_samples.het > ${OUTDIR}/heterozygosity/sample_Ho.txt

# 加入群体信息
awk '
BEGIN{OFS="\t"}
NR==FNR{
    pop[$1]=$2
    next
}
FNR==1{
    print $0,"Pop"
    next
}
{
    print $0,pop[$1]
}' "$POPMAP" ${OUTDIR}/heterozygosity/sample_Ho.txt > ${OUTDIR}/heterozygosity/sample_Ho_with_pop.txt

# 汇总每个群体 Ho
awk '
BEGIN{
    OFS="\t";
    print "Pop","Mean_Ho","SD_Ho","N_individuals"
}
NR>1{
    p=$7
    x=$6
    sum[p]+=x
    sumsq[p]+=x*x
    n[p]++
}
END{
    for(p in sum){
        mean=sum[p]/n[p]
        if(n[p]>1){
            sd=sqrt((sumsq[p]-sum[p]^2/n[p])/(n[p]-1))
        }else{
            sd="NA"
        }
        print p,mean,sd,n[p]
    }
}' ${OUTDIR}/heterozygosity/sample_Ho_with_pop.txt > ${OUTDIR}/summary/pop_mean_Ho.txt

########################
# 10. VCF 转 PLINK
########################

echo "VCF 转换为 PLINK 格式..."

plink \
    --vcf "$VCF" \
    --double-id \
    --allow-extra-chr \
    --make-bed \
    --out ${OUTDIR}/plink_roh/all.filtered

########################
# 11. 计算 ROH
########################

echo "开始计算 ROH..."

plink \
    --bfile ${OUTDIR}/plink_roh/all.filtered \
    --allow-extra-chr \
    --homozyg \
    --homozyg-kb ${ROH_KB} \
    --homozyg-snp ${ROH_SNP} \
    --homozyg-density ${ROH_DENSITY} \
    --homozyg-gap ${ROH_GAP} \
    --homozyg-window-snp ${ROH_WINDOW_SNP} \
    --homozyg-window-het ${ROH_WINDOW_HET} \
    --homozyg-window-missing ${ROH_WINDOW_MISSING} \
    --homozyg-window-threshold ${ROH_WINDOW_THRESHOLD} \
    --out ${OUTDIR}/plink_roh/all.ROH

########################
# 12. 计算 FROH
########################

echo "计算 FROH..."

# 用参考基因组 fai 计算基因组长度
# 如果你只分析染色体，可以提前把 REF_FAI 改成只包含染色体的 fai
GENOME_SIZE=$(awk '{sum+=$2} END{print sum}' "$REF_FAI")

echo "用于 FROH 的基因组长度为：${GENOME_SIZE} bp" > ${OUTDIR}/summary/genome_size_used_for_FROH.txt

# PLINK .hom.indiv 常见列：
# FID IID PHE NSEG KB KBAVG
# FROH = 总 ROH 长度 / 基因组长度
# KB 需要乘以 1000 转换为 bp

awk -v G=$GENOME_SIZE '
BEGIN{
    OFS="\t";
    print "Sample","N_ROH","Total_ROH_kb","Mean_ROH_kb","FROH"
}
NR>1{
    sample=$2
    nroh=$4
    total_kb=$5
    mean_kb=$6
    froh=(total_kb*1000)/G
    print sample,nroh,total_kb,mean_kb,froh
}' ${OUTDIR}/plink_roh/all.ROH.hom.indiv > ${OUTDIR}/plink_roh/sample_FROH.txt

# 加入群体信息
awk '
BEGIN{OFS="\t"}
NR==FNR{
    pop[$1]=$2
    next
}
FNR==1{
    print $0,"Pop"
    next
}
{
    print $0,pop[$1]
}' "$POPMAP" ${OUTDIR}/plink_roh/sample_FROH.txt > ${OUTDIR}/plink_roh/sample_FROH_with_pop.txt

# 汇总每个群体 FROH
awk '
BEGIN{
    OFS="\t";
    print "Pop","Mean_FROH","SD_FROH","Mean_N_ROH","Mean_Total_ROH_kb","N_individuals"
}
NR>1{
    p=$6
    froh=$5
    nroh=$2
    total=$3

    sum_froh[p]+=froh
    sumsq_froh[p]+=froh*froh
    sum_nroh[p]+=nroh
    sum_total[p]+=total
    n[p]++
}
END{
    for(p in sum_froh){
        mean=sum_froh[p]/n[p]
        if(n[p]>1){
            sd=sqrt((sumsq_froh[p]-sum_froh[p]^2/n[p])/(n[p]-1))
        }else{
            sd="NA"
        }
        print p,mean,sd,sum_nroh[p]/n[p],sum_total[p]/n[p],n[p]
    }
}' ${OUTDIR}/plink_roh/sample_FROH_with_pop.txt > ${OUTDIR}/summary/pop_mean_FROH.txt

########################
# 13. ROH 长度分级
########################

echo "统计不同长度等级 ROH..."

# PLINK .hom 常见列：
# FID IID PHE CHR SNP1 SNP2 POS1 POS2 KB NSNP DENSITY PHOM PHET
# 第 2 列 IID 是样本名，第 9 列 KB 是 ROH 长度

awk '
BEGIN{
    OFS="\t";
    print "Sample","CHR","POS1","POS2","ROH_kb","NSNP","Class"
}
NR>1{
    sample=$2
    chr=$4
    pos1=$7
    pos2=$8
    roh_kb=$9
    nsnp=$10

    if(roh_kb >= 100 && roh_kb < 500) class="100-500kb"
    else if(roh_kb >= 500 && roh_kb < 1000) class="500kb-1Mb"
    else if(roh_kb >= 1000 && roh_kb < 5000) class="1-5Mb"
    else if(roh_kb >= 5000) class=">5Mb"
    else class="<100kb"

    print sample,chr,pos1,pos2,roh_kb,nsnp,class
}' ${OUTDIR}/plink_roh/all.ROH.hom > ${OUTDIR}/plink_roh/ROH_length_class.txt

# 每个个体不同长度等级 ROH 统计
awk '
BEGIN{
    OFS="\t";
    print "Sample","Class","Total_ROH_kb","N_ROH"
}
NR>1{
    key=$1"\t"$7
    sum[key]+=$5
    n[key]++
}
END{
    for(k in sum){
        print k,sum[k],n[k]
    }
}' ${OUTDIR}/plink_roh/ROH_length_class.txt > ${OUTDIR}/plink_roh/sample_ROH_class_summary.txt

# 加群体信息
awk '
BEGIN{OFS="\t"}
NR==FNR{
    pop[$1]=$2
    next
}
FNR==1{
    print $0,"Pop"
    next
}
{
    print $0,pop[$1]
}' "$POPMAP" ${OUTDIR}/plink_roh/sample_ROH_class_summary.txt > ${OUTDIR}/plink_roh/sample_ROH_class_summary_with_pop.txt

# 每个群体不同长度等级 ROH 统计
awk '
BEGIN{
    OFS="\t";
    print "Pop","Class","Mean_Total_ROH_kb","Mean_N_ROH","N_individuals_with_this_class"
}
NR>1{
    sample=$1
    class=$2
    total=$3
    nroh=$4
    pop=$5

    key=pop"\t"class
    sum_total[key]+=total
    sum_nroh[key]+=nroh
    n[key]++
}
END{
    for(k in sum_total){
        print k,sum_total[k]/n[k],sum_nroh[k]/n[k],n[k]
    }
}' ${OUTDIR}/plink_roh/sample_ROH_class_summary_with_pop.txt > ${OUTDIR}/summary/pop_ROH_class_summary.txt

########################
# 14. 合并 π、Ho、FROH 群体结果
########################

echo "合并群体层面的 π、Ho、FROH..."

# 整理 π
awk 'NR==1{next}{pi[$1]=$2; pisites[$1]=$3} END{for(p in pi) print p,pi[p],pisites[p]}' OFS="\t" ${OUTDIR}/summary/pop_mean_pi.txt > ${OUTDIR}/summary/tmp_pi.txt

# 整理 Ho
awk 'NR==1{next}{ho[$1]=$2; hosd[$1]=$3; hon[$1]=$4} END{for(p in ho) print p,ho[p],hosd[p],hon[p]}' OFS="\t" ${OUTDIR}/summary/pop_mean_Ho.txt > ${OUTDIR}/summary/tmp_ho.txt

# 整理 FROH
awk 'NR==1{next}{froh[$1]=$2; frohsd[$1]=$3; nroh[$1]=$4; total[$1]=$5; fn[$1]=$6} END{for(p in froh) print p,froh[p],frohsd[p],nroh[p],total[p],fn[p]}' OFS="\t" ${OUTDIR}/summary/pop_mean_FROH.txt > ${OUTDIR}/summary/tmp_froh.txt

# 用 awk 合并
awk '
BEGIN{
    OFS="\t";
    print "Pop","Mean_pi","Pi_N_sites","Mean_Ho","SD_Ho","N_individuals_Ho","Mean_FROH","SD_FROH","Mean_N_ROH","Mean_Total_ROH_kb","N_individuals_FROH"
}
NR==FNR{
    pi[$1]=$2
    pisites[$1]=$3
    pops[$1]=1
    next
}
FILENAME==ARGV[2]{
    ho[$1]=$2
    hosd[$1]=$3
    hon[$1]=$4
    pops[$1]=1
    next
}
FILENAME==ARGV[3]{
    froh[$1]=$2
    frohsd[$1]=$3
    nroh[$1]=$4
    total[$1]=$5
    fn[$1]=$6
    pops[$1]=1
    next
}
END{
    for(p in pops){
        print p,pi[p],pisites[p],ho[p],hosd[p],hon[p],froh[p],frohsd[p],nroh[p],total[p],fn[p]
    }
}' ${OUTDIR}/summary/tmp_pi.txt ${OUTDIR}/summary/tmp_ho.txt ${OUTDIR}/summary/tmp_froh.txt > ${OUTDIR}/summary/pop_diversity_ROH_summary.txt

rm -f ${OUTDIR}/summary/tmp_pi.txt ${OUTDIR}/summary/tmp_ho.txt ${OUTDIR}/summary/tmp_froh.txt

########################
# 15. 完成
########################

echo "全部分析完成！"
echo "主要结果文件："
echo "1. 群体平均 π：${OUTDIR}/summary/pop_mean_pi.txt"
echo "2. 群体平均 Ho：${OUTDIR}/summary/pop_mean_Ho.txt"
echo "3. 个体 FROH：${OUTDIR}/plink_roh/sample_FROH_with_pop.txt"
echo "4. 群体平均 FROH：${OUTDIR}/summary/pop_mean_FROH.txt"
echo "5. ROH 长度分级：${OUTDIR}/summary/pop_ROH_class_summary.txt"
echo "6. 汇总表：${OUTDIR}/summary/pop_diversity_ROH_summary.txt"
