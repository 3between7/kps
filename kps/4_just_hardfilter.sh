
#参数传入
#$indir
#$outdir
#$refpath
#$samplelist
#$gatk
#$gatkthreads
#$bcftools


ARGS=`getopt -a  -o i:o:r:l:g:t:b: -l indir:,outdir:,refpath:,samplelist:,gatk:,gatkthreads:,bcftools: -n 'example.sh' -- "$@"`

usage()
{
    echo "Usage: `basename $0` -indir indir -outdir outdir -refpath refpath -samplelist  samplelist -gatk path_of_gatk -gatkthreads gatk_threads -bcftools path_of_bcftools" 
    exit 1
}

[ $? != 0 ] &&  usage
[ $# -eq 0 ] && usage
eval set -- "${ARGS}"



while true
do
    case "$1" in
        -i|--indir)
            echo "#Option i, argument $2";
            indir=$2
			shift 2
            ;;
        -o|--outdir)
            echo "#Option o, argument $2";
			outdir=$2
            shift 2
            ;;
        -r|--refpath)
            echo "#Option r, argument $2";
            refpath=$2
			shift 2
            ;;
		-l|--samplelist)
            echo "#Option l, argument $2";
            samplelist=$2
			shift 2
            ;;
		-g|--gatk)
            echo "#Option g, argument $2";
            gatk=$2
			shift 2
            ;;
		-t|--gatkthreads)
            echo "#Option t, argument $2";
            gatkthreads=$2
			shift 2
            ;;
		-b|--bcftools)
            echo "#Option a, argument $2";
            bcftools=$2
			shift 2
            ;;
		--)
            shift
            break
            ;;
        *)
            echo "Internal error!"
            exit 1
            ;;
    esac
done





#生成目录
mkdir -p $outdir

baseDirForScriptSelf=$(cd "$(dirname "$0")"; pwd)

if [ ! -s "$samplelist" ];then
  rm -f $samplelist
fi

if [ ! -f "$samplelist" ]; then
  gvcffiles=`ls $indir | grep "g.vcf" |grep -v "idx"|awk -v xxx=$indir '{print "--variant",xxx"/"$0}'|sed 'H;$!d;g;s/\n/ /g'`
  gvcfnum=`ls $indir | grep "g.vcf" |grep -v "idx"|wc -l`
else
  gvcffiles=`less $samplelist|awk -v xxx=$indir '{print "--variant", xxx"/"$0".raw_variants.g.vcf"}'|sed 'H;$!d;g;s/\n/ /g'`
  gvcfnum=`wc -l $samplelist|cut -f 1`
fi
echo $gvcffiles

#再生成samplelist
if [ ! -f "$samplelist" ]; then
  touch "$samplelist"
  for file in `ls $indir/*.raw_variants.g.vcf`
  do
      echo $file
      sample=$(basename $file .raw_variants.g.vcf)
      echo $sample >> "$samplelist"
  done
fi




#===========================================合并成vcf===========================#
echo "#Joint samples and call SNP###########"
java -Xmx15g -XX:ParallelGCThreads=1 -jar $gatk/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $refpath $gvcffiles  -o $outdir/HaplotypeCaller.raw.vcf --disable_auto_index_creation_and_locking_when_reading_rods 1>$outdir/../sh/1_GenotypeGVCFs.error 2>$outdir/../sh/1_GenotypeGVCFs.log
#===========================================严格过滤===========================#
#所有参数不根据项目改
#QD  QualByDepth 变异位点可信度除以非参考样品的未过滤的深度
#FS  FisherStrand 使用Fisher’s精确检验来检测strand bias而得到的Fhred格式的p值。该值越小越好。使用Fisher检验的p value 评估序列水平的偏好性。 如变异位点只在正链或负链上。偏好性暗示着 假阳性。
#MQ  RMSMappingQuality  所有样品的reads的比对质量的平方根。
#SOR This is another way to estimate strand bias using a test similar to the symmetric odds ratio test.SOR is a related annotation that applies a different statistical test (using the SB counts) that is better for high coverage data.

#MQRankSum  对比对质量的 等级秩和检验。比较有参考碱基的reads和 有alt碱基的reads.这个只能应用于杂合位点，需要有ref的reads和alt的reads。
#ReadPosRankSum 等级秩和检验  有Alt的reads的末尾的距离？如果Alt仅仅出现在reads的末尾，则说明可能有错误。只能应用于杂合位点。

hardfilter_snpExpression1="QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0"
hardfilter_snpExpression2="MQRankSum < -12.5 || ReadPosRankSum < -8.0"
hardfilter_indelExpression1="QD < 2.0 || FS > 200.0 || SOR > 10.0"
hardfilter_indelExpression2="ReadPosRankSum < -20.0"


#InbreedingCoeff   The InbreedingCoeff statistic is a population-level calculation that is only available with 10 or more samples. If you have fewer samples you will need to omit that particular filter statement.
#if [ $gvcfnum -gt 10 ]  #不做InbreedingCoeff < -0.8这个，因为如果是家系群体，则会导致位点过滤错误。
#then 
#	hardfilter_indelExpression2="ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8"
#else
#	hardfilter_indelExpression2="ReadPosRankSum < -20.0"
#fi


echo "#select snp&indel#"
java -Xmx15g -XX:ParallelGCThreads=10 -jar  $gatk/GenomeAnalysisTK.jar -T SelectVariants -R $refpath -V $outdir/HaplotypeCaller.raw.vcf -selectType SNP -o $outdir/HaplotypeCaller.rawSNP.vcf --disable_auto_index_creation_and_locking_when_reading_rods 1>$outdir/../sh/2_SelectSNP.error 2>$outdir/../sh/2_SelectSNP.log

java -Xmx15g -XX:ParallelGCThreads=10 -jar  $gatk/GenomeAnalysisTK.jar -T SelectVariants -R $refpath -V $outdir/HaplotypeCaller.raw.vcf -selectType INDEL -o $outdir/HaplotypeCaller.rawINDEL.vcf --disable_auto_index_creation_and_locking_when_reading_rods 1>$outdir/../sh/3_SelectINDEL.error 2>$outdir/../sh/3_SelectINDEL.log

echo "#filter snp&indel#"
java -Xmx15g -XX:ParallelGCThreads=10 -jar  $gatk/GenomeAnalysisTK.jar -T VariantFiltration -R $refpath -V $outdir/HaplotypeCaller.rawSNP.vcf --filterExpression "$hardfilter_snpExpression1" --filterName "filter1" --filterExpression "$hardfilter_snpExpression2" --filterName "filter2" -o $outdir/HaplotypeCaller.filterSNP.vcf --disable_auto_index_creation_and_locking_when_reading_rods  1>$outdir/../sh/4_FilterSNP.error 2>$outdir/../sh/4_FilterSNP.log


java -Xmx15g -XX:ParallelGCThreads=10 -jar  $gatk/GenomeAnalysisTK.jar -T VariantFiltration -R $refpath -V $outdir/HaplotypeCaller.rawINDEL.vcf --filterExpression "$hardfilter_indelExpression1" --filterName "filter1" --filterExpression "$hardfilter_indelExpression2" --filterName "filter2" -o $outdir/HaplotypeCaller.filterINDEL.vcf --disable_auto_index_creation_and_locking_when_reading_rods 1>$outdir/../sh/5_FilterINDEL.error 2>$outdir/../sh/5_FilterINDEL.log


less $outdir/HaplotypeCaller.filterSNP.vcf|awk -F '\t' '{if(($7 == "PASS")||($1 ~ /#/)) {print $0}}' > $outdir/HaplotypeCaller.filterSNP.ks.vcf

less $outdir/HaplotypeCaller.filterINDEL.vcf|awk -F '\t' '{if(($7 == "PASS")||($1 ~ /#/)) {print $0}}' > $outdir/HaplotypeCaller.filterINDEL.ks.vcf


#将这两个合并起来
java -Xmx15g -XX:ParallelGCThreads=10 -jar  $gatk/GenomeAnalysisTK.jar -T CombineVariants  -R $refpath --variant:indel $outdir/HaplotypeCaller.filterINDEL.ks.vcf --variant:snp $outdir/HaplotypeCaller.filterSNP.ks.vcf -o $outdir/HaplotypeCaller.hardfilter.vcf -genotypeMergeOptions PRIORITIZE -priority indel,snp --disable_auto_index_creation_and_locking_when_reading_rods  1>$outdir/../sh/6_CombineVariants.error 2>$outdir/../sh/6_CombineVariants.log





#===========================================严格过滤===========================#
#snp cluster过滤 （5bp内如果有2个snp则过滤掉） --clusterWindowSize 5 --clusterSize 2
#snp cluster过滤不做了，此处注释掉：2018-12-05-jyl
#echo "#filter snp1"
#java -Xmx15g -XX:ParallelGCThreads=10 -jar  $gatk/GenomeAnalysisTK.jar -T VariantFiltration -R $refpath -V $outdir/HaplotypeCaller.hardfilter.vcf --clusterWindowSize 5 --clusterSize 2 -o $outdir/vcf1.vcf --disable_auto_index_creation_and_locking_when_reading_rods 1>$outdir/../sh/7_VariantFiltration1.error 2>$outdir/../sh/7_VariantFiltration1.log

#less $outdir/vcf1.vcf |awk -F '\t' '{if(($7 !~ "SnpCluster")||($1 ~ /#/)) {print $0}}' > $outdir/vcf2.vcf


echo "#filter snp2" #Indel附近SNP过滤（indel附近5bp内的SNP过滤掉）相邻INDEL过滤（两个indel距离小于10bp过滤掉）"
#$bcftools  filter -g 5 -G 10 -o $outdir/vcf3.vcf -O v $outdir/vcf2.vcf 
#因为不做snp cluster的过滤，此处修改了输入文件，由vcf1改为原始vcf文件HaplotypeCaller.hardfilter.vcf：2018-12-05-jyl
$bcftools  filter -g 5 -G 10 -o $outdir/vcf3.vcf -O v $outdir/HaplotypeCaller.hardfilter.vcf


#SNP的质量值(MQ)不低于20，将质量值Q20，即1/100,测序错误率大于1%的snp过滤掉。
#SNP的reads支持数不低于4总覆盖深度，而不是每个样品的覆盖深度！  
#For filtering purposes it is better to use QD than either QUAL or DP directly.
#echo "#filter snp3" 
#hardfilter_Expression1="MQ < 20.0"
#hardfilter_Expression2="DP < 4.0"
#java -Xmx10g -jar  $gatk/GenomeAnalysisTK.jar -T VariantFiltration -R $refpath -V $outdir/vcf3.vcf --filterExpression "$hardfilter_Expression1" --filterName "filter1" --filterExpression "$hardfilter_Expression2" --filterName "filter2" -o $outdir/vcf4.vcf --disable_auto_index_creation_and_locking_when_reading_rods 1>$outdir/../sh/VariantFiltration2.error 2>$outdir/../sh/VariantFiltration2.log
#less $outdir/vcf4.vcf|awk -F '\t' '{if((\$7 !~ \"filter\")||(\$1 ~ /\#/)) {print \$0}}' > $outdir/vcf5.vcf 
#这段不要，如果上面直接接着hardfilter，hardfiler已经过滤了，而不是接着VQSR。



#每个样品的GQ不低于20;如果低于20，则分型就改变成缺失 ./.
#GQ 基因型的质量值(Genotype Quality)。Phred格式(Phred_scaled)的质量值，表示在该位点该基因型存在的可能性；该值越高，则Genotype的可能性越大；计算方法：Phred值 = -10 * log (1-p) p为基因型存在的概率。
echo "#Gfilter snp4"
Gfilter="GQ < 20.0"
#java -Xmx15g -XX:ParallelGCThreads=10 -jar $gatk/GenomeAnalysisTK.jar -T VariantFiltration -R $refpath -V $outdir/vcf3.vcf -G_filter "$Gfilter" -G_filterName "lowGQ" -o $outdir/vcf.final.vcf --setFilteredGtToNocall --disable_auto_index_creation_and_locking_when_reading_rods 1>$outdir/../sh/8_VariantFiltration3.error 2>$outdir/../sh/8_VariantFiltration3.log
#对于低质量的位点(标记为lowGQ，满足GQ < 20.0)，只标记lowGQ，不转换成缺失(./.)，即保留原始的基因型分型：2018-12-05-jyl
java -Xmx15g -XX:ParallelGCThreads=10 -jar $gatk/GenomeAnalysisTK.jar -T VariantFiltration -R $refpath -V $outdir/vcf3.vcf -G_filter "$Gfilter" -G_filterName "lowGQ" -o $outdir/vcf.final.vcf --disable_auto_index_creation_and_locking_when_reading_rods 1>$outdir/../sh/8_VariantFiltration3.error 2>$outdir/../sh/8_VariantFiltration3.log


#最后选择结果 #分别生成snp和indel的结果
echo "split snp indel"
java -Xmx15g -XX:ParallelGCThreads=10 -jar  $gatk/GenomeAnalysisTK.jar -T SelectVariants -R $refpath -V $outdir/vcf.final.vcf -selectType SNP -o $outdir/vcf.finalSNP.vcf --disable_auto_index_creation_and_locking_when_reading_rods 1>$outdir/../sh/9_SelectVariantsSNP.error 2>$outdir/../sh/9_SelectVariantsSNP.log

java -Xmx15g -XX:ParallelGCThreads=10 -jar  $gatk/GenomeAnalysisTK.jar -T SelectVariants -R $refpath -V $outdir/vcf.final.vcf -selectType INDEL -o $outdir/vcf.finalINDEL.vcf --disable_auto_index_creation_and_locking_when_reading_rods 1>$outdir/../sh/10_SelectVariantsINDEL.error 2>$outdir/../sh/10_SelectVariantsINDEL.log

#创建check文件用于表示程序结束
touch $outdir/../sh/hardfilter.sh.check
