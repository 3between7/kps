
#合并整理SNP和INDEL的mergeMapFrq文件
sort -n -t $'\t' -k 1n,1 -k 2n,2 -k 3n,3 /lustre/project/R2D/Rice60k-CAU/3-MergeSnpIndelInfo/AllSnpIndel_sort.txt > /lustre/project/R2D/Rice60k-CAU/3-MergeSnpIndelInfo/AllSnpIndel_sort.txt

#过滤MAF<0.01的位点
awk 'BEGIN{FS="\t";OFS="\t"}{if ($6>=0.01){print $0}}' /lustre/project/R2D/Rice60k-CAU/3-MergeSnpIndelInfo/AllSnpIndel_sort.txt > /lustre/project/R2D/Rice60k-CAU/4-FiltrateMAF/AllSnpIndel_sort_filteMAF.txt

#按染色体分割文件
int=1
while (($int <= 12)) 
do
	filename="C$int_sort_filtMAF.txt"
	awk -v i=$int 'BEGIN{FS="\t";OFS="\t"}{if ($1==i){print $0}}' /lustre/project/R2D/Rice60k-CAU/4-FiltrateMAF/AllSnpIndel_sort_filteMAF.txt > $filename
	let "int++"
