ChrLength=(158337067 137060424 121430405 120829699 121191424 119458736 112638659 113384836 105708250 104305016 107310763 91163125 84240350 84648390 85296676 81724687 75158596 66004023 64057457 72042655 71599096 61435874 52530062 62714930 42904170 51681464 45407902 46312546 51505224 148823899)
ChrPrefix="gnl|UMD3.1|GK0000"
C=1
N=10
while (($C<=30))
do
	if [ $C -lt $N ]
	then
		chromosome=$ChrPrefix"0$C.2"
		echo "$chromosome"
	else
		chromosome=$ChrPrefix"$C.2"
		echo "$chromosome"
	fi
	workpath="./BfCattleMAFSplitFile"
	filename="$workpath/BfCattle_SNP.maf.Chr$C"
	awk -v chrom=$chromosome -v C=$C 'BEGIN{FS="\t";OFS="\t"}{if ($1==chrom){print C,$2,$3,$4,$5}}' BfCattle.SNP30+.maf > $filename
	echo "the head 3 lines of $filename is "
	head -n 3 $filename
	let index=$C-1 
	echo "the length of chrom$C is ${ChrLength[$index]}"
#	echo "$C	0       NA      NA      NA" > $filename | cat $filename.temp >> $filename 
#	echo "$C	${ChrLength[$index]}	NA	NA      NA" >> $filename
#	sed -i 's/$chromosome/$C/g' $filename
#	echo "the head 3 lines of $filename" 
#	head -n 3 $filename
#	echo "the end 3 lines of $filename"
#	tail -n 3 $filename
	let "C++"
#	rm $filename.temp
done

cat BfCattle_SNP.maf.Chr1 BfCattle_SNP.maf.Chr2 BfCattle_SNP.maf.Chr3 BfCattle_SNP.maf.Chr4 BfCattle_SNP.maf.Chr5 BfCattle_SNP.maf.Chr6 BfCattle_SNP.maf.Chr7 BfCattle_SNP.maf.Chr8 BfCattle_SNP.maf.Chr9 BfCattle_SNP.maf.Chr10 BfCattle_SNP.maf.Chr11 BfCattle_SNP.maf.Chr12 BfCattle_SNP.maf.Chr13 BfCattle_SNP.maf.Chr14 BfCattle_SNP.maf.Chr15 BfCattle_SNP.maf.Chr16 BfCattle_SNP.maf.Chr17 BfCattle_SNP.maf.Chr18 BfCattle_SNP.maf.Chr19 BfCattle_SNP.maf.Chr20 BfCattle_SNP.maf.Chr21 BfCattle_SNP.maf.Chr22 BfCattle_SNP.maf.Chr23 BfCattle_SNP.maf.Chr24 BfCattle_SNP.maf.Chr25 BfCattle_SNP.maf.Chr26 BfCattle_SNP.maf.Chr27 BfCattle_SNP.maf.Chr28 BfCattle_SNP.maf.Chr29 BfCattle_SNP.maf.Chr30 >> BfCattle_SNP.maf

