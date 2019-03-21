ChrLength=(158337067 137060425 121430404 120829700 121191425 119458734 112638660 113384548 105708240 104304789 107310761 91163126 84240322 84648390 85296493 81724688 75158593 66004024 64057359 72040631 71599041 61435577 52529545 62714859 42904055 51681083 45407519 46312547 51505225 148823900)
C=1
while (($C<=30))
do
	filename="Cattle_UnionSites.Chr$C"
	awk -v chrom=$C 'BEGIN{FS="\t";OFS="\t"}{if ($2==chrom){print $0}}' /lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/20181124/4-filate35bp/Catlle_UnionSites.addpos2 > $filename.temp
	let index=$C-1 
	echo "the length of chrom$C is ${ChrLength[$index]}"
	echo "NA	$C	0       0	NA      NA      NA" > $filename | cat $filename.temp >> $filename 
	echo "NA	$C	${ChrLength[$index]}	${ChrLength[$index]}	NA	NA      NA" >> $filename 
	echo "the head 3 lines" 
	head -n 3 $filename
	echo "the end 3 lines"
	tail -n 3 $filename
	let "C++"
	rm $filename.temp
done

