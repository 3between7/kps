import re
f=open("/lustre/project/R2D/Rice60k-CAU/18-report/mergeVCF/result/raw/20190104-riceRaw-sort.vcf.tped")
#f=open("/lustre/project/R2D/Rice60k-CAU/18-report/mergeVCF/result/raw/20190104-riceRaw-sort_mordififed_test.tped")
for eachline in f:
	if "un" in eachline:
		ll=re.split("\t",eachline)
		ll[0]=ll[1].split(":")[0]
		ll[1]=ll[1].split(":")[1]
		#r1=re.sub(ll[0],ll[1].split(":")[0],eachline.strip())
		#r2=re.sub(ll[1],ll[1].split(":")[1],r1)
		new="\t".join(ll)
		print(new)
	else:
		print(eachline.strip())
