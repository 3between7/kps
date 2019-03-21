import re
f=open("tansheep_vcf.freq.frq")
f1=open("tansheep_final_SNP.frq","w")
f2=open("tansheep_final_INDEL.frq","w")
f.readline()
for eachline in f:
	linelist=re.split(r"\s+",eachline.strip())
	if len(linelist[2].strip())==1 and len(linelist[3].strip()) ==1:
		f1.write(eachline)
	else:
		f2.write(eachline)
f.close()
f1.close()
f2.close()
