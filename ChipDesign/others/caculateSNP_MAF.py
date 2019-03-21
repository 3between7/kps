import re
f=open("tansheep_final_SNP.frq")
f1=open("tansheep_final_SNP.frq.noesx","w")
NumOf_1=0
NumOf_2=0
NumOf_3=0
NumOf_4=0
NumOf_5=0
try:
	for eachline in f:
		linelist=re.split(r"\s+",eachline.strip())
		MAF=float(linelist[4].strip())
		if MAF >0 and MAF < 0.1:
			NumOf_1+=1
		elif MAF >= 0.1 and MAF <0.2:
			NumOf_2+=1
		elif MAF >= 0.2 and MAF < 0.3:
			NumOf_3+=1
		elif MAF >=0.3 and MAF < 0.4:
			NumOf_4 +=1
		elif MAF >=0.4 and MAF <= 0.5:
			NumOf_5 +=1
		else:
			f1.write(eachline)
except:
	f1.write(eachline)

print("0-0.1:",NumOf_1)
print("0.1-0.2:",NumOf_2)
print("0.2-0.3:",NumOf_3)
print("0.3-0.4:",NumOf_4)
print("0.4-0.5:",NumOf_5)
f.close()
f1.close()
