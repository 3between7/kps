import re

f=open("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/3-filate35bp/Cattle_Allpopulation_dataset1.info")
f1=open("Cattle_Allpopulation_dataset1.classfy","w")
for eachline in f:
	ll=re.split("\t",eachline.strip())
	try:
		MM=float(ll[2].split(":")[2])
	except:
		MM=0
	try:
		Bf=float(ll[3].split(":")[2])
	except:
		Bf=0
	try:
		CY=float(ll[4].split(":")[2])
	except:
		CY=0
	if MM >= 0.05 and Bf>=0.05 and CY>=0.05:
		f1.write(eachline.strip()+"\t"+"2"+"\n")
	elif (MM >=0.05 and Bf>=0.05 and CY<0.05) or (MM >=0.05 and Bf<0.05 and CY>=0.05) or (MM<0.05 and Bf>=0.05 and CY>=0.05):
		f1.write(eachline.strip()+"\t"+"3"+"\n")
	elif (MM >=0.05 and Bf <0.05 and CY<0.05) or (MM<0.05 and Bf>=0.05 and CY<0.05) or (MM<0.05 and Bf<0.05 and CY>=0.05):
		f1.write(eachline.strip()+"\t"+"4"+"\n")
	else:
		f1.write(eachline.strip()+"\t"+"5"+"\n")
f.close()
f1.close()
