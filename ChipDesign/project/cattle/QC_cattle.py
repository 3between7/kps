# -*- coding: UTF-8 -*-
import re,random,sys,os

def getCounts(file,filew):
	listID=[]
	f=open(file)
	fw=open(filew,"w")
	for eachline in f:
		print(eachline)
		ll=re.split(r"\s+",eachline.strip())
		listID.append(ll[0]);listID.append(ll[2])
	f.seek(0)
#	fw.write(title+"\t"+"id1_pairscounts"+"\t"+"id2_pairscounts")
	for eachline in f:
		ll=re.split(r"\s+",eachline.strip())
		fw.write(eachline.strip()+"\t"+str(listID.count(ll[0]))+"\t"+str(listID.count(ll[2]))+"\n")
	f.close()
	fw.close()

def delSamples(file,filew):
	f=open(file)
	fw=open(filew,"w")
	removeID=[]
	for eachline in f:
		ll=re.split(r"\s+",eachline.strip())
		if (ll[0] or ll[2]) in removeID:
			continue
		else:
			if int(ll[-2]) > int(ll[-1]):
				removeID.append(ll[0]+"\t"+ll[0]+"\n")
			elif int(ll[-2]) < int(ll[-1]):
				removeID.append(ll[2]+"\t"+ll[2]+"\n")
			else:
				choice=random.choice([ll[0],ll[2]])
				removeID.append(choice+"\t"+choice+"\n")
	uniqdelID=set(removeID)
	fw.writelines(uniqdelID)
	f.close()
	fw.close()
	
	
def qc(path,in_file,out_file,pihat=0.15):
		

	os.system("plink --bfile "+path+in_file+" --genome --chr-set 30 --allow-extra-chr --threads 5 --out "+path+in_file)
			
	os.system("sed -i '1d' "+path+in_file+".genome")
			
	os.system("sort -k 10nr "+path+in_file+".genome > "+path+in_file+".genome.sorted")
			
	f=open(path+in_file+".genome.sorted")
	fl=f.readline();
	f.close()
	li=re.split(r"\s+",fl.strip())
	if float(li[9])<=pihat:
		print("the current file is already statisfied your ask!")			
		exit()
	else:
		ft=open(path+in_file+".fam")	
		#N=os.popen("wc -l "+path+in_file+".fam")
		#print(str(N.read()).split(r"\s+")[0])
		nums=ft.read().count("\n");print(nums)
		if nums <48:
			print("the number of samples is:",nums,"stop running!")
			exit()
		elif nums <98 and nums >=48:
			print("the number of samples is:",nums,";the highest pi-hat is :",fl.split("\s+")[9]," go on or not (y/n)?")
			choice=input()
			if choice =="y":
				pass
			else:
				exit()
		else:
			print("there are ",nums,"samples remain,keep going")	
				
		os.system("awk '{if ($10>="+str(pihat)+"){print $0}}' "+path+in_file+".genome.sorted > "+path+in_file+".genome.sorted."+str(pihat))
				
		getCounts(path+in_file+".genome.sorted."+str(pihat),path+in_file+".genome.sorted."+str(pihat)+".classify")
				
		os.system("sort -k 10nr -k 15nr -k 16nr "+path+in_file+".genome.sorted."+str(pihat)+".classify > "+path+in_file+".genome.sorted."+str(pihat)+".classify.sorted")
				
		delSamples(path+in_file+".genome.sorted."+str(pihat)+".classify.sorted",path+in_file+"_"+str(pihat)+".del")

		os.system("plink --bfile "+path+in_file+" --remove "+path+in_file+"_"+str(pihat)+".del --chr-set 30 --allow-extra-chr --make-bed --out "+path +out_file)

if __name__=="__main__":
	
	#Cattle 600k SNP chip design
	'''
	os.system("plink --vcf /lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/1-mergeVCF/BCM_raw_merged_293pops.vcf --make-bed --allow-extra-chr --chr-set 29 --keep-allele-order --out /lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/2-QC/BCM_raw_merged_293pops") 		
	qc("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/2-QC/","BCM_raw_merged_293pops","BCM_thepops_v1")
	qc("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/2-QC/","BCM_thepops_v1","BCM_thepops_v2")
	qc("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/2-QC/","BCM_thepops_v2","BCM_thepops_v3")
	qc("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/2-QC/","BCM_thepops_v3","BCM_thepops_v4")
	qc("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/2-QC/","BCM_thepops_v4","BCM_thepops_v5")
	qc("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/2-QC/","BCM_thepops_v5","BCM_thepops_v6")
	qc("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/2-QC/","BCM_thepops_v6","BCM_thepops_v7")
	'''
	############################
	
	#rice 60k,20190211

	#qc("/lustre/home/zhanghuanhuan/RICE/1-QC/","NB_merge","NB_thepops_v1")
	#qc("/lustre/home/zhanghuanhuan/RICE/1-QC/","NB_thepops_v1","NB_thepops_v2")


	



