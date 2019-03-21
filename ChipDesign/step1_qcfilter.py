 # -*- coding: UTF-8 -*-
'''
update on 2019-02-18
@author: zhanghuanhuan

'''

import re,random,sys,os

def getCounts(file,filew):
	listID=[]
	f=open(file)
	fw=open(filew,"w")
	for eachline in f:
		#print(eachline)
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
	
	
def qc(rawVCF_prefix,chr_set,out_file,pihat=0.15):

	os.system("plink2 --vcf " + rawVCF_prefix + ".vcf --make-bed --allow-extra-chr --chr-set " + str(chr_set) + " --out " + rawVCF_prefix)
	in_file=rawVCF_prefix
	out_file_prefix=out_file
	t=1

	while True:

		os.system("plink --bfile "+in_file+" --genome --chr-set "+str(chr_set) + " --allow-extra-chr --threads 5 --out " + in_file)
		
		os.system("sed -i '1d' " + in_file+".genome")
				
		os.system("sort -k 10nr " + in_file+".genome > " + in_file + ".genome.sorted")
				
		f=open(in_file+".genome.sorted")
		fl=f.readline();
		f.close()
		li=re.split(r"\s+",fl.strip())
		if float(li[9])<=pihat:
			print("the current file is already statisfied your ask!")			
			break
		else:
			ft=open(in_file+".fam")	
			nums=ft.read().count("\n")

			if nums >=48:
				print("!!Notice:the number of samples BETTER greater than 98 and MUST greater than 48!")
				print("the number of samples is:",nums,";the highest pi-hat is : ",li[9]," go on or not (y/n)?")
				choice=input()
				if choice =="y":
					os.system("awk '{if ($10>="+str(pihat)+"){print $0}}' "+ in_file +".genome.sorted > "+ in_file +".genome.sorted."+str(pihat))					
					getCounts(in_file+".genome.sorted."+str(pihat),in_file+".genome.sorted."+str(pihat)+".classify")					
					os.system("sort -k 10nr -k 15nr -k 16nr "+ in_file +".genome.sorted."+str(pihat)+".classify > "+ in_file +".genome.sorted."+str(pihat)+".classify.sorted")					
					delSamples(in_file+".genome.sorted."+str(pihat)+".classify.sorted",in_file+"_"+str(pihat)+".del")
					os.system("plink2 --bfile "+ in_file +" --remove " + in_file + "_"+str(pihat)+".del --chr-set "+str(chr_set)+" --allow-extra-chr --make-bed --out "+out_file_prefix+"_v"+str(t))
				else:
					break
			else:
				print("the number of samples is:",nums,";the highest pi-hat is : ",li[9]," stop running!")
				break

		in_file=out_file_prefix+"_v"+str(t)
		t+=1
def filatehwemaf(prefix,chr_set,hwe=0.0001,maf=0.01):
	
	'''
	*.info : chr	pos	REF	ALT	MAF
	'''
	
	os.system("plink2 --bfile " + prefix + " --hwe " + str(hwe) + " --make-bed --out " + prefix + "_hwe --chr-set " + str(chr_set))
	os.system("plink2 --bfile " + prefix + "_hwe" + " --maf " + str(maf) + " --make-bed --out " + prefix + "_hwe_maf --chr-set " + str(chr_set))
	os.system("plink --bfile " + prefix + "_hwe_maf" + " --freq" + " --out " + prefix + "_hwe_maf --chr-set " + str(chr_set)) 
	os.system("sed -i '1d' " + prefix + "_hwe_maf.frq")
	os.system("paste " + prefix + "_hwe" + "_maf.bim " + prefix + "_hwe_maf.frq | awk '{print $1,$4,$6,$5,$11}' > " + prefix + "_hwe_maf.info" )
	os.system("head -3 " + prefix + "_hwe" + "_maf.info")
	
	
def getSitesPos2(input,output):	
	
	'''
	output:chr	pos1	pos2	REF	ALT	MAF,pos1:the real position of sites
	'''
	
	with open(input) as f1:
		f2=open(output,'w')
		f1.readline()
		for eachline in f1:
			linelist=re.split(r"\s+",eachline.strip())
			print(linelist)
			chrx=linelist[0]
			pos1=int(linelist[1])
			REF=linelist[2]
			ALT=linelist[3]
			maf=linelist[4]
			if len(ALT) > len(REF):#insert
				f2.write(chrx+"\t"+str(pos1)+'\t'+str(pos1+1)+'\t'+REF+'\t'+ALT+'\t'+maf+'\n')
			if len(ALT) < len(REF):# delete
				pos2=pos1+len(REF)-1
				f2.write(chrx+"\t"+str(pos1)+'\t'+str(pos2)+'\t'+REF+'\t'+ALT+'\t'+maf+'\n')
			if len(ALT) == len(REF):#snp
				f2.write(chrx+"\t"+str(pos1)+'\t'+str(pos1)+'\t'+REF+'\t'+ALT+'\t'+maf+'\n')
		f2.close()


def splitbychrom(input,specie,chrlength,path):
	
	'''
	split files by chrom, output :chr	pos1	pos2	REF	ALT	MAF
	'''

	chr=0
	f=open(input)
	flen=open(chrlength)
	fstep=0
		
	for line in flen:
		print(line)
		lf=re.split("\s+",line.strip())
		chr=lf[0]
		chr_length = lf[1]
		
		filename=path+specie+"_Chr" + chr + ".txt";print(filename)
		fw=open(filename,"w")
		fw.write(chr+"\t"+"0"+"\t"+"0"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\n")
		f.seek(fstep)
		for fl in f:
			ll = re.split("\s+",fl.strip())
			if ll[0] == chr:
				#print(fl)
				fstep+=len(fl)	
				fw.write(fl)
			if ll[0] > chr:
				print(fl)
				break
		fw.write(chr +"\t"+chr_length+"\t"+chr_length+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\n")	
		fw.close()
			

def filate35bp(input,filew,type="SNP",flanklength=35):
	
	'''
	input:chr	pos1	pos2	REF	ALT	MAF
	Output:snpid	chr	pos1	pos2	REF	ALT	MAF
	'''
	
	f1=open(input)
	f2=open(filew,"w")
	listfile=list(f1)
	N=len(listfile)
	for i in range(1,N-1):
		linel=listfile[i-1]
		linel_list=re.split(r"\s+",linel.strip())
		posl_right=int(linel_list[2])
	
		linem=listfile[i]
		linem_list=re.split(r"\s+",linem.strip())
		
		if linem_list[1]!=linem_list[2]: #if pos1!=pos2,the type of sites could be indel
			posm_left=int(linem_list[1])+1 #the real position that mutation occurred!
			posm_right=int(linem_list[2])
		else:
			posm_left=int(linem_list[1]) #if pos1==pos2,sites are SNP!
			posm_right=int(linem_list[2])
			
		liner=listfile[i+1]
		liner_list=re.split(r"\s+",liner.strip())
		if liner_list[1]!=liner_list[2]:	
			posr_left=int(liner_list[1])+1
		else:
			posr_left=int(liner_list[1])
			
		if int(linem_list[1]) == int(linel_list[1]) or int(linem_list[1]) == int(liner_list[1]):
			continue
					
		disl=posm_left-posl_right-1
		disr=posr_left-posm_right-1
		
		if (disl >= flanklength and disr < flanklength) or (disl < flanklength and disr >= flanklength) or (disl >= flanklength and disr >= flanklength):
			if type =="SNP":
				if linem_list[1] == linem_list[2]:
					f2.write("Chr"+linem_list[0]+"_"+linem_list[1]+"\t"+listfile[i])
			else:
				f2.write("Chr"+linem_list[0]+"_"+linem_list[1]+"\t"+listfile[i])			
	f1.close()
	f2.close()
	

if __name__=="__main__":
	
	#rice 60k
	qc("/lustre/home/zhanghuanhuan/RICE/1-QC/","NB_thepops_v2","NB_thepops")


	



