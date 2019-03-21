# -*- coding: UTF-8 -*-
'''
update on 2019-02-14
@author: zhanghuanhuan

'''

import re,os

####################################################################################################################

'''
Format file 
'''
	
def changeVCFformat(input,output,sexchr="XY",sexcode="9999"):
    f=open(input)
    fw=open(output,"w")
    for eachline in f:
        if "##" in eachline:
            fw.write(eachline)
        elif "#CHROM" in eachline:
            replacedline=re.sub(r"_","",eachline)
            fw.write(replacedline)
        elif sexchr in eachline:
        	replacedline=re.sub(sexchr,sexcode,eachline)
        	fw.write(replacedline)
        else:
        	fw.write(eachline)
    '''
    	#this code is for rice
    	
        if "gnl|UMD3.1|GK" in eachline:
            if "gnl|UMD3.1|GK000030.2" in eachline:
                chromsome="X"
                replacedline=re.sub(r"(gnl\|UMD3\.1\|GK0000[0-9]{2}\.2)",chromsome,eachline)
            else:
                  chromsome=str(int(eachline.strip().split("\t")[0][-4:-2]))
                  replacedline=re.sub(r"(gnl\|UMD3\.1\|GK0000[0-9]{2}\.2)",chromsome,eachline)
            fw.write(replacedline)      
        elif "31" == eachline.strip().split("\t")[0]:
            replacedline=re.sub(r"^(31){1}[\s+]","30\t",eachline)
            fw.write(replacedline)
        else:
            fw.write(eachline)
	'''
    f.close()
    fw.close()

###################################################################################################################

'''
The duplicate sample is the same sample
'''

def mergeVCF_way1(output,*files):

	'''
	Suitable for: (1) Samples are totally different; (2) Loci are totally different.  
	'''
	samples=[]
	dictINFO={}
	listchr=[];listpos=[];listID=[];listREF=[];listALT=[];listQUAL=[];listFILTER=[];listINFO=[];listFORMAT=[]
	fw=open(output+"_temp","w")
	N=0
	f_t=open("/home/zhanghuanhuan/Data/ChipDesign/example/checkdict.txt","w")
	for file in files:
		f=open(file)
		for eachline in f:
			ll=re.split(r"\s+",eachline.strip())
			if "##" in eachline:
				#pass annonated lines in vcf 
				continue
			elif "CHROM" in eachline:
				#get all samples in all files
				titlelist=re.split("\s+",eachline.strip())
				N=len(titlelist)
				sx=ll[9:]
				samples.extend(sx)
				thesample=set(samples)
			else:
				#get GT of samples
				value=[]
				for i in range(9,N):
					#pay attention:the FID of sample couldn't include "_"
					key=titlelist[i].split("_")[0]+"_"+ll[0]+"_"+ll[1]  #sampleName_chrom_pos,for example:MG-1_1_300 
					value=ll[i] #GT of the sample at chr/pos
					print(key,value,sep="\t",file=f_t)
					dictINFO[key]=value
				listchr.append(ll[0]);listpos.append(ll[1]);listID.append(ll[2]);
				listREF.append(ll[3]);listALT.append(ll[4]);listQUAL.append(ll[5]);
				listFILTER.append(ll[6]);listINFO.append(ll[7]);listFORMAT.append(ll[8])		
	#merge info from mutiple vcf	
	thesample=list(set(samples)) #MG-1_MG-1
	sampleinfo="\t".join(thesample)
	header="\t".join(("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",sampleinfo))
	fw.write(header+"\n")
	rowslen=len(listchr)
	for i in range(rowslen):
		line=listchr[i]+"\t"+listpos[i]+"\t"+listID[i]+"\t"+listREF[i]+"\t"+listALT[i]+"\t"+listQUAL[i]+"\t"+listFILTER[i]+"\t"+listINFO[i]+"\t"+listFORMAT[i]
		for samplex in thesample:
			try:
				line+="\t"+dictINFO[samplex.split("_")[0]+"_"+listchr[i]+"_"+listpos[i]]
			except:
				line+="\t"+"./."
		#print(line)
		fw.write(line+"\n")
	fw.close()
	os.system("cat " + output + "_temp | vcf-sort | uniq > " + output)
	
def mergeVCF_way2(temp_vcf,output,*files):
	
	'''
	Better for merging large files  
	'''
	samples=[]
	genotype=[]
	titlelist=[]
	sl=[]
	command_pre=["vcf-concat"]
	
	#remove annonated lines in vcf and collect all samples 
	for file in files:
		#fw=open(file.split(".")[0]+"_mordified.vcf","w")
		filenum=files.index(file)
		f=open(file)
		for eachline in f:
			ll=re.split(r"\s+",eachline.strip())
			if "##" in eachline:
				continue
			if "CHROM" in eachline:
				titlelist.append(eachline.strip())
				sl.append(set(ll[9:]))
				samples.extend(set(ll[9:]))
		f.close()
		
	#get Get the difference sample,make the GT of them "./."
	thesample=set(samples);#print(thesample)
	for i in range(len(files)):
		diffsamplex="\t".join(list(thesample.difference(sl[i])));print(diffsamplex)
		N=len(list(thesample.difference(sl[i])));print(N)
		for x in range(N):
			genotype.append("./.")
		genotype_string="\t".join(genotype)			
		f=open(files[i].split(".")[0]+".vcf")	
		fw=open(files[i].split(".")[0]+"_mor.vcf","w")
		titlex=titlelist[i]+"\t"+diffsamplex ################################
		for line in f:
			if "##" in line:
				fw.write(line)
			elif "#CHROM" in line:
				fw.write(titlex+"\n")
			else:
				newline=line.strip()+"\t"+genotype_string
				print(newline)
				fw.write(newline+"\n")
		f.close();fw.close()
		os.system("bgzip -c " +files[i].split(".")[0]+"_mor.vcf > "+files[i].split(".")[0]+"_mor.vcf.gz")
		genotype=[]
		
	#reorder the columns of vcf,something wrong in there!
	for  f in files:
		if f != temp_vcf:
			os.system("vcf-shuffle-cols -t " + temp_vcf.split(".")[0] + "_mor.vcf.gz" + f.split(".")[0]+"_mor.vcf.gz > " + f.split(".")[0]+"_mor_reor.vcf")
			os.system("bgzip -c " +f.split(".")[0]+"_mor_reor.vcf > "+f.split(".")[0]+"_mor_reor.vcf.gz")
			os.system("tabix -f "+f.split(".")[0]+"_mor_reor.vcf.gz")
			command_pre.append(f.split(".")[0]+"_mor_reor.vcf.gz")
		else:
			os.system("tabix -f "+f.split(".")[0]+"_mor.vcf.gz")
			command_pre.append(f.split(".")[0]+"_mor.vcf.gz")
	
	#concat vcf		
	command=" ".join(command_pre)+" > "+output
	print(command)
	os.system(command)

##################################################################################################################


if __name__ == '__main__':
	
	changeVCFformat("/home/zhanghuanhuan/test.vcf","X","30","/home/zhanghuanhuan/test_mor.vcf")

	
