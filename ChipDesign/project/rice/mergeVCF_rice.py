import re,os

def mergeVCF(*files):
	samples=[]
	dictINFO={}
	listchr=[];listpos=[];listID=[];listREF=[];listALT=[];listQUAL=[];listFILTER=[];listINFO=[];listFORMAT=[]
	fw=open("20181228-merged.vcf","w")
	N=0
	for file in files:
		f=open(file)
		for eachline in f:
			ll=re.split(r"\s+",eachline.strip());#print(ll)
			if "##" in eachline:
				continue
			elif "CHROM" in eachline:
				titlelist=re.split("\t",eachline.strip())
				N=len(titlelist)
				#print(N)
				sx=ll[9:]
				samples.extend(sx)
				thesample=set(samples)
			else:
				value=[]
				for i in range(9,N):
					key=titlelist[i].split("_")[0]+"_"+ll[0]+"_"+ll[1]
					value.append(ll[i])
					dictINFO[key]=value
				listchr.append(ll[0]);listpos.append(ll[1]);listID.append(ll[2]);listREF.append(ll[3]);listALT.append(ll[4]);listQUAL.append(ll[5]);listFILTER.append(ll[6]);listINFO.append(ll[7]);listFORMAT.append(ll[8])		
#	print(dictINFO.keys());print(dictINFO.values())			
	thesample=list(set(samples))
#	print(thesample)
	sampleinfo="\t".join(thesample);#print(sampleinfo)
	header="\t".join(("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",sampleinfo));#print(header)
	fw.write(header+"\n")
	rowslen=len(listchr)
	for i in range(rowslen):
		line=listchr[i]+"\t"+listpos[i]+"\t"+listID[i]+"\t"+listREF[i]+"\t"+listALT[i]+"\t"+listQUAL[i]+"\t"+listFILTER[i]+"\t"+listINFO[i]+"\t"+listFORMAT[i]
		#print(line)
		for samplex in thesample:
			#line+="\t"+dictINFO[samplex.split("_")[0]+"_"+listchr[i]+"_"+listpos[i]][0]
			try:
				line+="\t"+dictINFO[samplex.split("_")[0]+"_"+listchr[i]+"_"+listpos[i]][0]
			except:
				line+="\t"+"./."
		print(line)
		fw.write(line+"\n")
	fw.close()
	
def fillsamples(*files):
	samples=[]
	genotype=[]
	filewlist=["vcftools merge"]
	for file in files:
		fw=open(file+"_mordified.vcf","w")
		filenum=files.index(file)
		f=open(file)
		for eachline in f:
			ll=re.split(r"\s+",eachline.strip());print(ll)
			if "##" in eachline:
				continue
			elif "CHROM" in eachline:
				listname="list"+str(filenum)
				titlelist=re.split("\t",eachline.strip())
				listname=set(ll[9:])
				samples.extend(listname)
				fw.write(eachline)
			else:
				fw.write(eachline)
		fw.close()
		f.close()
	thesample=set(samples)
	for i in range(len(files)):
		diffsamplex="\t".join(list(thesample.difference("list"+str(i))))
		N=len(diffsamplex)
		for x in range(N):
			genotype.append("./.")
		genotype_string="\t".join(genotype)			
		f=open(files[i]+"_mordified.vcf")
		fw=open(files[i]+"_beforemerge.vcf","w")
		titlex=f.readline().strip()+"\t"+diffsamplex
		print(titlex)
		fw.write(titlex+"\n")
		for line in f:
			newline=line.strip()+"\t"+genotype_string
			print(newline)
			fw.write(newline+"\n")
		f.close();fw.close()
		os.system("bgzip -c " +files[i]+"_beforemerge.vcf > "+files[i]+"_beforemerge.vcf.gz")
		os.system("tabix -f "+files[i]+"_beforemerge.vcf.gz")
		filewlist.append(files[i]+"_beforemerge.vcf.gz")
	command="\s".join(filewlist)+" > rice_rawmerged.vcf.gz"
	print(command)
	os.system(command)

if __name__ == '__main__':
	
	'''
	mergeVCF("/home/zhanghuanhuan/TempDir/20181219/rice_allnotin.vcf",
			"/home/zhanghuanhuan/TempDir/20181219/Nipponbare_indel_mordified.vcf",
			"/home/zhanghuanhuan/TempDir/20181219/NB_final_snp_mordified.vcf",
			"/home/zhanghuanhuan/TempDir/20181219/JG7_mordified.vcf",
			"/home/zhanghuanhuan/TempDir/20181219/IG4_mordified.vcf",
			"/home/zhanghuanhuan/TempDir/20181219/IG1_mordified.vcf")	
	
	fillsamples("/home/zhanghuanhuan/TempDir/20181219/rice_allnotin.vcf",
			"/home/zhanghuanhuan/TempDir/20181219/Nipponbare_indel_mordified.vcf",
			"/home/zhanghuanhuan/TempDir/20181219/NB_final_snp_mordified.vcf",
			"/home/zhanghuanhuan/TempDir/20181219/JG7_mordified.vcf",
			"/home/zhanghuanhuan/TempDir/20181219/IG4_mordified.vcf",
			"/home/zhanghuanhuan/TempDir/20181219/IG1_mordified.vcf")

	###20181228
	mergeVCF("/lustre/home/zhanghuanhuan/Data/ChipDesign/Rice/18-report/mergeVCF/20181217-60828-RiceSitesInChip.sorted.vcf",
			"/home/zhanghuanhuan/TempDir/20181227/12-tile1-p8-indel.chrpos.recode.vcf",
			"/home/zhanghuanhuan/TempDir/20181227/6864-tile1-p8-snp.recode.vcf")
	'''
	
	###20190104 去掉Pan基因组与5个notin位点
	mergeVCF("/home/zhanghuanhuan/TempDir/20181219/Nipponbare_indel_mordified.vcf",
			"/home/zhanghuanhuan/TempDir/20181227/12-tile1-p8-indel.chrpos.recode.vcf",
			"/home/zhanghuanhuan/TempDir/20181219/NB_final_snp_mordified.vcf",
			"/home/zhanghuanhuan/TempDir/20181227/6864-tile1-p8-snp.recode.vcf")
	
