# -*- coding: UTF-8 -*-
import re

with open("Nipponbare_indel_mergeMapFrq_INFO.txt") as f1:
	f2=open("INDEL_INFO.txt",'w')
	for eachline in f1:
		linelist=re.split(r"\t",eachline.strip())
		pos1=int(linelist[1])
		ALT=linelist[3]
		REF=linelist[2]
		if len(ALT) > len(REF):#insert
			f2.write(linelist[0]+'\t'+str(pos1+1)+'\t'+str(pos1+1)+'\t'+REF+'\t'+ALT+'\t'+linelist[4]+'\n')
		if len(ALT) < len(REF):# delete
			pos2=pos1+len(REF)-1
			f2.write(linelist[0]+'\t'+str(pos1+1)+'\t'+str(pos2)+'\t'+REF+'\t'+ALT+'\t'+linelist[4]+'\n')
		if len(ALT) == len(REF):#snp
			f2.write(linelist[0]+'\t'+str(pos1)+'\t'+str(pos1)+'\t'+REF+'\t'+ALT+'\t'+linelist[4]+'\n')
	f2.close()
		

		
	
