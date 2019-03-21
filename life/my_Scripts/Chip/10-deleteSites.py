import random,sys,re

def deleteSites(file,ResultFile):
	f=open(file) # Order	id     priority
	f1=open(ResultFile,"w")
	priority_list=[]
	locus_temp=[]
	ts=0
	temp=[]
	filePOS=0

	f.seek(filePOS)
	for eachline in f:
		llist=re.split(r"\s+",eachline.strip())
		N=llist[0]
		if int(int(N)/30) == ts:
			locus_temp.append(llist[1])
			priority_list.append(int(llist[2]))
			linelength=len(eachline)
			filePOS+=linelength
		else:
			ts+=1
			minPriority=max(priority_list)
			for i in range(30):
				if priority_list[i] == minPriority:
					temp.append(i)
					index=random.choice(temp)
					#print(locus_temp[index],minPriority)
			f1.write(locus_temp[index]+"\n")
			temp=[]
			priority_list=[]
			locus_temp=[]
			f1.close()

			f.seek(filePOS-linelength)

	f.close()
	f1.close()

if __name__ == "__main__":	
	
	deleteSites("Soybean_Chr"+sys.argv[1]+"_new.txt","Soybean_AllDelSites"+sys.argv[1]+".txt")





