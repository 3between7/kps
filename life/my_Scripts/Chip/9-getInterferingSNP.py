import re,pandas,sys
f=open("/media/jason/Seagate Backup Plus Drive/111111111111111111/Cattle/2-mergeMAF/result_3population.c"+sys.argv[1]+".sort")
f1=open("/media/jason/Seagate Backup Plus Drive/111111111111111111/Cattle/10-getInterferingSNP/InterferingSNPGT.c"+sys.argv[1],"w")
f2=open("/media/jason/Seagate Backup Plus Drive/111111111111111111/Cattle/9-getRecord/CattlePickedRecord_v1.C"+sys.argv[1]+"id")
listid=[]
listGT=[]
listidPicked=[]

for eachline in f:
	ll=re.split(",",eachline.strip())
	id= ll[0]+"_"+ll[1]
	listid.append(id)
	
	if ll[2]:
		listGT.append(ll[2].split(":")[1])
	elif ll[3]:
		listGT.append(ll[3].split(":")[1])
	else:
		listGT.append(ll[4].split(":")[1])

print("listid and listGT is get!")

for eachline in f2:
	listidPicked.append(eachline.strip())

print("the list of Sites that picked is get!")
print("===============================================================")	


data3Pops={'1_id':listid,'2_GT':listGT}
data1id={'1_id':listidPicked}

a=pandas.DataFrame(data3Pops)
print("data3Pops[1:10]")
print(a[1:10])
print("===============================================================")

b=pandas.DataFrame(data1id)
print("ddata1id[1:10]")
print(b[1:10])
print("===============================================================")

a=a.append(b,sort="True")
result = a.drop_duplicates(subset=['1_id'],keep=False)
print("result[1:10]")
print(result[1:10])
print("===============================================================")

result.to_csv(f1,index=False,header=False)

print("the length of listid is :",len(listid),"the length of result is :",len(result),"the length of picked is :",len(listidPicked))

f.close()
f1.close()
f2.close()
	
		

