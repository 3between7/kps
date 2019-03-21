import pandas as pd
import sys,re



'''
import pandas as pd
a ={'a':[1,2,3,4], 'b':[5,6,7,8]}
b = {'c':[6,8], 'd':[11 ,12]}
aa = pd.DataFrame(a)
bb = pd.DataFrame(b)
cc = pd.merge(aa,bb, how='left' , left_on='b' , right_on='c' )
cc.loc[cc['c']>0,'b'] = cc[cc['c']>0]['d']

'''


f1=open("/media/jason/Seagate Backup Plus Drive/111111111111111111/Cattle_tiling1.id")
f2=open("/media/jason/Seagate Backup Plus Drive/111111111111111111/Cattle/fst/Cattle_Tiling2.id")


data=pd.read_table("/media/jason/Seagate Backup Plus Drive/111111111111111111/Cattle/6-formatDataSet1/dataset1_split/Cattle_Allpopulation_dataset1.lackGT.c"+sys.argv[1]+".step1",sep="\t")
print(data[1:10])
tiling1_id=[];tiling1_t=[]
tiling2_id=[];tiling2_t=[]
# make list of tiling1_id
for eachline in f1:
    if "Chr"+sys.argv[1]==eachline.strip().split("_")[0]:
        tiling1_id.append(eachline.strip())

# make list of tiling1_num
n1=len(tiling1_id)
for x in range(n1):
    tiling1_t.append(1)
print(n1,len(tiling1_t))
# make list of tiling2_id
for eachline2 in f2:
    if "Chr"+sys.argv[1]==eachline2.strip().split("_")[0]:
        tiling2_id.append(eachline2.strip())
# make list of tiling2_num
n2=len(tiling2_id)
for y in range(n2):
    tiling2_t.append(2)
print(n2,len(tiling2_t))
# make DataFrame of t1 and t2
data_t1=pd.DataFrame({'i':tiling1_id,'j':tiling1_t})
data_t2=pd.DataFrame({'i':tiling2_id,'j':tiling2_t})
# merge t2 and data
merge1=pd.merge(data,data_t2,how='left',left_on='b',right_on='i')
# replace tilingOrder2
merge1.loc[merge1['j']>0,'d'] = merge1[merge1['j']>0]['j']
print(merge1[1:10])
# merge t1 and new data  
data_new=merge1.loc[:,'a':'h']
merge2=pd.merge(data_new,data_t1,how='left',left_on='b',right_on='i')
# replace tilingOrder1
merge2.loc[merge2['j']>0,'d'] = merge2[merge2['j']>0]['j']
print(merge2[1:10])
merge2.to_csv("c"+sys.argv[1]+".step2",columns=['a','b','c','d','e','f','g','h'],sep="\t",header=0,index=0)
print("the chrom ",sys.argv[1],"is done!")



'''
for i in range(1,31):
	f1=open("Cattle_tiling1.id")
	f2=open("Cattle_Tiling2.id")
	tiling1=[]
	tiling2=[]
	for eachline in f1:
		if "Chr"+str(i)==eachline.strip().split("_")[0]:
			tiling1.append(eachline.strip())
	for eachline2 in f2:
		if "Chr"+str(i)==eachline2.strip().split("_")[0]:
			tiling2.append(eachline2.strip())	
	fx=open("./dataset1_split/Cattle_Allpopulation_dataset1.lackGT.c"+str(i)+".step1")
	fr=open("./dataset1_split/Cattle_Allpopulation_dataset1.lackGT.c"+str(i)+".step2","w")

	for eachline in fx:
		ll=re.split("\t",eachline.strip())
		if ll[1] in tiling1:
			fr.write(ll[0]+"\t"+ll[1]+"\t"+ll[2]+"\t"+"1"+"\t"+ll[4]+"\t"+ll[5]+"\t"+ll[6]+"\t"+ll[7]+"\n")
		elif ll[1] in tiling2:
			fr.write(ll[0]+"\t"+ll[1]+"\t"+ll[2]+"\t"+"2"+"\t"+ll[4]+"\t"+ll[5]+"\t"+ll[6]+"\t"+ll[7]+"\n")
		else:
			fr.write(eachline)
	print("the chrom ",str(i),"is done!")
	fx.close()
	fr.close()

'''

