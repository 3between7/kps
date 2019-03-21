#-*- coding: UTF-8 -*-
import re,pickle
listX=[]
for i in range(1,31):
	listname="C"+str(i)+"pickedlist"
	listname=pickle.load(open("/media/jason/Seagate Backup Plus Drive/111111111111111111/Cattle/8-slideWindow/Cattle_Chr"+str(i)+"picked.pickle","rb"))   
	listX.extend(listname)

# [['id1'],["id2"],[]...]

#try:
#    dictAllRecord=pickle.load(open("CattledictSites.pickle","rb"))
#    print(dictAllRecord.keys()[1:10])
#except:
#    dictAllRecord={}
#    f=open("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/6-formatDataSet1/Cattle_Allpopulation_dataset1.lackGT")
#    for eachline in f:
#        linelist=re.split(r"\t",eachline.strip())
#        key=linelist[1].strip()
#        dictAllRecord[key]=eachline.strip()
#        pickle.dump(dictAllRecord,open("CattledictSites.pickle","wb"))

dictAllRecord={}
f=open("/media/jason/Seagate Backup Plus Drive/Cattle_DataSet1_update.txt")
for eachline in f:
    linelist=re.split(r"\s+",eachline.strip())
    key=linelist[1].strip()
#    print(key)
    dictAllRecord[key]=eachline.strip()
#    pickle.dump(dictAllRecord,open("CattledictSites.pickle","wb"))


def getpick(listX):
    for eachid in listX:
        if eachid!=[]:
            K=eachid[0]
            print(dictAllRecord[K])
        
getpick(listX)
f.close()

        
