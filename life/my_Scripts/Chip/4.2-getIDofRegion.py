#-*- coding: UTF-8 -*-
import re

def getlist(file):
    f=open(file) # chr pos1 pos2
    chr=[]
    pos1=[]
    pos2=[]
    for each in f:
        linelist=re.split(r"\t",each.strip())
        chr.append(linelist[0])
        pos1.append(linelist[1])
        pos2.append(linelist[2])
    f.close()
    return chr,pos1,pos2
            
def getid(file,chr,pos1,pos2): # chr=[] pos1=[] pos2=[]
    n=len(chr)
    f1=open("Cattle_Allpopulation_dataset1.priority","w")
    filepos=0
    f=open(file)
    fl=len(f.read())
    f.seek(0)
    for i in range(n):
        chri=chr[i]
        pos1i=int(pos1[i])
        pos2i=int(pos2[i])
        print(chri,pos1i,pos2i)
#        print(filepos)
        f.seek(filepos)
        for line in f:
            filepos+=len(line)
            linelist1=re.split(r"\t",line.strip())
            if linelist1[0]==chri:
                left_dis=int(linelist1[1])-pos1i
                right_dis=pos2i-int(linelist1[1])
                if (left_dis >= 0) and (right_dis >=0):
#                    print(linelist1)
                    print("chr of SNP is ",linelist1[0],"pos of SNP is ",linelist1[1])
                    newline=linelist1[0]+"\t"+linelist1[1]+"\t"+linelist1[2]+"\t"+linelist1[3]+"\t"+linelist1[4]+"\t"+"1"+"\n"
                    print(newline)
                    f1.write(newline)
                else:
                    f1.write(line)
            else:
                filepos=filepos-len(line)
                break
        print((filepos/fl)*100,"%")    
    f.close()
    f1.close()  
    
    
chr,pos1,pos2=getlist("QTLregion.sort")
getid("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/4-classifyMAF/Cattle_Allpopulation_dataset1.classfy",chr,pos1,pos2)  #info format:id chr pos1 pos2

