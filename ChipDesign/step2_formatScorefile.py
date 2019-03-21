# -*- coding: UTF-8 -*-
'''
update on 2019-02-18
@author: zhanghuanhuan

'''
import pandas as pd
import os,pickle,re
from ChipDesign.tools import getSiteSeq


###############################################################################
#Step1：注释优先级


'''
annonate tiling order
'''

def annoTlingOrder(nums,input,output,tiling=[]):
    
    '''
    input file need title: snpid,chromosome,POS,REF,ALT,maf
    tiling file: [snpid tn,...,snpid t2,snpid t1]
    output file: snpid    chr    POS1    POS2    REF    ALT    maf    t(n)
    '''
    
    temp=output.split(".")[1]+"_temp.txt"
    #load file
    m=0
    dataset=pd.read_table(input,sep="\t") #title :snpid chromosome POS POS2 Ref Allele Alt Allele maf
    te=pd.read_table(tiling[m],sep="\t")
    result=pd.merge(dataset,te,how="left",left_on="snpid",right_on='snpid');print(result[:3])
    m=1
    
    while m <= int(nums-1):
        tc=int(nums)-m
        te=pd.read_table(tiling[m],sep='\t')
        result=pd.merge(result,te,how="left",left_on="snpid",right_on='snpid');print(result[:3])
        result.loc[result["t"+str(tc)]>0,"t"+str(nums)] = result.loc[result["t"+str(tc)]>0,"t"+str(tc)]
        m+=1
        
    result.to_csv(temp,columns=['snpid','chromosome','POS','POS2','Ref Allele','Alt Allele','maf',"t"+str(nums)],sep="\t",header=1,index=0)
    command= ''' awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,int($8)}' ''' + temp + " > " + output
    os.system(command) 
    os.remove(temp)
 
    
###############################################################################
#Step2:获得侧翼序列,并整理格式

def DataSet1(species,input,refFastaFileName,chrlength,output,vip=None,type="A",flanklength=35):
    
    '''
    input: snpid    chr    POS1    POS2    REF    ALT    maf    t(n)
    output: species    snpid    sequence    tiling_order    importance    chr    pos    REF    ALT    maf
    '''
    lengthdict={};poslist=[];vipsites=[]
    f=open(refFastaFileName)
    f1=open(input)
    f1.readline()
    f2=open(chrlength)
    f3=open(output,"w")
    
    try:
        refidx=pickle.load(open(refFastaFileName+".myfasteridx",'rb'))
    except:
        getSiteSeq.generateFasterRefIndex(refFastaFileName, refFastaFileName+".myfasteridx",mapname=None,startchar=">",chrsignal=None,romanSignal=False)
        
    try:
        vip=open(vip)
        for s in vip:
            vipsites.append(s.strip())
    except:
        pass

    for eachline in f2:
        linelist1=re.split(r"\s+",eachline.strip())
        lengthdict[linelist1[0]]=linelist1[1]
    
    f3.write("Organism\tsnpid\tseventyonemer\tTiling Order\timportance\tchromosome\tPOS\tRef Allele\tAlt Allele\tMAF\n")
    if type == "A":
        for line in f1:
            ll=re.split(r"\t",line.strip());print(ll)
            chrom=ll[1]
            REF=ll[4]
            ALT=ll[5]
            lengthREF=len(ll[4])
            lengthALT=len(ll[5])
            chrlength=int(lengthdict[ll[1]])
            importance="VIP" if ll[0] in vipsites else "Standard"
            
            if lengthREF == lengthALT: #SNP
                RefSeqMap = getSiteSeq.getRefSeqBypos_faster(f,refidx,chrom,int(ll[2])-flanklength,int(ll[2])+flanklength,chrlength)
                f3.write(species+"\t"+ll[0]+"\t"+"".join(RefSeqMap[chrom][1:36]).upper()+"["+REF+"/"+ALT+"]"+"".join(RefSeqMap[chrom][37:]).upper()+"\t"+ll[7]+"\t"+importance+"\t"+chrom+"\t"+ll[2]+"\t"+REF+"\t"+ALT+"\t"+ll[6]+"\n")    
            elif lengthREF < lengthALT: # insert
                lengthinsert=len(ll[5])-1
                RefSeqMap = getSiteSeq.getRefSeqBypos_faster(f,refidx,chrom,int(ll[2])-flanklength,int(ll[3])+flanklength-1,chrlength)
                f3.write(species+"\t"+ll[0]+"\t"+"".join(RefSeqMap[chrom][1:36]).upper()+"["+"-"+"/"+ALT[1:]+"]"+"".join(RefSeqMap[chrom][36:]).upper()+"\t"+ll[7]+"\t"+importance+"\t"+chrom+"\t"+ll[2]+"\t"+REF+"\t"+ALT[1:]+"\t"+ll[6]+"\n")
            else: # delete
                lengthdel= len(ll[4])-1
                RefSeqMap = getSiteSeq.getRefSeqBypos_faster(f,refidx,chrom,int(ll[2])-flanklength,int(ll[3])+flanklength,chrlength)
                p2=36+lengthdel
                f3.write(species+"\t"+ll[0]+"\t"+"".join(RefSeqMap[chrom][1:36]).upper()+"["+REF[1:]+"/"+"-"+"]"+"".join(RefSeqMap[chrom][p2:]).upper()+"\t"+ll[7]+"\t"+importance+"\t"+chrom+"\t"+ll[2]+"\t"+REF[1:]+"\t"+ALT+"\t"+ll[6]+"\n")
                
    if type=="I":
        '''
        Yet to be developed!
        '''
        pass        
              
    f.close();f1.close();f2.close();f3.close()
    
    