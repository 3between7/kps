# -*- coding: UTF-8 -*-
'''
update on 2019-02-21
@author: zhanghuanhuan，liurui

'''

import random,re,pickle,copy
import pandas as pd
from ChipDesign.tools import uniqFlankSeq,window
#from ChipDesign.tools.Caculators import Caculator_multiselect



def slidewindow(input,chrnum,output,tilingnums,windowWidth,slideSize,chrlengthfile,Caculator):
	
    '''
    input/output:species    snpid    sequence    tiling_order    importance    chr    pos    REF    ALT    maf
    '''
    
    #get chr length
    dict_chrlen={}
    pickedSites=[]
    f_len=open(chrlengthfile)
    M=window.prepareforwindow(input) #DataSet1,need title
	#get length of chr
    for line in f_len:
        ll=re.split(r"\s+",line.strip())
        chr=ll[0]
        length=int(ll[1])
        dict_chrlen[chr]=length
	#silde window	
    for chr in range(chrnum):
        chrx=str(chr+1)	
        ListChr=M.getSitesListByChrom(chrx,startpos=9236)
        Mywindow=window.Window()
        c = Caculator(tilingnums)
        Mywindow.slidWindowOverlap(ListChr, dict_chrlen[chrx], windowWidth, slideSize, c)
        ss = copy.deepcopy(Mywindow.winValueL);#print("testnow",ss[:5])

        for eachwindow in ss:
            if eachwindow[3]:
                pickedSites.append(eachwindow[3][0])
                			
    dictresult = {'id':pickedSites}
    sites = pd.DataFrame(dictresult);print(sites[:3])
    dataset1=pd.read_csv(input,sep="\t");print(dataset1[:3])
    result=pd.merge(dataset1,sites,how="inner",left_on="snpid",right_on="id")		
    result.to_csv(output,columns=['Organism','snpid','seventyonemer','Tiling Order','importance','chromosome','POS','Ref Allele','Alt Allele','MAF'],header=1,index=0,sep="\t")		
    #get uniq sites
    uniqFlankSeq.getUniqSites(output)
    
    
    