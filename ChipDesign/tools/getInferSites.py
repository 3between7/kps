# -*- coding: UTF-8 -*-
'''
update on 2019-02-21
@author: zhanghuanhuan

'''

import pandas as pd
import os

def getInfer(species,input,DataSet1,output):

    #load DATASET1
    os.system('''awk 'BEGIN{FS="\t";OFS="\t"}{print $2}' ''' + input + " > PickedSitesid.txt")
    dataset1=pd.read_csv(DataSet1,sep="\t")
    #load final sites
    fs=pd.read_csv("PickedSitesid.txt",sep="\t")
    dup=dataset1.append(fs)
    dups=dup.append(fs)
    result=dups.drop_duplicates(subset=['snpid'],keep=False)
    result.to_csv(output+"_tmp",columns=['snpid','seventyonemer'],header=1,index=0,sep="\t")
    os.system('''awk 'BEGIN{FS="\t";OFS="\t"}{print ''' + species + ",$1,$2}' " + output + "_tmp > " + output)
    os.remove(output+"_tmp")
    os.remove("PickedSitesid.txt")