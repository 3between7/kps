# -*- coding: UTF-8 -*-
'''
update on 2019-02-19
@author: zhanghuanhuan，liurui

'''

import os,re,sys

def complementary(seq):
    """
    ['tg', 'a', 't', 'g', 'c', 'acacacgatg', 'ctttttcccccccc', 'c', 'c', 'a', 'a', 'aaagagagagacagaaaaaggc', 'atatcgactg', 'catcga']
    reverse to 
    ['ac', 't', 'a', 'c', 'g', 'tgtgtgctac', 'gaaaaagggggggg', 'g', 'g', 't', 't', 'tttctctctctgtctttttccg', 'tatagctgac', 'gtagct']
    """
    newseq = []
    for i in range(0, len(seq)):
        if seq[i].lower() == 'a':
            newseq.insert(i, 't')
        elif seq[i].lower() == 't':
            newseq.insert(i, 'a')
        elif seq[i].lower() == 'c':
            newseq.insert(i, 'g')
        elif seq[i].lower() == 'g':
            newseq.insert(i, 'c')
        elif len(seq[i]) > 1:
            newseq.insert(i, complementary(seq[i]))
        else:
            newseq.insert(i, seq[i])
    if isinstance(seq, str):
        newseq = "".join(newseq)
    return newseq


def getUniqSites(input,seqidx=2):
    
    '''
    input:species    snpid    sequence    tiling_order    importance    chr    pos    REF    ALT    maf
    '''
    
    infile1map={};dupmap={};dupids=[] 
    
    tilingorder_idx=2
    f=open(input,'r');f.readline();ofo=open(input+"dup",'w')
    
    for line in f:
        recl=re.split(r"\s+",line.strip());#print(recl)
        seql=re.split(r"\[.+\]",recl[seqidx].strip());print(seql)
        flankseq=list(seql[0]+seql[1]);flankseq.reverse()
        compflankseq=complementary(flankseq)
        print(seql[0]+seql[1],compflankseq)
        if seql[0]+seql[1] in infile1map :
            print(*recl,sep="\t",file=ofo)
            print(*infile1map[seql[0]+seql[1]],sep="\t",file=ofo)
        elif ("".join(compflankseq)).upper() in infile1map:
            print(*recl,"recomp",sep="\t",file=ofo)
            print(*infile1map[("".join(compflankseq)).upper()],sep="\t",file=ofo)
        else:
            infile1map[seql[0]+seql[1]]=recl
    f.seek(0);ofo.close()

    print("notice! when you use grep -w/vFf to get the same/different file, notice some ID may not only occur in the ID column,but also duplicates cols ")
    print("some times may print one records multiple times")
     
    print("##############################")    
    print("remove dup >>>>>>>>>>")
        
    dupf=open(input+"dup",'r');remdupf=open(input+"redup",'w');print(f.readline().strip(),file=remdupf)
    for line in dupf:
        recl=re.split(r"\t",line.strip())
        if recl[1] in dupmap:
            dupmap[recl[1]].append(recl)
        else:
            dupids.append(recl[1])
            dupmap[recl[1]]=[recl]
    print(len(dupids))
    for line in f:
        recl=re.split(r"\t",line.strip())
        if  recl[1] not in dupids:
            print(*recl,sep="\t",file=remdupf)

    f.close();remdupf.close();dupf.close()

