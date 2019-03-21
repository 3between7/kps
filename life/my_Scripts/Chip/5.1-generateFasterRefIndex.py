# -*- coding: UTF-8 -*-
import copy
import re
import pickle,os,configparser


def generateFasterRefIndex(refFastaFileName, indexFileName,mapname=None,startchar=">",chrsignal=None,romanSignal=False):
    refFastaFile = open(refFastaFileName, 'r')
    refChromIndex = {}
    refline = refFastaFile.readline()
    while refline:
        if re.search(r'^['+startchar+']', refline) != None:
            basecount=1
            m=1
            if mapname == "transcript:":
                currentChromNo = re.search(r'transcript:(.*?)\s+', refline).group(1).strip()
            else:
                if not chrsignal or chrsignal not in refline:
                    a = re.search(r'^'+startchar+'([^'+startchar+'|]+)', (re.split(r'\s+', refline))[0]).group(1)
                else:
                    linelist=re.split(r'\s+', refline)
                    a=re.sub('[��!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~]+',"", linelist[linelist.index(chrsignal)+1])#for example chromosome 1,
                a=transform_roman_num2_alabo(a,romanSignal)
                currentChromNo=a.replace("chr","").replace("CHR", "")
                print(currentChromNo,type(currentChromNo))
            refChromIndex[currentChromNo] = [(basecount,int(refFastaFile.tell()))]# (no of base befor,cur file pos)
        else:
            basecount+=len(refline.strip())
            if basecount>=6000*m:
                refChromIndex[currentChromNo].append((basecount,int(refFastaFile.tell())))
                m+=1
        refline = refFastaFile.readline()
    pickle.dump(refChromIndex,open(indexFileName, 'wb'))
    refFastaFile.close()
    
def transform_roman_num2_alabo(one_str,changesignal=True):  
    ''''' 
    将罗马数字转化为阿拉伯数字 
    '''
    if re.search('^M{0,4}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})$',one_str)!=None and changesignal:  
        define_dict={'I':1,'V':5,'X':10,'L':50,'C':100,'D':500,'M':1000}  
        if one_str=='0':  
            return 0  
        else:  
            res=0  
            for i in range(0,len(one_str)):  
                if i==0 or define_dict[one_str[i]]<=define_dict[one_str[i-1]]:  
                    res+=define_dict[one_str[i]]  
                else:  
                    res+=define_dict[one_str[i]]-2*define_dict[one_str[i-1]]  
            return str(res)
    else:
        return one_str

generateFasterRefIndex("Bos_taurus.UMD3.1.dna.chromosomes.fa","Bos_taurus.UMD3.1.dna.chromosomes.myfasteridx")
