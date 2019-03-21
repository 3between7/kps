# -*- coding: UTF-8 -*-
'''
update on 2019-02-18
@author: liurui

'''

import copy,re,sys,os,copy,pickle,configparser


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
                    a=re.sub('[!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~]+',"", linelist[linelist.index(chrsignal)+1])#for example chromosome 1
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
    Converting Roman numerals into Arabic numerals  
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
    
def getRefSeqBypos_faster(refFastahandle, fasterrefindex, currentChromNO, startpos, endpos, currentChromNOlen=None, seektuple=()):
    '''
    pos start at 1
    seektuple=(filepos,basesbeforefilepos)
    the refSeqMap has only one chromosome's sequence
    There is no restriction on refFastahander
    refindex must indexed by generateFasterRefIndex func
    return {'1': [0, 'A']}
    '''    
    refSeqMap = {}
    if startpos <= 0:
        startpos = 1
    if currentChromNOlen != None and endpos > currentChromNOlen:
        endpos = currentChromNOlen

    filehandle = refFastahandle
    if not seektuple or seektuple[1] > startpos:
        refSeqMap[currentChromNO] = [startpos - 1]
        low=0;high=len(fasterrefindex[currentChromNO])-1
        while low<=high:
            mid=(low + high)>>1
            if fasterrefindex[currentChromNO][mid][0]<startpos:
                low=mid+1
            elif fasterrefindex[currentChromNO][mid][0]>startpos:
                high=mid-1
            else:
                perivouspos_idx=mid
#                 filehandle.seek(refindex[currentChromNO][mid][1])
                break
        else:
            if fasterrefindex[currentChromNO][high][0]>startpos:
                high-=1
                perivouspos_idx=high
            else:
                perivouspos_idx=high
        filehandle.seek(fasterrefindex[currentChromNO][perivouspos_idx][1])
          # seekmap is empty so go to the first bases of the currentChromNO
        if fasterrefindex[currentChromNO][perivouspos_idx][0]<startpos:
            preseq = filehandle.read(startpos- fasterrefindex[currentChromNO][perivouspos_idx][0])
            dn = preseq.count('\n')
            while dn != 0:
                preseq = filehandle.read(dn)
                dn = preseq.count('\n')
        elif fasterrefindex[currentChromNO][perivouspos_idx][0]==startpos:
            pass
        else:
            pass
            
        # now filehander is right stay at the startpos
        myseqline = filehandle.read(endpos - startpos + 1)
        myseqn = myseqline.count('\n')
#        if len(myseqline)>200:
#            print(myseqn)
#            exit(-1)
#        print("myseqline=",myseqline,"myseqn", myseqn)
        while myseqn != 0:  # fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehandle.read(myseqn)
            myseqn = myseqline.count('\n')
            
#            print(currentChromNO,myseqline, myseqn)
        if myseqline.count('>') >= 1:
            exit(-1)
        refSeqMap[currentChromNO].extend(list(myseqline))
    else:
        filehandle.seek(seektuple[0])  # seekmap is not empty
        refSeqMap[currentChromNO] = [startpos - 1]
        preseq = filehandle.read(startpos - seektuple[1] - 1)
        dn = preseq.count('\n')
        while dn != 0:
            preseq = filehandle.read(dn)
            dn = preseq.count('\n')
        # now filehander is right stay at the startpos
        myseqline = filehandle.read(endpos - startpos + 1)
        myseqn = myseqline.count('\n')
        while myseqn != 0:  # fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehandle.read(myseqn)
            myseqn = myseqline.count('\n')
        refSeqMap[currentChromNO].extend(list(myseqline))
    plus = myseqline.count('>')
    if plus != 0:
        return -1
    return refSeqMap 

