# -*- coding: UTF-8 -*-
import sys
import re
import pickle
import os
import configparser
import pandas as pd
from pandas.core.frame import DataFrame
from collections import defaultdict

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
                    a=re.sub('[’!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~]+',"", linelist[linelist.index(chrsignal)+1])#for example chromosome 1,
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
    Convert Roman numerals to Arabic numerals 
    s.remove
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
        while myseqn != 0:  # fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehandle.read(myseqn)
            myseqn = myseqline.count('\n')
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

if __name__ == '__main__':
    
    print("Please check the format of Reference genome sequence file")
    fasta=input("Please input the path of Reference genome sequence :")
    ChrLength=input("Please input the path of the file of chromosome's length :")
    sites=input("Please input the path of sites file(chr and pos) :") # the input file
    seqlength=input("Please input the length of Flanking sequence of sites:")
#   os.system("sort -t $'\t' -k 1n -k 2n +sites -o +sites ") #notice: the sites need to be sorted!
#    output=input("Please input the path of sequences file:")
    try:
        refidx=pickle.load(open(fasta+".myfasteridx",'rb'))
        print("Sucessfully load "+fasta+".myfasteridx")
    except:
        generateFasterRefIndex(fasta,fasta+".myfasteridx")
        print("The myfasteridx file of Reference genome sequence was generated!")

    ##################################################################
    f=open(fasta)
    refidx=pickle.load(open(fasta+".myfasteridx",'rb'))
    f1=open(ChrLength)
    f2=open(sites)
    f3=open("./sequences_normalSites.txt","w")
    f4=open("./sequence_closeSites.txt","w")
    
    # get the length of chromsomes
    lengthdict={}
    for eachline in f1:
        ll=re.split(r"\t",eachline.strip())
        lengthdict[ll[0]]=ll[1]
    #################################################################

    # get dict of pos of chromsomes
    listX=[]
    chrom=0
    dictC={}
    for eachline in f2: #chr pos
        ll=re.split("\t",eachline.strip())
        chromC=ll[0]
        pos=ll[1]
        if chromC !=chrom:
            if listX:
                dictC[chrom]=listX
            chrom=chromC
            listX=[pos]
        else:
            listX.append(pos)
    dictC[chrom]=listX
    # get fasta sequence of the pos
    for c in dictC.keys():
        listN=dictC[c] # the POS of the chronsome,the format is:[pos1,pos2,...posN]
        n=len(listN)  # the nums of sites of the chromsome
        listPos=[listN[0]]
        for i in range(n):
            if i+1==n:
                listrt=[]
                leftPos=int(listPos[0])-int(seqlength)
                rightPos=int(listPos[-1])+int(seqlength)              
                chr=c
                N=leftPos-1
                RefSeqMap = getRefSeqBypos_faster(f,refidx,chr,leftPos,rightPos,int(lengthdict[chr]))
                for p in listPos:
                    listrt.append(str(int(p)-N))
                if len(listPos)==1:
                    f3.write(">Chr"+chr+";absolute_pos:"+",".join(listPos)+";relative_pos:"+",".join(listrt)+"\n")
                    f3.write("".join(RefSeqMap[chr][1:]).upper()+"\n")
                else:
                    f4.write(">Chr"+chr+";absolute_pos:"+",".join(listPos)+";relative_pos:"+",".join(listrt)+"\n")
                    f4.write("".join(RefSeqMap[chr][1:]).upper()+"\n")

            elif int(listN[i+1])-int(listN[i])<=int(seqlength)*2:
                listPos.append(listN[i+1])
            else:
                listrt=[]
                leftPos=int(listPos[0])-int(seqlength)
                rightPos=int(listPos[-1])+int(seqlength)
                chr=c
                N=leftPos-1
                RefSeqMap = getRefSeqBypos_faster(f,refidx,chr,leftPos,rightPos,int(lengthdict[chr]))
                # write the title of fasta,the format is : >Chr1;absolute_pos:pos1,pos2,pos3;relative_pos:p1,p2,p3
                for p in listPos:
                    listrt.append(str(int(p)-N))
                if len(listPos)==1:
                    f3.write(">Chr"+chr+";absolute_pos:"+",".join(listPos)+";relative_pos:"+",".join(listrt)+"\n")
                    f3.write("".join(RefSeqMap[chr][1:]).upper()+"\n")
                else:
                    f4.write(">Chr"+chr+";absolute_pos:"+",".join(listPos)+";relative_pos:"+",".join(listrt)+"\n")
                    f4.write("".join(RefSeqMap[chr][1:]).upper()+"\n")
                listPos=[listN[i+1]]
    f.close();f1.close();f2.close();f3.close();f4.close()
    ##################################################################

    f3 = open("./sequences_normalSites.txt","r")
    f4 = open("primer3_settingfile.txt","w")
    fw = open("primer3_inputfile.txt","w")
    
    ## the primer3 inputfile
    ##########################################################################################
    print("Please set up Primers_core!")
    primersfile=input("Please set the outputfile :") # the output file
    SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=input("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST:")
    PRIMER_PRODUCT_SIZE_RANGE=input("PRIMER_PRODUCT_SIZE_RANGE:")
    PRIMER_NUM_RETURN=input("PRIMER_NUM_RETURN:")
    PRIMER_MIN_SIZE=input("PRIMER_MIN_SIZE:")
    PRIMER_OPT_SIZE=input("PRIMER_OPT_SIZE:")
    PRIMER_MAX_SIZE=input("PRIMER_MAX_SIZE:")
    PRIMER_MIN_TM=input("PRIMER_MIN_TM:")
    PRIMER_OPT_TM=input("PRIMER_OPT_TM:")
    PRIMER_MAX_TM=input("PRIMER_MAX_TM:")
    PRIMER_PAIR_MAX_DIFF_TM=input("PRIMER_PAIR_MAX_DIFF_TM:")
    PRIMER_MAX_END_GC=input("PRIMER_MAX_END_GC:")
    PRIMER_MIN_GC=input("PRIMER_MIN_GC:")
    PRIMER_MAX_GC=input("PRIMER_MAX_GC:")
    PRIMER_THERMODYNAMIC_PARAMETERS_PATH=input("PRIMER_THERMODYNAMIC_PARAMETERS_PATH:")
    ########################################################################################
    
    ## Primer3_settingfile
    f4.write("Primer3 File - http://primer3.org"+"\n")
    f4.write("P3_FILE_TYPE=settings"+"\n")
    f4.write("\n")
    f4.write("PRIMER_FIRST_BASE_INDEX=1"+"\n")
    f4.write("PRIMER_TASK=generic"+"\n")
    f4.write("P3_FILE_ID=Settings for PCR amplification followed by Sanger sequencing on both strands using PCR primers"+"\n")
    f4.write("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=" + SEQUENCE_PRIMER_PAIR_OK_REGION_LIST +"\n")
    f4.write("PRIMER_PRODUCT_SIZE_RANGE=" + PRIMER_PRODUCT_SIZE_RANGE +"\n")
    f4.write("PRIMER_NUM_RETURN=" + PRIMER_NUM_RETURN +"\n")
    f4.write("PRIMER_MIN_SIZE=" + PRIMER_MIN_SIZE +"\n")
    f4.write("PRIMER_OPT_SIZE=" + PRIMER_OPT_SIZE +"\n")
    f4.write("PRIMER_MAX_SIZE=" + PRIMER_MAX_SIZE +"\n")
    f4.write("PRIMER_MIN_TM=" + PRIMER_MIN_TM +"\n")
    f4.write("PRIMER_OPT_TM=" + PRIMER_OPT_TM +"\n")
    f4.write("PRIMER_MAX_TM=" + PRIMER_MAX_TM +"\n")
    f4.write("PRIMER_PAIR_MAX_DIFF_TM=" + PRIMER_PAIR_MAX_DIFF_TM +"\n")
    f4.write("PRIMER_MAX_END_GC=" + PRIMER_MAX_END_GC +"\n")
    f4.write("PRIMER_MIN_GC=" + PRIMER_MIN_GC +"\n")
    f4.write("PRIMER_MAX_GC="+ PRIMER_MAX_GC +"\n")
    ss=("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=",PRIMER_THERMODYNAMIC_PARAMETERS_PATH,"\n")
    path="".join(ss)
    #f4.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$PATH"+"\n")
    f4.write(path)
    f4.write("PRIMER_MIN_THREE_PRIME_DISTANCE=3 "+"\n")
    f4.write("PRIMER_MAX_LIBRARY_MISPRIMING=12.00"+"\n")
    f4.write("PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00"+"\n")
    f4.write("PRIMER_MAX_END_STABILITY=9.0"+"\n")
    f4.write("PRIMER_MAX_SELF_ANY_TH=45.00"+"\n")
    f4.write("PRIMER_MAX_SELF_END_TH=35.00"+"\n")
    f4.write("PRIMER_PAIR_MAX_COMPL_ANY_TH=45.00"+"\n")
    f4.write("PRIMER_PAIR_MAX_COMPL_END_TH=35.00"+"\n")
    f4.write("PRIMER_MAX_HAIRPIN_TH=24.00"+"\n")
    f4.write("PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00"+"\n")
    f4.write("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00"+"\n")
    f4.write("PRIMER_TM_FORMULA=1"+"\n")
    f4.write("PRIMER_SALT_CORRECTIONS=1"+"\n")
    f4.write("PRIMER_SALT_MONOVALENT=50.0"+"\n")
    f4.write("PRIMER_INTERNAL_SALT_MONOVALENT=50.0"+"\n")
    f4.write("PRIMER_SALT_DIVALENT=1.5"+"\n")
    f4.write("PRIMER_INTERNAL_SALT_DIVALENT=1.5"+"\n")
    f4.write("PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1"+"\n")
    f4.write("PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=1"+"\n")
    f4.write("PRIMER_LOWERCASE_MASKING=0"+"\n")
    f4.write("yPRIMER_MAX_NS_ACCEPTED=0"+"\n")
    f4.write("PRIMER_MAX_POLY_X=4"+"\n")
    f4.write("PRIMER_OUTSIDE_PENALTY=0"+"\n")
    f4.write("PRIMER_GC_CLAMP=0"+"\n")
    f4.write("PRIMER_LIBERAL_BASE=1"+"\n")
    f4.write("PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0"+"\n")
    f4.write("PRIMER_PICK_ANYWAY=1"+"\n")
    f4.write("PRIMER_WT_TM_LT=1.0"+"\n")
    f4.write("PRIMER_WT_TM_GT=1.0"+"\n")
    f4.write("PRIMER_WT_SIZE_LT=1.0"+"\n")
    f4.write("PRIMER_WT_SIZE_GT=1.0"+"\n")
    f4.write("PRIMER_WT_GC_PERCENT_LT=0.0"+"\n")
    f4.write("PRIMER_WT_GC_PERCENT_GT=0.0"+"\n")
    f4.write("PRIMER_WT_SELF_ANY_TH=0.0"+"\n")
    f4.write("PRIMER_WT_SELF_END_TH=0.0"+"\n")
    f4.write("PRIMER_WT_HAIRPIN_TH=0.0"+"\n")
    f4.write("PRIMER_WT_NUM_NS=0.0"+"\n")
    f4.write("PRIMER_WT_LIBRARY_MISPRIMING=0.0"+"\n")
    f4.write("PRIMER_WT_SEQ_QUAL=0.0"+"\n")
    f4.write("PRIMER_WT_END_QUAL=0.0"+"\n")
    f4.write("PRIMER_WT_POS_PENALTY=0.0"+"\n")
    f4.write("PRIMER_WT_END_STABILITY=0.0"+"\n")
    f4.write("PRIMER_WT_TEMPLATE_MISPRIMING_TH=0.0"+"\n")
    f4.write("PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0"+"\n")
    f4.write("PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0"+"\n")
    f4.write("PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0"+"\n")
    f4.write("="+"\n")
    f4.close()
    
    ## Primer3_inputfile
    dictSeq = {}
    for eachline in f3:
        if ">" in eachline:
            key = eachline.strip()[1:]
        else:
            value = eachline.strip()
            dictSeq[key] = value
    for item in dictSeq.keys():
        fw.write("SEQUENCE_ID="+item+"\n")
        fw.write("SEQUENCE_TEMPLATE="+dictSeq[item]+"\n")
        fw.write("="+"\n")
    f3.close();fw.close()
    
    ## run Primer3
    os.system("primer3_core < primer3_inputfile.txt --format_output -p3_settings_file=primer3_settingfile.txt -output=./primer3_outfile.txt")
    
    ## extract Primers
    re_tag = re.compile(r'^PRIMER PICKING RESULTS FOR (\w+;\w+:\w+;\w+:\w+)')
    re_left = re.compile(r'LEFT PRIMER')
    re_right = re.compile(r'RIGHT PRIMER')
    
    infile = 'primer3_outfile.txt'

    tags = []
    left = None
    right = None
    
    primers = defaultdict(list)
    
    with open(infile, 'r') as f,open(primersfile,'w') as out:
        for line in f:
            if re_tag.search(line):
                tags = re_tag.search(line).groups()[0].split(',')
                left = None
                right = None
                continue
            
            if re_left.findall(line):
                left = line.rstrip().split()[-1]
                continue
    
            if re_right.findall(line):
                right = line.rstrip().split()[-1]
                for i in tags:
                    primers[i].append((left, right))
                left = None
                right = None
                continue
        k=list(primers.keys())
        v=list(primers.values())
        result=pd.DataFrame(list(zip(k,v)),columns=['k','v'])
        result.to_csv(out,sep=',', header=False, index=False)
    
    ## remove the  intermediate files
    os.remove("./primer3_inputfile.txt")
    os.remove("./primer3_settingfile.txt")
    os.remove("./primer3_outfile.txt")
