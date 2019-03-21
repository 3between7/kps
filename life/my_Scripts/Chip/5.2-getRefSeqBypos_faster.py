

import re,sys,os,copy,pickle

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


#if __name__=="__main__":

	
def getSeq(file1,file2):
	f=open("Bos_taurus.UMD3.1.dna.chromosomes.fa")
	refidx=pickle.load(open("Bos_taurus.UMD3.1.dna.chromosomes.myfasteridx",'rb'))
	f1=open(file1)
	f2=open(file2,"w")
	f3=open("CattleChrLength.txt")
	lengthdict={}
	for eachline in f3:
		linelist1=re.split(r"\t",eachline.strip())
		lengthdict[linelist1[0]]=linelist1[1]
	for line in f1:
		ll2=re.split(r"\s+",line.strip())
		RefSeqMap = getRefSeqBypos_faster(f,refidx,ll2[0],int(ll2[1])-35,int(ll2[1])+35,int(lengthdict[ll2[0]]))
		f2.write(ll2[0]+"\t"+ll2[1]+"\t"+"".join(RefSeqMap[ll2[0]][36]).upper()+"\n")		
	f.close()
	f1.close()
	f2.close()
	f3.close()

#chrom:['1', '2',...,'30']

getSeq("../9-getRecord/CattlePickedRecord_v1.ChrPos","../10-getInterferingSNP/CattleSitesREF.txt")	


    


    
    
    
    
    
    
    

