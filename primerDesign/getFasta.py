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



if __name__ == '__main__':

	'''
	dictC:{'1':[pos1,pos2,...posN],'2':[pos1,pos2,...posN],'i':[,,,,,]..}

	'''
	
	f1=open("/home/zhanghuanhuan/GaoHong/gao_closeSites.txt") # the input file
	f2=open("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/6-formatDataSet1/Bos_taurus.UMD3.1.dna.chromosomes.fa")
	f3=open("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/6-formatDataSet1/CattleChrLength.txt")
	f4=open("/home/zhanghuanhuan/GaoHong/gao_closeSites.seq","w") # the output file
	refidx=pickle.load(open("/lustre/home/zhanghuanhuan/Data/ChipDesign/Cattle/6-formatDataSet1/Bos_taurus.UMD3.1.dna.chromosomes.myfasteridx",'rb'))
	
	# get the length of chromsomes
	lengthdict={}
	for eachline in f3:
		linelist1=re.split(r"\t",eachline.strip())
		lengthdict[linelist1[0]]=linelist1[1]

	# get dict of pos of chromsomes
	listX=[]
	chrom=0
	dictC={}
	for eachline in f1:
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
	print(dictC)
	# get fasta sequence of the pos
	for c in dictC.keys():
		listN=dictC[c] # the POS of the chronsome,the format is:[pos1,pos2,...posN]
		n=len(listN)  # the nums of sites of the chromsome
		listPos=[listN[0]]
		for i in range(n):
			if i+1==n:
				listrt=[]
				leftPos=int(listPos[0])-400
				rightPos=int(listPos[-1])+400
				chr=c
				N=leftPos-1
				for p in listPos:
					listrt.append(str(int(p)-N))
				print("chromsome:",chr,"leftPos:",leftPos,"rightPos:",rightPos,"N:",N)
				RefSeqMap = getRefSeqBypos_faster(f2,refidx,chr,leftPos,rightPos,int(lengthdict[chr]))
				f4.write(">Chr"+chr+";absolute_pos:"+",".join(listPos)+";relative_pos:"+",".join(listrt)+"\n")
				f4.write("".join(RefSeqMap[chr][1:]).upper()+"\n")
			elif int(listN[i+1])-int(listN[i])<=400:
				listPos.append(listN[i+1])
			else:
#				print(listPos)
				listrt=[]
				leftPos=int(listPos[0])-400
				rightPos=int(listPos[-1])+400
				chr=c
				N=leftPos-1
				print("chromsome:",chr,"leftPos:",leftPos,"rightPos:",rightPos,"N:",N)
				RefSeqMap = getRefSeqBypos_faster(f2,refidx,chr,leftPos,rightPos,int(lengthdict[chr]))
				#print(RefSeqMap)
				# write the title of fasta,the format is : >Chr1;absolute_pos:pos1,pos2,pos3;relative_pos:p1,p2,p3
				for p in listPos:
					listrt.append(str(int(p)-N))
#				print(listrt)
				print("===========================================")
				f4.write(">Chr"+chr+";absolute_pos:"+",".join(listPos)+";relative_pos:"+",".join(listrt)+"\n")
				f4.write("".join(RefSeqMap[chr][1:]).upper()+"\n")
				listPos=[listN[i+1]]

