#-*- coding: UTF-8 -*-

import random,re,pickle,copy


class RiceCaculator():

	def __init__(self):
		self.n=0
		self.id=[]
		self.tiling=[]
		self.minMAFofPop=[]
		self.indexID=None

	def process(self,T):
		self.n+=1
		self.id.append(T[1])
		self.tiling.append(T[2])
		self.minMAFofPop.append(T[3])

	def getMAF(self,l):
		maflist=[]
		for each in l:
			maflist.append(self.minMAFofPop[each])
		maxMAF=max(maflist)
		indexMAF=maflist.index(maxMAF)
		self.indexID=l[indexMAF]
		print("getMAF ",self.indexID)
		      
	def getResult(self):
		print(self.id)
		print(self.tiling)
		print(self.minMAFofPop)
		print("old :",self.indexID)
		order1=[];order6=[]
		order2=[];order7=[]
		order3=[]
		order4=[]
		order5=[]
		value=[]
		totalNums=self.n
		N=-1
		for each in self.tiling:
			N+=1
			if each == "1.0":
				order1.append(N)
			elif each == "2.0":
				order2.append(N)
			elif each == "3.0":
				order3.append(N)
			elif each == "4.0":
				order4.append(N)
			elif each == "5.0":
				order5.append(N)
			elif each == "6.0":
				order6.append(N)
			else:
				order7.append(N)
#		print(order1);print(order2);print(order3);print(order4);print(order5)
		try:
			if order1:
				print(self.indexID)
				n=self.getMAF(order1)
				value.append(self.id[self.indexID])
				print("new",self.indexID)
			elif order2:
				print(self.indexID)
				n=self.getMAF(order2)
				value.append(self.id[self.indexID])
				print("new",self.indexID)
			elif order3:
				print(self.indexID)
				n=self.getMAF(order3)
				value.append(self.id[self.indexID])
				print("new",self.indexID)
			elif order4:
				print(self.indexID)
				n=self.getMAF(order4)
				value.append(self.id[self.indexID])
				print("new",self.indexID)
			elif order5:
				print(self.indexID)
				n=self.getMAF(order5)
				value.append(self.id[self.indexID])
				print("new",self.indexID)
			elif order6:
				print(self.indexID)
				n=self.getMAF(order6)
				value.append(self.id[self.indexID])
				print("new",self.indexID)
			else:
				print(self.indexID)
				n=self.getMAF(order5)
				value.append(self.id[self.indexID])
				print("new",self.indexID)
		except:
			print("空窗口！")
		self.id=[]
		self.minMAFofPop=[]
		self.tiling=[]
		self.n=0
 
		return totalNums,value
				
		
class Window():
    def __init__(self):
        super().__init__()
        self.winValueL = []  # [(startPos,lastPos,value),(),,,,,,]

    def slidWindowOverlap(self, L, L_End_Pos, windowWidth, slideSize, Caculator,L_Start_Pos=0):
        print("L_End_Pos",L_End_Pos,L_Start_Pos)
        """
        require at least one SNP pos in L in sliding window
        window slide from L_Start_Pos to L_End_Pos
        L = [(pos,p1,p2,p3,A_base_idx),(pos,"a,b","c,d","e,f",0),(pos,"a,b","c,d","e,f",1),....] for D-statistics wihtout "no covered"
        or 
        L = [(pos, REF, ALT, INFO,FORMAT,sampleslist),(pos, REF, ALT, INFO,FORMAT,sampleslist),(),...........] for any score need one vcf,
        or 
        L = [(pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist)),(pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist)),(),...........] for any score need two or more vcf,for example two vcf's compare,eg. fst, one or multiple vcf caculate hp
        like the two situation upside,return a value
        or
        L = [(pos,samples1dp,samples2dp,samples3dp,,,),(pos,samples1dp,samples2dp,samples3dp,,,),(),(),......]#in this situation ,value formation like this ([sample1_pecentage,sample2_pecentage,,,],[sample1_average_depth,sample2_average_depth,,,])
        
        """
#         del self.winValueL[:]
        self.winValueL = []  # notice here
        nextIdx = -1  # always be -1 if windowWidth == slideSize
        currentIdx = 0
        winStart = L_Start_Pos
        FoundNextIdx = False
        firstComeInWin = True
        notjustforsnp = True
        for findfirstidx in range(len(L)):
            if L[findfirstidx][0]>winStart:
                currentIdx=findfirstidx
                break
        while currentIdx != len(L) and L[currentIdx][0]<= L_End_Pos:
#             print(L[currentIdx][0],L[currentIdx])
            if L[currentIdx][0] > winStart and  L[currentIdx][0] <= (winStart + windowWidth):
#                if notjustforsnp or (len(L[currentIdx][1])==1 and re.search(r'[^a-zA-Z]', L[currentIdx][2]) != None and len(L[currentIdx][2])==1):# it's not a snp? indel or cnv
                if firstComeInWin:
                    startPos = L[currentIdx][0]
                    firstComeInWin = False
                lastPos = L[currentIdx][0]
                Caculator.process(L[currentIdx])
                if FoundNextIdx == False and L[currentIdx][0] > (winStart + slideSize):  # always go to |currentIdx+=1|
                    nextIdx = currentIdx
                    FoundNextIdx = True
            else:
                noofsnps, value = Caculator.getResult()
                try:
                    self.winValueL.append((startPos, lastPos, noofsnps, value))
#                     print(startPos, lastPos, noofsnps, value)
                except:
                    print("no snp in first severl wins", len(L), currentIdx, value, L[currentIdx])
                    self.winValueL.append((0, 0, noofsnps, value))
                    winStart += slideSize
                    continue
#                 self.winValueL.append((0, 0, value))
                winStart += slideSize
                firstComeInWin = True
                
                FoundNextIdx = False
                if nextIdx == -1:
                    if slideSize >= windowWidth:
                        while not (L[currentIdx][0] > winStart and  L[currentIdx][0] <= (winStart + windowWidth)) and L[currentIdx][0] > winStart + windowWidth:
                            winStart += slideSize
                            noofsnps, value = Caculator.getResult()
                            self.winValueL.append((0, 0, noofsnps, value))
                        if L[currentIdx][0] < winStart:
                            while currentIdx != len(L):
                                if L[currentIdx][0] > winStart and L[currentIdx][0] <= (winStart + windowWidth):
                                    break
                                elif L[currentIdx][0] < winStart:
                                    winStart += slideSize
                                    noofsnps, value = Caculator.getResult()
                                    self.winValueL.append((0, 0, noofsnps, value))
                                currentIdx += 1
#                             self.winValueL.append((0,0,'NA'))
#                             winStart += slideSize
                    continue  # go to |if L[currentIdx][0] > winStart and L[currentIdx][0] < (winStart + windowWidth):| in upside block
                else:
                    currentIdx = nextIdx
                    nextIdx = -1
                    continue
                
            currentIdx += 1
        else:
            noofsnps, value = Caculator.getResult()
            try :
                self.winValueL.append((startPos, lastPos, noofsnps, value))
#                 print(startPos, lastPos, noofsnps, value)
            except UnboundLocalError:
                self.winValueL.append((0, 0, noofsnps, value))
#             if nextIdx!=-1:
#                 currentIdx = nextIdx
#                 nextIdx = -1
#                 while currentIdx != len(L):
#                     lastPos = L[currentIdx][0]
#                     Caculator.process(L[currentIdx])
#                     currentIdx += 1
#                 else:
#                     noofsnps, value = Caculator.getResult()
#                     try:
#                         self.winValueL.append((startPos, lastPos, noofsnps, value))
#                     except:
#                         self.winValueL.append((0, 0, noofsnps, value))
#            
        
        n = int((L_End_Pos-L_Start_Pos - (len(self.winValueL) * slideSize + windowWidth)) / slideSize) + 1
        for i in range(n):
            noofsnps, value = Caculator.getResult()
            self.winValueL.append((0, 0, noofsnps, value))
		
class Frq_Data():

	def __init__(self, freqfileName):
		super().__init__()
		self.FrqMap_AllChrom = {}
		self.FrqIndexMap = {}
		self.chromOrder = []
		self.frqfileName=freqfileName
		self.NumOfRecbychromOrder = []
		try:
			self.FrqIndexMap = pickle.load(open(freqfileName + ".myindex", 'rb'))
		except:
			Frq_Data.indexFRQ(frqName=freqfileName, indexFileName=(freqfileName + ".myindex"))
			self.FrqIndexMap = pickle.load(open(freqfileName + ".myindex", 'rb'))
		self.chromOrder = self.FrqIndexMap["chromOrder"]
		self.NumOfRecbychromOrder = self.FrqIndexMap["NumOfRecbychromOrder"]
		
	@staticmethod
	def indexFRQ(frqName, indexFileName):
		"""
		{chrom:position_in_file_of_first_SNP_of_this_chrom,chrom:position,,,,,,}
		"""
		frqfile = open(frqName, 'r')
		frqChromIndex = {}
		chromOrder = []
		NumOfRecbychromOrder = []
		line = frqfile.readline()   
		if re.search(r'id', line) != None:
			frqChromIndex["title"] = re.split(r'\t', line.strip()) 
		else:                                                      
			print("need title'#Organism	id    SEQ    Tiling Order	importance	CHR	POS   minMAFofPop'")
			exit(-1)       
		currentChrom = "temptodele"
		lastPosition = frqfile.tell()
		lastChromend_currentChromstartPostion = lastPosition
		print(line) 
		line = frqfile.readline()
		print("first line:",line)
		i = 0
		while line:   
			linelist = re.split(r"\t", line) 
			if currentChrom != linelist[5]:
				chromOrder.append(currentChrom)
				NumOfRecbychromOrder.append(i);i = 0  # collect the  number of snp recs  of the last chrom
				frqChromIndex[currentChrom] = (lastChromend_currentChromstartPostion, lastPosition)
				lastChromend_currentChromstartPostion = lastPosition
				currentChrom = linelist[5].strip()
			lastPosition = frqfile.tell()
			line = frqfile.readline()
			i += 1
		else:
			chromOrder.append(currentChrom.strip())
			NumOfRecbychromOrder.append(i - 1)  # collect the  number of snp recs  of the lastest chrom of all chroms
			frqChromIndex[currentChrom] = (lastChromend_currentChromstartPostion, lastPosition)
		frqChromIndex.pop("temptodele")
		i = chromOrder.index("temptodele")
		if i != 0:
			print("wrong indexFRQ")
			exit(-1)
		b = chromOrder.pop(i)
		a = NumOfRecbychromOrder.pop(i)
		frqChromIndex["chromOrder"] = chromOrder
		frqChromIndex["NumOfRecbychromOrder"] = NumOfRecbychromOrder
		pickle.dump(frqChromIndex, open(indexFileName, 'wb'))
		frqfile.close()

	def getFrqListByChrom(self, chrom,startpos=1,endpos=9999999999999999999999999999999999999999999999999999):
		"""
			although dilute and dilutetodensity can exist at the same time,but it not make sense and may final produce a bug.
			return a list that contain all vcf record of a chrom
		"""
		print("getFrqListByChrom",chrom,startpos,endpos)

		FrqList_A_Chrom = []
		if (chrom not in self.chromOrder) or startpos>=endpos:
			return []
			print(chrom + "didn't find in " + self.frqfileName)
		i = self.chromOrder.index(chrom.strip())
		frqfile = open(self.frqfileName, 'r')                
		frqfile.seek(self.FrqIndexMap[chrom][0])#        #find the first line
		filepos=frqfile.tell()
		line=frqfile.readline()
		while line:
#         find the startpos
			linelist = re.split(r'\t', line.strip())
			print(linelist)
			c_chrom = linelist[5].strip() #zhanghuanhuan
			pos=int(linelist[6].strip()) #zhanghuanhuan
			if chrom.strip()!=c_chrom:
				return FrqList_A_Chrom            
			if pos >= startpos : 
				break
			filepos=frqfile.tell() 
			line=frqfile.readline()
		else:
			if (linelist is  None) or (pos < startpos):
				return FrqList_A_Chrom
        #from startpos to collected recs    
		frqfile.seek(filepos)
		linescontent=frqfile.read(self.FrqIndexMap[chrom][1]-filepos)
		print("need check encoding is utf-8")
		frqlineslist=re.split(r"\n",linescontent.strip())
		for line in frqlineslist:
			linelist = re.split(r'\t', line.strip())
			c_chrom = linelist[5].strip()
			pos = int(linelist[6].strip())
			if pos>endpos or chrom!=c_chrom:
				break
			id=linelist[1].strip()
			TilingOrder = linelist[3].strip()
			minMAFofPop = float(linelist[7].strip())
			FrqList_A_Chrom.append((pos,id,TilingOrder,minMAFofPop))
		frqfile.close()
		print("getFrqListByChrom",chrom,len(FrqList_A_Chrom),self.FrqIndexMap[chrom],self.NumOfRecbychromOrder[i],"total recs in this file belong to this chrom") 
		return FrqList_A_Chrom



if __name__ == "__main__":
	
	chrlenm=[158337067,137060424,121430405,120829699,121191424,119458736,112638659,113384836,105708250,104305016,107310763,91163125,84240350,84648390,85296676,81724687,75158596,66004023,64057457,72042655,71599096,61435874,52530062,62714930,42904170,51681464,45407902,46312546,51505224,148823899]

	def pick(num,length):
		MyData=Frq_Data("/media/jason/Seagate Backup Plus Drive/Cattle_DataSet1_update.txt")
		ListName="FrqList_"+num+"_Chrom"
		ListName=MyData.getFrqListByChrom(num)
		Mywindow=Window()
		c=RiceCaculator()
		Mywindow.slidWindowOverlap(ListName, length, 1100, 1100,c)
		ss=copy.deepcopy(Mywindow.winValueL)
		ChromName="Chrom"+num+"id"
		ChromName=[]
		for eachwindow in ss:
			try:
				if eachwindow[3]!="NA":
					ChromName.append(eachwindow[3])
			except:
				continue
		pickle.dump(ChromName, open("Cattle_Chr"+num+"picked.pickle", 'wb'))

	for i in range(1,31):		
		pick(str(i),chrlenm[i-1])


	

 




