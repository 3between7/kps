# -*- coding: UTF-8 -*-
'''
update on 2019-03-07
@author: liurui,zhanghuanhuan

'''

import random,re,pickle,copy
import pandas as pd

class Window():
    def __init__(self):
        super().__init__()
        self.winValueL = []  # [(startPos,lastPos,value),(),,,,,,]

    def slidWindowOverlap(self, L, L_End_Pos, windowWidth, slideSize,Caculator,L_Start_Pos=0):
        print("L_End_Pos",L_End_Pos,L_Start_Pos)
        
        """
        require at least one SNP pos in L in sliding window
        window slide from L_Start_Pos to L_End_Pos
		L=[(pos,id,TilingOrder,importance,MAF),(pos,id,TilingOrder,importance,MAF),,,(pos,id,TilingOrder,importance,MAF)]
        """
        
        self.winValueL = []  # notice here
        nextIdx = -1  # always be -1 if windowWidth == slideSize
        currentIdx = 0
        winStart = L_Start_Pos
        FoundNextIdx = False
        firstComeInWin = True
        #notjustforsnp = True
        
        for findfirstidx in range(len(L)):
            if L[findfirstidx][0]>winStart:
                currentIdx=findfirstidx
                break
        while currentIdx != len(L) and L[currentIdx][0]<= L_End_Pos:

            if L[currentIdx][0] > winStart and  L[currentIdx][0] <= (winStart + windowWidth):

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
                    self.winValueL.append((startPos, lastPos, noofsnps, value)) #value is a list,[(,,,[]),(,,,[]),...]
                except:
                    print("no snp in first severl wins", len(L), currentIdx, value, L[currentIdx])
                    self.winValueL.append((0, 0, noofsnps, value))
                    winStart += slideSize
                    continue

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
            except UnboundLocalError:
                self.winValueL.append((0, 0, noofsnps, value))
      
        n = int((L_End_Pos-L_Start_Pos - (len(self.winValueL) * slideSize + windowWidth)) / slideSize) + 1
        for i in range(n):
            noofsnps, value = Caculator.getResult()
            self.winValueL.append((0, 0, noofsnps, value))
		
class prepareforwindow():
    
    def __init__(self, DataSet1):
        super().__init__()
        self.SitesMap_AllChrom = {}
        self.SitesIndexMap = {}
        self.chromOrder = []
        self.SitesfileName=DataSet1
        self.NumOfRecbychromOrder = []
        try:
            self.SitesIndexMap = pickle.load(open(DataSet1 + ".myindex", 'rb'))
        except:
            prepareforwindow.indexScoreFile(SitesName=DataSet1, indexFileName=(DataSet1 + ".myindex"))
            self.SitesIndexMap = pickle.load(open(DataSet1 + ".myindex", 'rb'))
        self.chromOrder = self.SitesIndexMap["chromOrder"]
        self.NumOfRecbychromOrder = self.SitesIndexMap["NumOfRecbychromOrder"]
	
    @staticmethod
    def indexScoreFile(SitesName, indexFileName):
        """
        return SitesChromIndex:
        {'title':the header of the dataset1,'chr':[lastChromend_currentChromstartPostion, lastPosition],
        "chromOrder":['1','2','3',...,'N'],'NumOfRecbychromOrder':[9302,111,23344,...,9999] }
        """
        Sitesfile = open(SitesName, 'r')
        SitesChromIndex = {}
        chromOrder = []
        NumOfRecbychromOrder = []
        line = Sitesfile.readline()   
        if re.search(r'snpid', line) != None:
            SitesChromIndex["title"] = re.split(r'\t', line.strip()) 
        else:                                                      
            print("need title'#Organism	snpid    seventyonemer    Tiling Order	importance	chromosome	POS	Ref Allele	Alt Allele	MAF'")
            exit(-1)       
        currentChrom = "temptodele"
        lastPosition = Sitesfile.tell()
        lastChromend_currentChromstartPostion = lastPosition
        print("title:",line) 
        line = Sitesfile.readline()
        print("first line:",line)
        i = 0
        while line:   
            linelist = re.split(r"\t", line) 
            if currentChrom != linelist[5]:
                chromOrder.append(currentChrom)
                NumOfRecbychromOrder.append(i);i = 0  # collect the  number of snp recs  of the last chrom
                SitesChromIndex[currentChrom] = (lastChromend_currentChromstartPostion, lastPosition)
                lastChromend_currentChromstartPostion = lastPosition
                currentChrom = linelist[5].strip()
            lastPosition = Sitesfile.tell()
            line = Sitesfile.readline()
            i += 1
        else:
            chromOrder.append(currentChrom.strip())
            NumOfRecbychromOrder.append(i - 1)  # collect the  number of snp recs  of the lastest chrom of all chroms
            SitesChromIndex[currentChrom] = (lastChromend_currentChromstartPostion, lastPosition)
        SitesChromIndex.pop("temptodele")
        i = chromOrder.index("temptodele")
        if i != 0:
            print("wrong indexScoreFile")
            exit(-1)
        b = chromOrder.pop(i)
        a = NumOfRecbychromOrder.pop(i)
        SitesChromIndex["chromOrder"] = chromOrder
        SitesChromIndex["NumOfRecbychromOrder"] = NumOfRecbychromOrder
        pickle.dump(SitesChromIndex, open(indexFileName, 'wb'))
        Sitesfile.close()

    def getSitesListByChrom(self, chrom,startpos=1,endpos=9999999999999999999999999999999999999999999999999999):
        
        """
        return:[(pos,id,TilingOrder,importance,MAF),(pos,id,TilingOrder,importance,MAF),,,(pos,id,TilingOrder,importance,MAF)]
        
        """
        print("getSitesListByChrom",chrom,startpos,endpos)
        SitesList_A_Chrom = []
        if (chrom not in self.chromOrder) or startpos>=endpos:
            return []
            print(chrom + "didn't find in " + self.frqfileName)
        i = self.chromOrder.index(chrom.strip())
        Sitesfile = open(self.SitesfileName, 'r')                
        Sitesfile.seek(self.SitesIndexMap[chrom][0])##find the first line
        filepos=Sitesfile.tell()
        line=Sitesfile.readline()
        while line:
#         find the startpos
            linelist = re.split(r'\t', line.strip())
            #print(linelist)
            c_chrom = linelist[5].strip() 
            pos=int(linelist[6].strip()) 
            if chrom.strip()!=c_chrom:
                return SitesList_A_Chrom            
            if pos >= startpos : 
                break
            filepos=Sitesfile.tell() 
            line=Sitesfile.readline()
        else:
            if (linelist is  None) or (pos < startpos):
                return SitesList_A_Chrom
        #from startpos to collected recs    
        Sitesfile.seek(filepos)
        linescontent=Sitesfile.read(self.SitesIndexMap[chrom][1]-filepos)
        print("need check encoding is utf-8")
        Siteslineslist=re.split(r"\n",linescontent.strip())
        for line in Siteslineslist:
            linelist = re.split(r'\t', line.strip())		
            c_chrom = linelist[5]
            pos = int(linelist[6])
            if pos>endpos or chrom!=c_chrom:
                break
            id=linelist[1].strip()
            TilingOrder = linelist[3].strip()
            importance=linelist[4]
            try:
                MAF = float(linelist[9].strip())
            except:
                MAF=0
            SitesList_A_Chrom.append((pos,id,TilingOrder,importance,MAF))
        Sitesfile.close()
        print("getSitesListByChrom",chrom," has total ",len(SitesList_A_Chrom)," filepos is :",self.SitesIndexMap[chrom]," but ",self.NumOfRecbychromOrder[i],"total recs in this file belong to this chrom") 
        return SitesList_A_Chrom


