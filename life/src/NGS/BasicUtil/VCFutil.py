# -*- coding: UTF-8 -*-
'''
Created on 2013-6-30

@author: rui
'''


import random, sys
import re, pickle, copy

# from pip._vendor.requests.packages import chardet


class VCF_Data():
    def __init__(self, vcffileName):
        super().__init__()
        self.VcfMap_AllChrom = {}
        self.VcfIndexMap = {}
        self.chromOrder = []
        self.vcfFileName=vcffileName
        self.NumOfRecbychromOrder = []
        try:
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
        except:
            VCF_Data.indexVCF(VCFName=vcffileName, indexFileName=(vcffileName + ".myindex"))
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
        self.chromOrder = self.VcfIndexMap["chromOrder"]
        self.NumOfRecbychromOrder = self.VcfIndexMap["NumOfRecbychromOrder"]
    @staticmethod
    def indexVCF(VCFName, indexFileName):
        """
        {chrom:position_in_file_of_first_SNP_of_this_chrom,chrom:position,,,,,,}
        """
        vcffile = open(VCFName, 'r')#以可读模式打开VCF文件
        vcfChromIndex = {}
        chromOrder = []
        NumOfRecbychromOrder = []
        line = vcffile.readline()#读取vcffile文件的第一行
        vcfChromIndex["header"] = [line]#创建一个字典，value是一个列表！
        while re.search(r'^##', line) != None:#^匹配开始位置，^##匹配注释？？
            line = vcffile.readline()
            vcfChromIndex["header"].append(line)#把所有的注释全部存在key为header的value中
        
        if re.search(r'^#CHROM', line) != None:#re.search(pattern, string, flags=0) re.search函数会在字符串内查找模式匹配,只要找到第一个匹配然后返回，如果字符串没有匹配，则返回None。
            vcfChromIndex["title"] = re.split(r'\s+', line.strip()) #匹配任何空白字符:[<空格>\t\r\n\f\v]#按照能够匹配的子串将string分割后返回列表。可以使用re.split来分割字符串，如：re.split(r'\s+', text)；将字符串按空格分割成一个单词列表。格式：
                                                                    #re.split(pattern, string[, maxsplit])maxsplit用于指定最大分割次数，不指定将全部分割。
                                                                    #Python strip() 方法用于移除字符串头尾指定的字符（默认为空格或换行符）或字符序列。
                                                                    #把vcf的title的每一列分开未单独的元素存在key为title的value中
        else:                                                      
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
            exit(-1)  #????      
        currentChrom = "temptodele"
        lastPosition = vcffile.tell()#判断文件指针位置
        lastChromend_currentChromstartPostion = lastPosition
#         vcfChromIndex[currentChrom]=(lastPosition,0)
        print(line)
        line = vcffile.readline()
        print("first line:",line)
        i = 0
        while line:  #当while条件不成立，直接跳到else处输出    
            linelist = re.split(r"\s+", line) #把每一行具体数据也分隔开存储为一个列表
            if currentChrom != linelist[0]:#初始currentChrom = "temptodele"肯定与linelist[0]不等
                chromOrder.append(currentChrom)
                NumOfRecbychromOrder.append(i);i = 0  # collect the  number of snp recs  of the last chrom
                vcfChromIndex[currentChrom] = (lastChromend_currentChromstartPostion, lastPosition)
                lastChromend_currentChromstartPostion = lastPosition
                currentChrom = linelist[0]
                
#                 vcfChromIndex[currentChrom] = (lastPosition,0)
            
            lastPosition = vcffile.tell()

            line = vcffile.readline()
            i += 1
        else:
            chromOrder.append(currentChrom)
            NumOfRecbychromOrder.append(i - 1)  # collect the  number of snp recs  of the lastest chrom of all chroms
            vcfChromIndex[currentChrom] = (lastChromend_currentChromstartPostion, lastPosition)
#'8': (25757955057, 29473584687), '9': (29473584687, 33784528369)
# lastPosition and the next lastChromend_currentChromstartPostion are the same , modfiy it
        vcfChromIndex.pop("temptodele")
        i = chromOrder.index("temptodele")
        if i != 0:
            print("wrong indexVCF")
            exit(-1)
        b = chromOrder.pop(i)
        a = NumOfRecbychromOrder.pop(i)
#         print(i, a, b)
        vcfChromIndex["chromOrder"] = chromOrder
        
        vcfChromIndex["NumOfRecbychromOrder"] = NumOfRecbychromOrder
        pickle.dump(vcfChromIndex, open(indexFileName, 'wb'))
        vcffile.close()
#     def extractVcfRecByChroms(self,vcfFileName,chromlist,replacechromlist,outfile):
#         vcfFile = open(vcfFileName, 'r')
#         if len(chromlist)!=len(replacechromlist):
#             print("the length of the chromlist and replacechromlist should be the same")
#             return
#         for chrom in chromlist:
#             vcfFile.seek(self.VcfIndexMap[chrom][0])
    @staticmethod        
    def Vcf2geno_snp_ind(vcfFileName, sampleID_to_popmap, outputfileprefix, software, Morgenperbp, outputfilepart, VcfIndexMap=None, genotypesep=" ", withheader=False):
        vcffile = open(vcfFileName, "r")
        Morgenperbp=float(Morgenperbp)
        genofile = open(outputfileprefix + str(outputfilepart) + ".geno", "w")
        snpfile = open(outputfileprefix + str(outputfilepart) + ".snp", "w")
        indfile = open(outputfileprefix + str(outputfilepart) + ".ind", 'w')
        snpPositionlist = []
        pedmap = {}
        genolistOrderbySamplelist = []
        if withheader:
            line = vcffile.readline()
            while re.search(r"^##", line) != None:
                line = vcffile.readline()
            title = re.split(r"\s+", line.strip())
            total_individ = len(title) - 9
            print("Vcf2geno_snp_ind", title, len(title), total_individ)
            for sampleName in title[len(title) - total_individ:]:
                print(sampleName, "M", sampleID_to_popmap[sampleName], sep="\t", file=indfile)
        else:
            title = VcfIndexMap["title"]
            total_individ = len(VcfIndexMap["title"]) - 9
            print("Vcf2geno_snp_ind", VcfIndexMap["title"], len(VcfIndexMap["title"]), total_individ)
            for sampleName in VcfIndexMap["title"][len(VcfIndexMap["title"]) - total_individ:]:
                print(sampleName, "M", sampleID_to_popmap[sampleName], sep="\t", file=indfile)
        for line in vcffile:
            linelist = re.split(r"\s+", line)
            if linelist[3].strip().upper() == 'N' or len(linelist[3].strip()) > 1 or len(linelist[4].strip()) > 1:  # when ref is N ,or INDEL ,or multiple allels 
                continue
            print("\t" + linelist[0] + "_" + linelist[1] + "\t" + linelist[0] + "\t" + str(round(Morgenperbp * int(linelist[1]),6)) + "\t" + linelist[1] + "\t" + linelist[3] + "\t" + linelist[4], file=snpfile)
            if software.upper() == "GATK":
                GT_idx = (re.split(":", linelist[8])).index("GT")
                PL_idx = (re.split(":", linelist[8])).index("PL")
                genolistOrderbySamplelist = []
                for i in range(total_individ):
                    sample = linelist[i + 9]
                    if len(re.split(":", sample)) == 1 or re.split(":", sample)[GT_idx] == "./." or  len(re.split(r",", re.split(":", sample)[PL_idx])) != 3:  # ./.
                        genolistOrderbySamplelist += ['9']
                    else:
                        pl = re.split(":", sample)[PL_idx]
                        genotype = re.split(":", sample)[GT_idx]
                        a1 = int(re.search(r"(\d)/(\d)", genotype).group(1))
                        a2 = int(re.search(r"(\d)/(\d)", genotype).group(2))
                        if a1 == a2 and a1 == 1:
                            genolistOrderbySamplelist += ['0']
                        elif a1 == a2 and a1 == 0:
                            genolistOrderbySamplelist += ['1']
                        elif a1 != a2:
                            genolistOrderbySamplelist += ['2']
            print(*genolistOrderbySamplelist, sep=genotypesep, file=genofile)
        indfile.close()
        genofile.close()
        snpfile.close()
    @staticmethod
    def Vcf2Ped_WapperForpoolthreads(a):
        positionlist,pedmap=VCF_Data.Vcf2Ped(**a)
        print(a["chrom"],"finished")
        return positionlist,pedmap
    @staticmethod
    def Vcf2Ped(vcfFileName, software,VcfIndexMap=None, withheader=False, chrom=None,affectedlist=[],unaffectedlist=[],chromchangemap={},geno=1,maf=0,excludesits=[]):
        #when chromchangemap is not empty, then change if chrom in to the value corresponding  the k in chromchangemap
        vcffile = open(vcfFileName, "r")
#         mapfile = open(outputfileprefix + ".map", "w")
#         pedfile = open(outputfileprefix + ".ped", "w")
        positionlist = []
        pedmap = {}
        if withheader:
            line = vcffile.readline()
            while re.search(r"^##", line) != None:
                line = vcffile.readline()
            title = re.split(r"\s+", line.strip())
            total_individ = len(title) - 9
            print(title, len(title), total_individ)
            for outName in title[len(title) - total_individ:]:
                pedmap[outName] = []
        else:
            title = VcfIndexMap["title"]
            total_individ = len(VcfIndexMap["title"]) - 9
            print(VcfIndexMap["title"], len(VcfIndexMap["title"]), total_individ)
            for outName in VcfIndexMap["title"][len(VcfIndexMap["title"]) - total_individ:]:
                pedmap[outName] = []

        currentChromSome = None
        excludesitsForAchr=[]
        if chrom!=None:
            vcffile.seek(VcfIndexMap[chrom][0])
            excludesitsForAchr=excludesits
        for line in vcffile:
            linelist = re.split(r"\s+", line)
            if linelist[0].strip()!=chrom:
                break
            if linelist[3].strip().upper() == 'N' :#or len(linelist[3].strip()) > 1 or len(linelist[4].strip()) > 1:  # when ref is N ,or INDEL ,or multiple allels 
                print("exclude those sites which ref is not N ,INDEL ,or multiple alleles")
                continue
            AF=float(re.search(r"AF=([\d\.]+)[;,]", linelist[7]).group(1))
            AN = float(re.search(r"AN=([\d]+);", linelist[7]).group(1))
            if (total_individ*2-AN)/(total_individ*2)>geno or AF<maf or int(linelist[1].strip()) in excludesitsForAchr:
#                 print(line,total_individ*2-AN)
                if int(linelist[0])>50000000:
                    print(line);sys.stdout.flush()
                continue
            if linelist[0] in chromchangemap:
                positionlist.append((chromchangemap[linelist[0]],chromchangemap[linelist[0]]+"_"+linelist[1],0,linelist[1]))
            else:
                positionlist.append((re.sub(r"[\D\W]+","",linelist[0]),re.sub(r"\..*$","",linelist[0])+"_"+linelist[1],0,linelist[1]))
#             positionlist.append((linelist[0].replace("scaffold", ""), linelist[0] + "_" + linelist[1], 0, linelist[1]))
            if software.upper() == "GATK":

                GT_idx = (re.split(":", linelist[8])).index("GT")  # gatk GT:AD:DP:GQ:PL
                PL_idx = (re.split(":", linelist[8])).index("PL")
                for i in range(total_individ):
                    sample = linelist[i + 9]
                    if len(re.split(":", sample)) == 1 or re.split(":", sample)[GT_idx] == "./." or  len(re.split(r",", re.split(":", sample)[PL_idx])) != 3:  # ./.
                        pl = "0,0,0"
                    else:
                        pl = re.split(":", sample)[PL_idx]
                        genotype = re.split(":", sample)[GT_idx]                
                    if pl != "0,0,0":
                        a1 = int(re.search(r"(\d)/(\d)", genotype).group(1))
                        a2 = int(re.search(r"(\d)/(\d)", genotype).group(2))
                        alle1 = linelist[a1 + 3].strip()
                        alle2 = linelist[a2 + 3].strip()
                        pedmap[title[9 + i]] += [alle1, alle2]
                    else:
                        pedmap[title[9 + i]] += ['0', '0']
            
            elif software.upper() == "SAMTOOLS":
                pass
                for i in range(total_individ):
                    Sample = linelist[i + 9]
                    genotype = re.search(r"([^:]+):([^:]+):([^:]+)", Sample.strip()).group(1)
                    pl = re.search(r"([^:]+):([^:]+):([^:]+)", Sample.strip()).group(2)
                    if pl != "0,0,0":
                        a1 = int(re.search(r"(\d)/(\d)", genotype).group(1))
                        a2 = int(re.search(r"(\d)/(\d)", genotype).group(2))
                        alle1 = linelist[a1 + 3].strip()
                        alle2 = linelist[a2 + 3].strip()
                        pedmap[title[9 + i]] += [alle1, alle2]
                    else:
                        pedmap[title[9 + i]] += ['0', '0']
#         pickle.dump(pedmap,open(chrom+".pedmap", 'wb'))
#         pickle.dump(positionlist,open(chrom+".positionlist", 'wb'))
        return positionlist,pedmap
#         for elem in positionlist:
#             print(elem[0], elem[1], elem[2], elem[3], sep='\t', file=mapfile)
#         i = 1
#         for name in sorted(pedmap.keys()):
#             print(i, name, "0", "0", "1", "1", "\t".join(pedmap[name]), sep='\t', file=pedfile)
#             i += 1       
#         mapfile.close()
#         pedfile.close()   
#         vcffile.close()
    def getVcfListByChrom(self, chrom,startpos=1,endpos=9999999999999999999999999999999999999999999999999999, dilute=1, dilutetodensity="noofsnpperkb", posUniq=True, considerINDELandmultpleallele=False,MQfilter=28):
        """
            although dilute and dilutetodensity can exist at the same time,but it not make sense and may final produce a bug.
            return a list that contain all vcf record of a chrom
        """
        print("getVcfListByChrom",chrom,startpos,endpos,dilute,dilutetodensity)

        VcfList_A_Chrom = []
        if (chrom not in self.chromOrder) or startpos>=endpos:
            return []
            print(chrom + "didn't find in " + self.vcfFileName)
        i = self.chromOrder.index(chrom.strip())
        if dilute != 1:
            VcfRecRandomSelectIdxlist = random.sample([j for j in range(self.NumOfRecbychromOrder[i])], int(dilute * self.NumOfRecbychromOrder[i]))
            VcfRecRandomSelectIdxlist.sort()


        vcfFile = open(self.vcfFileName, 'r') #打开真正的vcf源文件
                    
        vcfFile.seek(self.VcfIndexMap[chrom][0])#？
        #find the first line
        filepos=vcfFile.tell()
        line=vcfFile.readline()
        
        while line:
#         find the startpos

            linelist = re.split(r'\s+', line.strip())
            samples = linelist[9:len(linelist)]
            c_chrom = linelist[0].strip()
            pos = int(linelist[1].strip())
            REF = linelist[3].strip()
            ALT = linelist[4].strip()
            if chrom.strip()!=c_chrom:
                return VcfList_A_Chrom            
            if pos>=startpos : 
                break

            filepos=vcfFile.tell()
            line=vcfFile.readline()
        else:
            if linelist is  None or pos<startpos:
                return VcfList_A_Chrom
        #from startpos to collected recs    
        vcfFile.seek(filepos)
        linescontent=vcfFile.read(self.VcfIndexMap[chrom][1]-filepos)
        print("need check encoding is utf-8")
#         encoding=chardet.detect(linescontent)
#         print(chardet.detect(linescontent))
#         if encoding["encoding"]!='utf-8':
#             exit(-1)
#         print("in it",filepos,self.VcfIndexMap[chrom][1]-filepos-1,encoding["encoding"])
        vcflineslist=re.split(r"\n",linescontent.strip())
#         line = vcfFile.readline().strip()
        recidx = 0
#         print(line)
        for line in vcflineslist:
#         while line and (re.split(r'\s+', line))[0] == chrom:
            if dilute != 1 and len(VcfRecRandomSelectIdxlist) == 0:
                break
            elif dilute != 1 and recidx != VcfRecRandomSelectIdxlist[0]:
                recidx += 1#line = vcfFile.readline();
                continue
            elif dilute != 1 and recidx == VcfRecRandomSelectIdxlist[0]:
                VcfRecRandomSelectIdxlist.pop(0)
                
            # the code block blow will be excute only when dilute==1,or    dilute!=1 and recidx == VcfRecRandomSelectIdxlist[0]:
            linelist = re.split(r'\s+', line.strip())
            samples = linelist[9:len(linelist)]
            c_chrom = linelist[0].strip()
            pos = int(linelist[1].strip())
            if pos>endpos or chrom!=c_chrom:
                break
            REF = linelist[3].strip()
            ALT = linelist[4].strip()
            recidx += 1#line = vcfFile.readline();
            if not considerINDELandmultpleallele and (len(REF) > 1 or len(ALT)> 1) :
                continue
            INFO = linelist[7]
            FORMAT = linelist[8]
            if linelist[6].strip()=="LowQual":
                continue
            if MQfilter!=None and re.search(r"MQ=([\d\.]+);",INFO)!=None:
                MQvalue=float(re.search(r"MQ=([\d\.]+);",INFO).group(1))
                if MQvalue<MQfilter:
                    continue
           # recidx += 1#line = vcfFile.readline();
            if posUniq and VcfList_A_Chrom and pos == VcfList_A_Chrom[-1][0]:
#                 print("VCFutil unique the vcf pos", line[:110].replace("\t"," "), VcfList_A_Chrom[-1][0:3],VcfList_A_Chrom[-1][3][:30])
                continue
            VcfList_A_Chrom.append((pos, REF, ALT, INFO, FORMAT, samples))
        vcfFile.close()
        if dilutetodensity != "noofsnpperkb":
            print("need to be fixed")
            dilute = dilutetodensity * (VcfList_A_Chrom[-1][0] - VcfList_A_Chrom[0][0]) / (1000 * len(VcfList_A_Chrom))
#             print(dilute)
            if dilute < 1:
                VcfRecRandomSelectIdxlist = random.sample([j for j in range(len(VcfList_A_Chrom))], int(dilute * len(VcfList_A_Chrom)))
                VcfRecRandomSelectIdxlist.sort()
            new_VcfList_A_Chrom = []
            for recidx in VcfRecRandomSelectIdxlist:
                new_VcfList_A_Chrom.append(VcfList_A_Chrom[recidx])
            VcfList_A_Chrom = copy.deepcopy(new_VcfList_A_Chrom)
        print("getVcfListByChrom",chrom, len(VcfList_A_Chrom), self.VcfIndexMap[chrom], self.NumOfRecbychromOrder[i], "total recs in this vcf belong to this chrom,dilute to", int(dilute * self.NumOfRecbychromOrder[i]),sep="\t") 
        if VcfList_A_Chrom!=[]:
            print(VcfList_A_Chrom[-1][0])
        return VcfList_A_Chrom
        

    def getVcfMap(self, vcfFileName):
        """
        this func is from bio\test\posAroundGene\func.py ,and did some improvement,that is add  INFO = collist[7],and add INFO into
        read the vcffile into a map which keys are chrom,values are a list of tuple
        {chrNo:[(pos,REF,ALT,INFO),(pos,REF,ALT,INFO),,,,,],chrNo:[],,,,,,},the order of the tuples in the list,is according pos,
        you we can search a record by  binary chop search
        no matter self.VcfMap_AllChrom has or has not value,the value will be clean
        """
        vcfMap = {}
        vcfFile = open(vcfFileName, 'r')
        
        line = vcfFile.readline()
        while re.search(r'^##', line) != None:
    #        print(line)
            line = vcfFile.readline()
        if re.search(r'^#', line) != None:
            lineslist = vcfFile.readlines()
        else:
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'\n" + line)
            exit(-1)
    #    print("pass")
        currentLine = 0
        totalRecs = len(lineslist)
        while currentLine != totalRecs:
    #        print(currentLine)
            collist = re.split(r'\s+', lineslist[currentLine])
            samples = collist[9:len(collist)]
            chrom = collist[0].strip()
            pos = int(collist[1].strip())
            REF = collist[3].strip()
            ALT = collist[4].strip()
            INFO = collist[7]
            FORMAT = collist[8]
            if chrom in vcfMap:
                vcfMap[chrom].append((pos, REF, ALT, INFO, FORMAT, samples))
            else:
                vcfMap[chrom] = [(pos, REF, ALT, INFO, FORMAT, samples)]
            currentLine += 1
        vcfFile.close()
        self.VcfMap_AllChrom = vcfMap
        








