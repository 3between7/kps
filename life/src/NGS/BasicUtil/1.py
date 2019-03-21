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

MyVcfData=VCF_Data("01_Nipponbare_B001.snp.vcf")
print(MyVcfData.chromOrder)