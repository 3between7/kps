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
        vcffile = open(VCFName, 'r')#�Կɶ�ģʽ��VCF�ļ�
        vcfChromIndex = {}
        chromOrder = []
        NumOfRecbychromOrder = []
        line = vcffile.readline()#��ȡvcffile�ļ��ĵ�һ��
        vcfChromIndex["header"] = [line]#����һ���ֵ䣬value��һ���б�
        while re.search(r'^##', line) != None:#^ƥ�俪ʼλ�ã�^##ƥ��ע�ͣ���
            line = vcffile.readline()
            vcfChromIndex["header"].append(line)#�����е�ע��ȫ������keyΪheader��value��
        
        if re.search(r'^#CHROM', line) != None:#re.search(pattern, string, flags=0) re.search���������ַ����ڲ���ģʽƥ��,ֻҪ�ҵ���һ��ƥ��Ȼ�󷵻أ�����ַ���û��ƥ�䣬�򷵻�None��
            vcfChromIndex["title"] = re.split(r'\s+', line.strip()) #ƥ���κοհ��ַ�:[<�ո�>\t\r\n\f\v]#�����ܹ�ƥ����Ӵ���string�ָ�󷵻��б�����ʹ��re.split���ָ��ַ������磺re.split(r'\s+', text)�����ַ������ո�ָ��һ�������б���ʽ��
                                                                    #re.split(pattern, string[, maxsplit])maxsplit����ָ�����ָ��������ָ����ȫ���ָ
                                                                    #Python strip() ���������Ƴ��ַ���ͷβָ�����ַ���Ĭ��Ϊ�ո���з������ַ����С�
                                                                    #��vcf��title��ÿһ�зֿ�δ������Ԫ�ش���keyΪtitle��value��
        else:                                                      
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
            exit(-1)  #????      
        currentChrom = "temptodele"
        lastPosition = vcffile.tell()#�ж��ļ�ָ��λ��
        lastChromend_currentChromstartPostion = lastPosition
#         vcfChromIndex[currentChrom]=(lastPosition,0)
        print(line)
        line = vcffile.readline()
        print("first line:",line)
        i = 0
        while line:  #��while������������ֱ������else�����    
            linelist = re.split(r"\s+", line) #��ÿһ�о�������Ҳ�ָ����洢Ϊһ���б�
            if currentChrom != linelist[0]:#��ʼcurrentChrom = "temptodele"�϶���linelist[0]����
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