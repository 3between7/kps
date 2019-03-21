# -*- coding: UTF-8 -*-
import re
def replace(file1,file2,file3):
    neutraldict={}
    contentdict={}
    with open('600ATCG.txt',"r") as f0:
            file1list=list(f0)
            listrecordATCG=[]
            for snpid in file1list:
                listrecordATCG.append(snpid.replace('\n', ''))
            print(len(listrecordATCG))    
    with open(file1) as f1:
        filelist=list(f1)
        n=len(filelist)
        NumOfATCG=0
        NumOfATCGCanBeReplaced=0
        for i in range(n-3):
            listline=re.split(r"\n",filelist[i].strip())
            if listline[0] in listrecordATCG:
                NumOfATCG +=1
                listlineb1=re.split(r"\s+",filelist[i-1].strip())
                listlineb2=re.split(r"\s+",filelist[i-2].strip())
                listlineb3=re.split(r"\s+",filelist[i-3].strip())
                listlinea1=re.split(r"\s+",filelist[i+1].strip())
                listlinea2=re.split(r"\s+",filelist[i+2].strip())
                listlinea3=re.split(r"\s+",filelist[i+3].strip())
                if (int(listline[1]) < 3700) and (-int(listline[2]) > 3700):
                    if listlineb1[0] not in listrecordATCG:
                        neutraldict[listlineb1[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif (listlineb2[0] not in listrecordATCG) and (int(listline[2])+int(listlineb1[2]) <=3700):
                        neutraldict[listlineb2[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif (listlineb3[0] not in listrecordATCG) and (int(listline[2])+int(listlineb1[2]) + int(listlineb2[2])<=3700):
                        neutraldict[listlineb3[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    else:
                        continue
                elif (int(listline[1]) > 3700) and (-int(listline[2]) < 3700):
                    if listlinea1[0] not in listrecordATCG:
                        neutraldict[listlinea1[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif (listlinea2[0] not in listrecordATCG) and -(int(listline[2])-int(listlinea1[2])) <=3700:
                        neutraldict[listlinea2[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif (listlinea3[0] not in listrecordATCG) and -(int(listline[2])-int(listlinea1[2])-int(listlinea2[2]))<=3700:
                        neutraldict[listlinea3[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    else:
                        continue
                elif int(listline[1]) <= -int(listline[2]) <= 3700:
                    if listlineb1[0] not in listrecordATCG:
                        neutraldict[listlineb1[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif listlinea1[0] not in listrecordATCG:
                        neutraldict[listlinea1[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif (listlineb2[0] not in listrecordATCG) and (int(listline[1])+int(listlineb1[1]) <=3700):
                        neutraldict[listlineb2[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif (listlinea2[0] not in listrecordATCG) and -(int(listline[2])-int(listlinea1[2])) <=3700:
                        neutraldict[listlinea2[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif (listlineb3[0] not in listrecordATCG) and (int(listline[1])+int(listlineb1[1]) + int(listlineb2[1])<=3700):
                        neutraldict[listlineb3[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif (listlinea3[0] not in listrecordATCG) and -(int(listline[2])-int(listlinea1[2])-int(listlinea2[2]))<=3700:
                        neutraldict[listlinea3[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    else:
                        continue
                elif -int(listline[2]) < int(listline[1]) < 3700:
                    if listlinea1[0] not in listrecordATCG:
                        neutraldict[listlinea1[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif listlineb1[0] not in listrecordATCG:
                        neutraldict[listlineb1[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif (listlinea2[0] not in listrecordATCG) and -(int(listline[2])-int(listlinea1[2])) <=3700:
                        neutraldict[listlinea2[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif (listlineb2[0] not in listrecordATCG) and (int(listline[1])+int(listlineb1[1]) <=3700):
                        neutraldict[listlineb2[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    elif (listlinea3[0] not in listrecordATCG) and -(int(listline[2])-int(listlinea1[2])-int(listlinea2[2]))<=3700:
                        neutraldict[listlinea3[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1                    
                    elif (listlineb3[0] not in listrecordATCG) and (int(listline[1])+int(listlineb1[1]) + int(listlineb2[1])<=3700):
                        neutraldict[listlineb3[0]]=listline[0]
                        NumOfATCGCanBeReplaced+=1
                    else:
                        continue
            else:
                continue
                                                      
    print("There are "+str(NumOfATCG)+" ATCG records in "+file1+" , "+str(NumOfATCGCanBeReplaced)+" can be repalced.")

    with open(file2) as f2:
        NumOfRecordbePicked=0
        for eachline in f2:
            eachlinelist=re.split(r"\s+",eachline.strip())
            if eachlinelist[0] in neutraldict.keys():
                NumOfRecordbePicked+=1
                value=neutraldict[eachlinelist[0]]
                contentdict[value]=eachline
            else:
                continue
    print("There are "+str(NumOfRecordbePicked)+"  records in "+file2+" that picked sucessfully.")
    
    with open(file3) as f3:
        NumOfRecordbeReplaced=0
        f4=open("resultATCG","w")
        for eachline in f3:
            snpid=re.split(r"\s+",eachline.strip())[0]
            if snpid in contentdict.keys():
                NumOfRecordbeReplaced+=1
                recordcontent=contentdict[snpid]
                f4.write(recordcontent)
            else:
                f4.write(eachline)
    print("There are "+str(NumOfRecordbeReplaced)+" records in "+file3+" be replaced in fact.")

replace("file1new.txt","file2.txt","Axiom_KPSmilet_redo_scored.txt.sorted3700withtitle.miodifylastcol1_7netural")
                
        
                 
        
                