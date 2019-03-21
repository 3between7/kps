# -*- coding: UTF-8 -*-
import re
'''file1 is zongwenjian,file2 is paixuwenjian '''
def replace(file1,file2):
    Totalfiledict={}
    snpidused=[]
    
    with open(file1,"r") as f1:
        for eachline in f1:
            linelist=re.split(r"\s+",eachline.strip())
            key=linelist[0]
            value=eachline.strip()
            Totalfiledict[key]=value

    with open(file2) as f2:
        '''file2 format : snpid    qiandis    houdis    seq'''
        filelist=list(f2)
        n=len(filelist)
        NumOfATCG=0
        NumOfATCGCanBeReplaced=0
        for i in range(n-2):
            listline=re.split(r"\s+",filelist[i].strip())
            line0content=Totalfiledict[listline[0]]
            line0priority=line0content.strip()[-1]
            if ("A/T" or "C/G") in listline[3]:
                NumOfATCG +=1
                listlineb1=re.split(r"\s+",filelist[i-1].strip())
                listlineb2=re.split(r"\s+",filelist[i-2].strip())
                listlinea1=re.split(r"\s+",filelist[i+1].strip())
                listlinea2=re.split(r"\s+",filelist[i+2].strip())
                
                if (int(listline[1]) < 3700) and (-int(listline[2]) > 3700):
                    
                    if ("A/T" or "C/G") not in listlineb1[3]:
                        b1content=Totalfiledict[listlineb1[0]]
                        b1priority=b1content.strip()[-1]
                        Totalfiledict[listlineb1[0]]=b1content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+b1priority
                        NumOfATCGCanBeReplaced+=1
                        
                    elif ("A/T" or "C/G") not in listlineb2[3] and (int(listline[1])+int(listlineb1[1]) <=3700):
                        b2content=Totalfiledict[listlineb2[0]]
                        b2priority=b2content.strip()[-1]                       
                        Totalfiledict[listlineb2[0]]=b2content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+b2priority
                        NumOfATCGCanBeReplaced+=1        
                    else:
                        continue
                    
                elif (int(listline[1]) > 3700) and (-int(listline[2]) < 3700):
                    if ("A/T" or "C/G") not in listlinea1[3]:                       
                        a1content=Totalfiledict[listlinea1[0]]
                        a1priority=a1content.strip()[-1]                        
                        Totalfiledict[listlinea1[0]]=a1content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+a1priority
                        NumOfATCGCanBeReplaced+=1
                    elif ("A/T" or "C/G") not in listlinea2[3] and -(int(listline[2])-int(listlinea1[2])) <=3700:
                        a2content=Totalfiledict[listlinea2[0]]
                        a2priority=a2content.strip()[-1]                        
                        Totalfiledict[listlinea2[0]]=a2content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+a2priority
                    else:
                        continue
                    
                elif int(listline[1]) <= -int(listline[2]) <= 3700:
                    if ("A/T" or "C/G") not in listlineb1[3]:
                        b1content=Totalfiledict[listlineb1[0]]
                        b1priority=b1content.strip()[-1]
                        Totalfiledict[listlineb1[0]]=b1content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+b1priority
                        NumOfATCGCanBeReplaced+=1
                    elif ("A/T" or "C/G") not in listlinea1[3]:
                        a1content=Totalfiledict[listlinea1[0]]
                        a1priority=a1content.strip()[-1]                        
                        Totalfiledict[listlinea1[0]]=a1content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+a1priority
                        NumOfATCGCanBeReplaced+=1
                    elif ("A/T" or "C/G") not in listlineb2[3] and (int(listline[1])+int(listlineb1[2]) <=3700):
                        b2content=Totalfiledict[listlineb2[0]]
                        b2priority=b2content.strip()[-1]                        
                        Totalfiledict[listlineb2[0]]=b2content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+b2priority
                    elif ("A/T" or "C/G") not in listlinea2[3] and -(int(listline[2])-int(listlinea1[2])) <=3700:
                        a2content=Totalfiledict[listlinea2[0]]
                        a2priority=a2content.strip()[-1]                        
                        Totalfiledict[listlinea2[0]]=a2content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+a2priority
                    else:
                        continue
                    
                elif -int(listline[2]) < int(listline[1]) < 3700:
                    if ("A/T" or "C/G") not in listlinea1[3]:
                        a1content=Totalfiledict[listlinea1[0]]
                        a1priority=a1content.strip()[-1]                        
                        Totalfiledict[listlinea1[0]]=a1content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+a1priority
                        NumOfATCGCanBeReplaced+=1
                    elif ("A/T" or "C/G") not in listlineb1[3]:
                        b1content=Totalfiledict[listlineb1[0]]
                        b1priority=b1content.strip()[-1]
                        Totalfiledict[listlineb1[0]]=b1content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+b1priority
                    elif ("A/T" or "C/G") not in listlinea2[3] and -(int(listline[2])-int(listlinea1[2])) <=3700:
                        a2content=Totalfiledict[listlinea2[0]]
                        a2priority=a2content.strip()[-1]                        
                        Totalfiledict[listlinea2[0]]=a2content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+a2priority
                    elif ("A/T" or "C/G") not in listlineb2[3] and (int(listline[1])+int(listlineb1[1]) <=3700):
                        b2content=Totalfiledict[listlineb2[0]]
                        b2priority=b2content.strip()[-1]                        
                        Totalfiledict[listlineb2[0]]=b2content[:-1]+line0priority
                        Totalfiledict[listline[0]]=line0content[:-1]+b2priority
                    else:
                        continue
                    
            else:
                continue
    f3=open("resultATCG.log","w")
    f3.write("There are "+str(NumOfATCG)+" Neutral records in "+file2+" , "+str(NumOfATCGCanBeReplaced)+" can be repalced.")
    for eachline in Totalfiledict.values():
        print(eachline.strip())
    f3.close()

replace("result1.txt","newtihuanwenjian.txt")
                
        
                 
        
                