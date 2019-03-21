# -*- coding: UTF-8 -*-
import re

''' the format of file1.txt is :snpid    best_recommendation     dishou-qian    disqian-hou'''

def tihuan(file0,file1,file2):
    netrualdict={}  # {"x1":x2}  x1: the SNP ID of 'netrual' that can be replace ;x2:the SNP ID of  "recomended" used be replaced
    tihuandict={}   # {"x1":x3}  x1: the SNP ID of 'netrual' that can be replace ;x3:the record used to be replaced 
    with open(file0,'r') as f:
        filelist=list(f)
        n=len(filelist)
        Numofnetrual=0
        for i in range(n-2): # notice the range is n-2 ,not n!
            linelist=re.split(r'\s+', filelist[i].strip())
            if linelist[1]=="neutral":
                Numofnetrual+=1
                linelistb1=re.split(r'\s+', filelist[i-1].strip())
                linelista1=re.split(r'\s+', filelist[i+1].strip())
                linelistb2=re.split(r'\s+', filelist[i-2].strip())
                linelista2=re.split(r'\s+', filelist[i+2].strip())
                if (int(linelist[2]) <= 3700) and linelistb1[1] != "neutral": 
                    netrualdict[linelist[0]]=linelistb1[0]
                elif (-int(linelist[3]) <= 3700) and linelista1[1] != "neutral": 
                    netrualdict[linelist[0]]=linelista1[0]                  
                elif (int(linelist[2])+int(linelistb1[2])) <=3700 and linelistb2[1] != "neutral":
                    netrualdict[linelist[0]]=linelistb2[0]                               
                elif ((-int(linelist[3]))+(-int(linelista1[3])))<=3700 and linelista2[1] != "neutral":
                     netrualdict[linelist[0]]=linelista2[0] 
                else:
                    continue   
            else:
                continue
        print(file0+" has "+str(Numofnetrual)+' neutral records totally,has '+str(len(netrualdict))+" records can be used to replace")
    with open(file1,'r') as f1:
        n1=0
        for eachline in f1:
            linelist1=re.split(r"\t",eachline.strip())
            if linelist1[0] in netrualdict.keys():
                tihuandict[linelist1[0]]=eachline
                n1+=1
        print(file1+" picked "+str(n1)+" records")
    
    with open(file2,"r") as f2:
        n2=0
        f3=open("result","w")
        for eachline in f2:
            linelist2=re.split(r"\t",eachline.strip())
            key=linelist2[0]
            if key in tihuandict.keys(): 
                f3.write(tihuandict[key])
                n2+=1
            else:
                f3.write(eachline)
        f3.close()
        print(file2+" has "+str(n2)+" be replaced.")
                
                
                
tihuan("file1.txt","Axiom_KPSmilet_redo_scored.txt.sorted3700withtitle.miodifylastcol","Axiom_KPSmilet_redo_scored.txt.sorted3700withtitle.miodifylastcol")
            
    
                        
                
        