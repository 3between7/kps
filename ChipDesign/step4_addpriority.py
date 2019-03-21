'''
updated on 2019.03.07

@author: Dr.liu,zhanghuanhuan
'''

import os,re
from ChipDesign.tools import window,Caculators

                     
def addpriority(scoredfile,mustincludefile,winsize,thresholdfor_I,seqid=1,type="A"):
    
        print("python chip_design.py scoredfile winsize I/A thresholdfor_I")
        print("add prority accord tiling_order and win; i.e repriority")
        print("thresholdfor_I is for I only, but should give any value for A ")
        
        if type.strip().upper()=="A":
            
            #if there are sites must include in chip but they are not VIP sites!
            musctincludes=set()
            mustincludef=open(mustincludefile,'r');
            for line in mustincludef:
                mustinl=re.split(r"\t",line.strip())
                musctincludes.add(mustinl[1])
            mustincludef.close()
                   
            dupseqmap={}
            os.system("sed -n '1,1p' " +scoredfile+" > titlelinetemp")            
            print(scoredfile)
            f=open(scoredfile,'r');title=re.split(r"\t",f.readline().strip())
            tidx=title.index("tiling_order");curchr=title[4];print("curchr:",curchr,"tilingorder idx:",tidx)
            win = window.Window()
            ofo=open(scoredfile+str(winsize),'w')
            print(*title,"priority",sep="\t",file=ofo)
            addprortycaculator = Caculators.Caculator_addpriority(of=ofo,seqidx=seqid,tilingorderidx=tidx,best_recommendation=title.index("best_recommendation"),rmdupmap=dupseqmap,mustin=musctincludes)
            addprortycaculator.addVIPidx=title.index("importance")
            recs=[];count=0;tcount=0
            for line in f:
                tcount+=1
                recl=re.split(r"\t",line.strip())
                
                if curchr==recl[4]:
                    currentchrLen=int(recl[5])
                    recs.append([int(recl[5])]+recl)
                elif curchr!="cust_chr":                
                    count+=len(recs)
                    if count!=tcount-1:
                        print(line,count,tcount);exit(-1)
                    win.slidWindowOverlap(recs, currentchrLen, int(winsize), int(winsize), addprortycaculator)
                    recs=[[int(recl[5])]+recl];curchr=recl[4]
                    currentchrLen=int(recl[5])
                else:#first line
                    recs=[[int(recl[5])]+recl];curchr=recl[4];currentchrLen=int(recl[5])
            else:
                win.slidWindowOverlap(recs, currentchrLen, int(winsize), int(winsize), addprortycaculator)
            ofo.close()
            f.close();addprortycaculator.temp.close()
            
        elif type.strip().upper()=="I":
 
            win = window.Window()
            ofo=open(scoredfile+winsize,'w')
            f=open(scoredfile,'r')
            
            recs=[];count=0;tcount=0
            curchr=None
            for line in f:
                tcount+=1
                recl=re.split(r",",line.strip())
    #             print(recl)
                if recl[3]!="Chromosome" and curchr is None:
                    continue
                elif curchr is None:
                    curchr=recl[3];title=recl#"Chromosome"
                    tidx=title.index("Source_Version");print("curchr:",curchr,"tilingorder idx:",tidx);best_recommendation=title.index("Final_Score")
                    addprortycaculator = Caculators.Caculator_selectsitesForillumina(of=ofo,tilingorderidx=tidx,best_recommendation=best_recommendation,T=float(sys.argv[4]),VIPidx=title.index("Source"))
                    addprortycaculator.spf=open(scoredfile+winsize+".originscored",'w')
                    continue
                if curchr==recl[3]:
                    currentchrLen=int(recl[4])
                    recs.append([int(recl[4])]+recl)
                elif curchr!="Chromosome":
                    if curchr=="0":
                        for e in recs:
                            if e[title.index("Source")+1]=="VIP" or( int(e[tidx+1])<=3 and float(e[best_recommendation+1])>0.6): 
                                print(*e[1:],sep="\t",file=addprortycaculator.spf)
                                print(e[1],"SNP",e[2],e[4],e[5],e[3],e[6],e[7],e[8],addprortycaculator.m[e[8].upper()],"Soybean","FALSE",sep=",",file=ofo)
                    else:
                        win.slidWindowOverlap(recs, currentchrLen, int(winsize), int(winsize), addprortycaculator)
                    recs=[[int(recl[4])]+recl];curchr=recl[3]
                    currentchrLen=int(recl[4])
                else:#first line
                    recs=[[int(recl[4])]+recl];curchr=recl[3];currentchrLen=int(recl[4])
            else:
                win.slidWindowOverlap(recs, currentchrLen, int(winsize), int(winsize), addprortycaculator)

        exit() 