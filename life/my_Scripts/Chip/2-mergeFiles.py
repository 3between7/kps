import re,pickle
import pandas as pd

MM=open("MMCattle_SNP.info")
Bf=open("BfCattle_SNP.info")
CY=open("CYCattle_SNP.info")

MMlist=[];Vm=[]
Bflist=[];Vb=[]
CYlist=[];Vc=[]

try:
	MMlist=pickle.load(open("MMlist.pickle","rb"))
	Vm=pickle.load(open("Vm.pickle","rb"))
	Bflist=pickle.load(open("Bflist.pickle","rb"))
	Vb=pickle.load(open("Vb.pickle","rb"))
	CYlist=pickle.load(open("CYlist.pickle","rb"))
	Vc=pickle.load(open("Vc.pickle","rb"))
except:
	for lm in MM:
		llm=re.split(r"\s+",lm.strip())
		MMlist.append(llm[0])
		Vm.append(llm[1])	
	pickle.dump(MMlist,open("MMlist.pickle","wb"))
	pickle.dump(Vm,open("Vm.pickle","wb"))
	print(len(MMlist))
	print(len(Vm))

	for lb in Bf:
		llb=re.split(r"\t",lb.strip())
		Bflist.append(llb[0])
		Vb.append(llb[1])
	pickle.dump(Bflist,open("Bflist.pickle","wb"))
	pickle.dump(Vb,open("Vb.pickle","wb"))
	print(len(Bflist))
	print(len(Vb))

	for lc in CY:
		llc=re.split(r"\s+",lc.strip())
		CYlist.append(llc[0])
		Vc.append(llc[1])
	pickle.dump(CYlist,open("CYlist.pickle","wb"))
	pickle.dump(Vc,open("Vc.pickle","wb"))
	print(len(CYlist))
	print(len(Vc))


dfm=pd.DataFrame({'key':MMlist,'value1':Vm})
dfb=pd.DataFrame({'key':Bflist,'value2':Vb})
dfc=pd.DataFrame({'key':CYlist,'value3':Vc})

data=pd.merge(dfm,dfb,how='outer') 
data_merge3=pd.merge(data,dfc,how='outer')
data_merge3.to_csv("result_3population.txt",index=False,header=False)	
