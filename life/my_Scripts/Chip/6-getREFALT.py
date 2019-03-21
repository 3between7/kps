import re

f=open("../7-getGenotype/CattlePickedRecord_v1.GTSeq")
f1=open("../7-getGenotype/CattlePickedRecord_v1.REFALT","w")
#f2=open("../7-getGenotype/CattlePickedRecord_v1.REFALT.wrong","w")
f1.write("Chr	Pos	REF	ALT\n")
for eachline in f:
	ll=re.split("\t",eachline.strip())
	if ll[2] in ll[3]:
		if ll[2]==ll[3][0]:
			f1.write(ll[0]+"\t"+ll[1]+"\t"+ll[2]+"\t"+ll[3][2]+"\n")	
		elif ll[2]==ll[3][2]:
			f1.write(ll[0]+"\t"+ll[1]+"\t"+ll[2]+"\t"+ll[3][0]+"\n")
	else:
		f1.write(ll[0]+"\t"+ll[1]+"\t"+ll[3][0]+"\t"+ll[3][2]+"\n")

f.close()
f1.close()
#f2.close()
