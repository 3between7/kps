filename=input("please set input file name:")
f=open(filename)
outputfile=input("please set output filename:")
f1=open(outputfile,"w")

dictSeq={}
for eachline in f:
	if ">" in eachline:
		key=eachline.strip()[1:]
	else:
		value=eachline.strip()
		dictSeq[key]=value
for item in dictSeq.keys():
	print("please set:",item)
	f1.write("SEQUENCE_ID="+item+"\n")
	print("SEQUENCE_ID="+item)
	f1.write("SEQUENCE_TEMPLATE="+dictSeq[item]+"\n")
	print("SEQUENCE_TEMPLATE="+dictSeq[item][:10]+"..."+dictSeq[item][-10:])
	f1.write("PRIMER_TASK=generic"+"\n")
	print("PRIMER_TASK=generic")
	PRIMER_PICK_LEFT_PRIMER=input("PRIMER_PICK_LEFT_PRIMER :")
	f1.write("PRIMER_PICK_LEFT_PRIMER"+"="+PRIMER_PICK_LEFT_PRIMER+"\n")
	PRIMER_PICK_INTERNAL_OLIGO=input("PRIMER_PICK_INTERNAL_OLIGO :")
	f1.write("PRIMER_PICK_INTERNAL_OLIGO="+PRIMER_PICK_INTERNAL_OLIGO+"\n")
	PRIMER_PICK_RIGHT_PRIMER=input("PRIMER_PICK_RIGHT_PRIMER :")
	f1.write("PRIMER_PICK_RIGHT_PRIMER="+PRIMER_PICK_RIGHT_PRIMER+"\n")
	PRIMER_OPT_SIZE=input("PRIMER_OPT_SIZE : ")
	f1.write("PRIMER_OPT_SIZE="+PRIMER_OPT_SIZE+"\n")
	PRIMER_MIN_SIZE=input("PRIMER_MIN_SIZE :")
	f1.write("PRIMER_MIN_SIZE="+PRIMER_MIN_SIZE+"\n")
	PRIMER_MAX_SIZE=input("PRIMER_MAX_SIZE :")
	f1.write("PRIMER_MAX_SIZE="+PRIMER_MAX_SIZE+"\n")
	PRIMER_PRODUCT_SIZE_RANGE=input("PRIMER_PRODUCT_SIZE_RANGE :")
	f1.write("PRIMER_PRODUCT_SIZE_RANGE="+PRIMER_PRODUCT_SIZE_RANGE+"\n")
	PRIMER_PRODUCT_SIZE_RANGE=input("PRIMER_PRODUCT_SIZE_RANGE(for example:75-150) :")
	f1.write("PRIMER_EXPLAIN_FLAG="+PRIMER_EXPLAIN_FLAG+"\n")
	print("=======================================================")
print("done!")
f.close()
f1.close()
