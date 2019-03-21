import os
def clean(indir,outdir,lengthofreads):
	#allfiles = os.walk(indir);print(allfiles[2])
	for item in os.walk(indir):
		for f in item[2]:
			if "SRR" or "ERR" in item[2]:
				SRR=f.split("/")[-1].split("_")[0];print(SRR)
				if ".gz" in f:
					os.system("fastp -i " + indir + "/"+SRR + "_1.fastq.gz -I " + indir + "/"+SRR + "_2.fastq.gz -o " + outdir + "/clean_" + SRR + "_1.fastq.gz -O " + outdir + "/clean_" + SRR + "_2.fastq.gz -u 50 -n 1 -q 15 -l " + lengthofreads)
				else:
					os.system("fastp -i "  + indir + "/"+SRR + "_1.fastq -I " + indir + "/"+SRR + "_2.fastq -o " + outdir + "/clean_" + SRR + "_1.fastq -O " + outdir + "/clean_" + SRR + "_2.fastq -u 50 -n 1 -q 15 -l " + lengthofreads)

if __name__ == '__main__':
	clean("/lustre/project/R2D/Sheep50k-CAAS-zhangli/0-rawdata/fastq.gz","/home/zhanghuanhuan","126")	
