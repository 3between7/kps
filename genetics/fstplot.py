import csv,sys
import matplotlib.pyplot as plt
from geneview.gwas import manhattanplot
import numpy as np

def fstplot(infile,outfile):
	with open(infile) as f:
		f_csv = csv.reader(f)
		headers=next(f_csv)
		data=[[row[0],int(row[1]),float(row[2])] for row in f_csv]

	fig = plt.figure(figsize=(15,2)) #设置图片的宽和高
	ax1 = fig.add_subplot(1,1,1)
	ax1.grid(False) #网格线开关
	ax = manhattanplot(data, xlabel="Chromosome", ylabel="fst",
					xtick_label_set=('2','4','6','8','10','12','14','16','18','20','22','24','26','28','30'),
					s=1,mlog10=False)
	plt.xticks(rotation=90) # x轴下面的文字旋转90度
	plt.savefig(outfile)
	plt.show()

if __name__ == "__main__":
	fstplot(sys.argv[1],sys.argv[2])
