'''
updated on 2019.03.11

@author: zhanghuanhuan
'''

import re,sys
import pandas as pd

orif =pd.read_csv(sys.argv[1])
addprif=pd.read_csv(sys.argv[2],sep="\t")
result=pd.merge(orif,addprif,how="left",left_on="snpid",right_on="snpid")
result.to_csv(sys.argv[3],sep="\t",column=['priority'],header=1,index=0)
os.system("paste -t $'\t' " + orif + " " + sys.argv[3])
