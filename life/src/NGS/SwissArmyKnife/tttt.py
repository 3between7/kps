# -*- coding: UTF-8 -*-
'''
Created on 2018年6月1日

@author: Dr.liu
'''

if __name__ == '__main__':
    f=open("chr11list",'r');vf=open('zhongxinyihao.sorted.format.vcf','r');l=[]
    chr11f=open("chr11f.vcf",'w');chr10f=open("chr10f.vcf",'w')
    for line in f:
        l.append(line.strip())
    for line in vf:
        for e in l:
            if e in line:
                print(line.strip(),file=chr11f)
                break
        else:
            print(line.strip(),file=chr10f)
    chr11f.close();chr10f.close()