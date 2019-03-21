# -*- coding: UTF-8 -*-
import re
'''file1 is whoisbepicked.txt, file2 is sheet5'''
def uniq(file1,file2):
    f1=open(file1)
    f2=open(file2)
    f3=open(file2.split('.')[0]+"filter.txt",'w')
    file1list=list(f1)
    listrecordpicked=[]
    for snpid in file1list:
        listrecordpicked.append(snpid.replace('\n', ''))    
    for line in f2:
        if re.split(r"\s+",line.strip())[0] in listrecordpicked:
            continue
        else:
            f3.write(line)
    f1.close()
    f2.close()       
    f3.close()

uniq("whoisbepicked.txt","sheet5.txt")
            
            
    
    
    
    
    
