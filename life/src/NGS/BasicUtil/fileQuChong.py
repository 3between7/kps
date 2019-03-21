# -*- coding: UTF-8 -*-
'''file1 is sheet5, file2 is reslut1'''
def uniq(file1,file2):
    f1=open(file1)
    f2=open(file2)
    f3=open(file1.split('.')[0]+"quchong.txt",'w')
    file2list=list(f2)
    for line in f1:
        if line in file2list:
            continue
        else:
            f3.write(line)
    f1.close()
    f2.close()       
    f3.close()

uniq("file1.txt","file2.txt")
            
            
    
    
    
    
    
