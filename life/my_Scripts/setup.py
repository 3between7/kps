#!python
#cython: language_level=3, boundscheck=False
'''
Created on 2018��9��10��

@author: Dr.liu
'''
# �������£�������ĿĿ¼�´������л���shell��������ֻ�ܱ���һ���ļ�������֮��ᷢ�ֳ��������ļ���yourmod.c��yourmod.html��yourmod-win_amd64.pyd����ʱ��c��html��ԭpy�ļ�ɾ������pyd�ļ���������Ϊyourmod�Ϳ��ԣ�
#cythonize -3 -a -i yourmod.pyx

import pyximport

from distutils.core import setup
from Cython.Build import cythonize
pyximport.install(pyimport=True,language_level =3)
import shutil
import os
import re

'''
���ļ���ִ����Ҫ����Terminal������   python setup.py build_ext --inplace ������
ʹ��Cpython ����python�ļ����ؼ����������pyd�ļ����൱��dll��
'''
# ��Զ��ļ�������ã����ļ���ֻдһ������
key_funs = ["/lustre/home/liurui/software/life/src/pipelinecontrol/Util.py", "/lustre/home/liurui/software/life/src/test.py", "/lustre/home/liurui/software/life/src/NGS/BasicUtil/VCFutil.py","/lustre/home/liurui/software/life/src/NGS/BasicUtil/Util.py",
            "/lustre/home/liurui/software/life/src/NGS/BasicUtil/DBManager.py", "/lustre/home/liurui/software/life/src/NGS/BasicUtil/geneUtil.py", "/lustre/home/liurui/software/life/src/NGS/BasicUtil/Caculators.py",
            "/lustre/home/liurui/software/life/src/web/views.py", "/lustre/home/liurui/software/life/src/web/Service.py","/lustre/home/liurui/software/life/src/web/forms.py", "/lustre/home/liurui/software/life/src/web/dba.py", "/lustre/home/liurui/software/life/src/web/DBA.py","/lustre/home/liurui/software/life/src/web/Entity.py"]

#setup(
 #   name="life project", 
#    ext_modules = cythonize(key_funs,annotate=True,exclude_failures=True,compiler_directives={"language_level":3})
#   )

'''
1����������so�ļ����������ĳ���ԭpy�ļ�һ��
2��ɾ�������õ���c�ļ���ԭpy�ļ�
'''
rootDir=os.getcwd()
print("������������", rootDir, "������������")


def list_all_files(rootdir):
#     import os
    _files = []
    list = os.listdir(rootdir) #�г��ļ��������е�Ŀ¼���ļ�
    for i in range(0,len(list)):
           path = os.path.join(rootdir,list[i])
           if os.path.isdir(path):
              _files.extend(list_all_files(path))
           if os.path.isfile(path):
              _files.append(path)
    return _files

files = list_all_files(rootDir)
for fi in files:
    if fi.__contains__(".so"):
        re_name=fi.split(".")[0]+".so"
        print(re_name)
        os.rename(fi,re_name)
        try:
            fi_path=re.sub(r"/\w*\.\w*\.so","",fi.split("src/life/")[0]+fi.split("src/life/")[1])
            shutil.copy(fi,fi_path)
        except:
            pass
    elif fi.__contains__(".c") or fi in key_funs:
        os.remove(fi)
#print(files)
"gcc -pthread -shared /lustre/home/liurui/software/life/src/build/temp.linux-x86_64-3.6/NGS/BasicUtil/VCFutil.o -o /lustre/home/liurui/software/life/src/NGS/BasicUtil/VCFutil.so"
