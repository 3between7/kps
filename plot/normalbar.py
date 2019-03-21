import pandas as pd
from matplotlib import pyplot as plt

def normalbar(input,output,title_x,title_y,xlabel,ylabel,set_title,set_xlim_l,set_xlim_r):
    
    data=pd.read_csv(input,sep="\s+")
    fig,ax=plt.subplots()
    
    ax.bar(data[title_x],data[title_y])
    ax.set_xlabel(xlabel)  #设置x轴标签
    ax.set_ylabel(ylabel)  #设置y轴标签
    ax.set_title(title)  #设置标题
    ax.set_xlim(set_xlim_l,set_xlim_r)  #设置x轴数据限值
    plt.savefig(output)
    plt.show()  #显示图像
    
    
def Acchistogram(input,output):

    data=pd.read_csv(input)
    
    key=[int(i) for i in data.columns]  #年份从header中提取
    value=hot_dog.T.values   #将冠亚季军所吃热狗的数量转化成matrix，也就是[[25,24,22],[50.0,31.0,23.5],...]
    v1=[i[0]+i[1]+i[2] for i in value]  #第一次画的柱形图y值为冠亚季军所吃热狗数量的总和
    v2=[i[1]+i[2] for i in value]  #第二次画的柱形图y值为亚军所吃热狗的数量+季军所吃热狗的数量
    v3=[i[2] for i in value]  #第三次画的柱形图y值为季军所吃热狗的数量
    
    ax.bar(year,v1,color="green")
    ax.bar(year,v2,color="red")
    ax.bar(year,v3,color="blue")
    ax.set(xlabel="Year",title="Hotdog game scores 2000-2010")
    ax.text(1998,184,"(HDB)")  #设置文字
    ax.legend(["first place","second place","third place"])  #设置图例
    plt.show()     