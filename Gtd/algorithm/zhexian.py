# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 19:17:00 2019

@author: hsc
"""

# encoding=utf-8
import matplotlib.pyplot as plt
from pylab import *                 #支持中文
mpl.rcParams['font.sans-serif'] = ['SimHei']
 
time = [4,5,6,7,8,9,10,11,12,13,16,18,19]
x = range(len(time))
APOL1 = [80,70,59,66,51,52,55,65,71,62,75,69,72]
ATM=[88,71,50,60,35,50,46,40,64,55,57,56,47]


plt.plot(x, APOL1, marker='o', color='r',mec='r', mfc='w',label=u'ras')
#可用MarkerEdgeColor或mec设置标记边缘颜色；MarkerFaceColor或mfc设置标记填充颜色;MarkerSize设置标记大小
 # r：红色 w：白色 g：绿色  b：蓝色  c：青绿色  m：洋红色  k：黑色
plt.plot(x, ATM, marker='o',color='g', mec='g', mfc='w',label=u'GBP')


plt.legend() # 让图例生效
plt.xticks(x, time, rotation=0)
plt.margins(0.1)#居中显示
plt.subplots_adjust(bottom=0)#调整子图布局
plt.xlabel(u"time(week)") #X轴标签
plt.ylabel("GeneExpression") #Y轴标签
plt.title("correlation") #标题
 
plt.show()