# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 19:47:55 2019

@author: hsc
"""

# encoding=utf-8
import matplotlib.pyplot as plt
from pylab import *                 #支持中文
mpl.rcParams['font.sans-serif'] = ['SimHei']

#import openpyxl
 
# 打开excel文件,获取工作簿对象
#wb = openpyxl.load_workbook(u'C:/Users/hsc/Desktop/liver.xlsx')
# 从表单中获取单元格的内容
#ws = wb.active  # 当前活跃的表单
 
# 从表单中获取行和列
 
#colC = ws['C'] # 获取整个C列
#print(colC)
#APOL1 = ws[2]   # 获取第6行
#ATM=ws[3]
 
time = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
x = range(len(time))
A = [108.9,70.6,76,78.5,71.9,65.8,74.8,93.2,80.9,78.6,137.1,244.9,88.5,43.9,79.6,80.9,43.7,102.8,204.6,182,123.2]
B = [79.9,91.9,55.7,66.7,104,86.1,79.8,85.7,72.1,66.1,103.2,174,53,56.1,81.9,132.6,55.8,156.7,191.2,210.3,140.5]
#C=[134.1,253.8,279.8,269.3,379.9,233,254.7,153,223.2,177.1,365.4,283.8,247.5,167.1,105.2,116.4,272.9,430.3,276.5,228.1,278.8]


plt.plot(x, A, marker='o', color='r',mec='r', mfc='w',label=u'A')
#可用MarkerEdgeColor或mec设置标记边缘颜色；MarkerFaceColor或mfc设置标记填充颜色;MarkerSize设置标记大小
 # r：红色 w：白色 g：绿色  b：蓝色  c：青绿色  m：洋红色  k：黑色
plt.plot(x, B, marker='o',color='g', mec='g', mfc='w',label=u'B')
#plt.plot(x, C, marker='o',color='b', mec='b', mfc='w',label=u'C')


plt.legend() # 让图例生效
plt.xticks(x, time, rotation=0)
plt.margins(0.1)#居中显示
plt.subplots_adjust(bottom=0)#调整子图布局
plt.xlabel(u"time(week)") #X轴标签
plt.ylabel("GeneExpression") #Y轴标签
plt.title("correlation") #标题
 
plt.show()