# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 19:14:20 2019

@author: hsc
"""

import xlrd

temp=xlrd.open_workbook('C:\\Users\\hsc\\Desktop\\liver.xlsx')

table = temp.sheet()[0]

x=table.row_values(0)

APOL1=x[1]
   