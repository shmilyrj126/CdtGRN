# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 11:19:33 2019

@author: hsc
"""

import matplotlib.pyplot as plt
import numpy as np

y=np.loadtxt("C:/Users/hsc/Desktop/liver1.csv",
               delimiter=',',
               unpack=True)

       
x=[4,13]   


fig  = plt.figure()

ax = fig.add_subplot(1,1,1)

ax.plot(x,y)
    
ax.set_title("Correlation")
ax.set_xlabel("x-week")
ax.set_ylabel("expression")
plt.show()

print(y)
