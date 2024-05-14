#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:20:50 2024

@author: heitor
"""


#-----------


import bezier
import numpy as np
import matplotlib.pyplot as plt


#-----------




nodes = np.asfortranarray([
    [0.0, 0.625, 1.0],
    [0.0, 0.5  , 0.5],
])


curve = bezier.Curve(nodes, degree=2)

curve.evaluate(0.75)


s_vals = np.linspace(0.0, 1.0, 10)

curve.evaluate_multi(s_vals)



#plot
f, ax = plt.subplots(figsize=(10, 5))


ax.plot(curve.evaluate_multi(s_vals)[0],curve.evaluate_multi(s_vals)[1])

ax.plot(curve.evaluate(0.75)[0],curve.evaluate(0.75)[1],'o')

ax.plot(nodes[0],nodes[1],'o')









































































#