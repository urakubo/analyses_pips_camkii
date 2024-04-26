
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot, patches

import utils
import parameters as p
import colormap as c

x1 = [0, 1,2,3,4,5,6]
x2 = [2,4,6,8,10,12]
X1, X2 = np.meshgrid(x1, x2)


dir_edited_data, prefix, suffix = 'data3/valency_length', 'tension', 'cluster'
d   = utils.load(dir_edited_data, prefix, suffix)

angles_interface = d['angles_interface']
angles_all = d['angles_all']
cluster_coefficients = d['cluster_coefficients']

angles_interface = angles_interface.reshape(6,7)
angles_all = angles_all.reshape(6,7)
cluster_coefficients = cluster_coefficients.reshape(6,7)


'''
fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X1, X2, angles_interface, cmap='bwr', linewidth=0)
plt.show()

fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X1, X2, cluster_coefficients, cmap='bwr', linewidth=0)
plt.show()
'''

'''
fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X1, X2, angles_all, cmap='bwr', linewidth=0)
plt.show()
'''

fig = plt.figure()
ax  = fig.add_subplot(111)
surf = ax.pcolormesh(X1, X2, angles_all, cmap='bwr', linewidth=0)
plt.show()


