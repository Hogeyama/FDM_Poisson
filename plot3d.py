#!/usr/bin/env python3.5

import csv
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Util
def chunks(l, n):
    return [l[i:i + n] for i in range(0, len(l), n)]

def cull(l,n):
    return [l[i] for i in range(0, len(l), n)]

def unzip(l):
    return tuple(map(list, zip(*l)))

def readCSV(filepath):
    """
    x1,y1,z1,...
    x2,y2,z2,... => ([x1,x2,...],[y1,y2,...],[z1,z2,...],...)
    x3,y3,z3,...
    ...
    """
    with open(filepath, newline='') as csvfile:
        return unzip(list(map(lambda w: list(map(float,w)), csv.reader(csvfile, delimiter=','))))

def to2D(x):
    n = int(math.sqrt(len(x)))
    return chunks(x,n)

def to2Dcull(x,n):
    return cull(list(map(lambda w: cull(w,n), to2D(x))),n)

# Plot
fig = plt.figure()
xs,ys,zs,ws = map(lambda w: to2Dcull(w,2) , readCSV('./plot/plot.csv'))
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(xs, ys, ws)
ax.scatter(xs, ys, zs)

plt.show()

