"""Visualize the slices of Appositionizer
"""
import glob
from itertools import product
import math
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import seaborn as sns

try:
    filename, outname = sys.argv[1:]
except ValueError:
    print(f"usage: {os.path.basename(sys.argv[0])} filename outputname")
    sys.exit(1)

names = 'xmin xmax ymin ymax zmin zmax'.split()
df = pd.read_csv(filename, sep=' ', header=None, names=names, comment='#')

size = (len(df) - 1) // 2

universe = df.iloc[0]
slices = df.iloc[1:1 + size]
boxes = df.iloc[1 + size:]


def draw_boxes(ax, df, color):
    """Draw boxes from a dataframe
    """
    for _, (xmin, xmax, ymin, ymax) in df.iterrows():
        r = Rectangle((xmin, ymin),
                      xmax - xmin,
                      ymax - ymin,
                      edgecolor=color,
                      fill=False,
                      linewidth=1)
        ax.add_patch(r)
    return r


def tack(a, b, scale: float) -> float:
    """Convert a number to a point on a grid between to limits
    """
    start = scale * (math.floor(a / scale) + 0.5)
    end = scale * (math.floor(b / scale) + 1)
    while start <= end:
        yield start
        start += scale


def pixelize(df, cols, scale: float=50.):
    """Convert bounding boxes to points

    :param scale: size of the grid to pizelize on
    """
    for _, row in df.iterrows():
        limits = [tack(getattr(row, col + 'min'),
                       getattr(row, col + 'max'),
                       scale) for col in cols]
        for x, y in product(*limits):
            yield x, y


def draw_heat(ax, df, cols, scale: float=50.):
    """Draw a heatmap of pixelized boxes
    """
    data = pd.DataFrame(pixelize(df, cols, scale), columns=cols)
    xs = getattr(data, cols[0])
    xbins = xs.unique() - 0.5 * scale
    xbins.sort()
    ys = getattr(data, cols[1])
    ybins = ys.unique() - 0.5 * scale
    ybins.sort()
    plt.hist2d(xs, ys, bins=[xbins, ybins], cmap='Blues')
    cb = plt.colorbar()
    # sns.kdeplot(, getattr(data, cols[1]), shade=True, ax=ax)
    return

    data = pd.DataFrame(pixelize(df, cols), columns=cols).groupby(cols).size()
    print(data)
    # print(data.to_frame().pivot(*cols))
    data = data.to_frame(name='v').reset_index()
    data = data.pivot('x', 'z', 'v').fillna(0)
    sns.heatmap(data, ax=ax)


sns.set()

ax = plt.subplot(111)
cols = ['x', 'z']

cpubox = draw_heat(ax, boxes, cols)

draw = []
for c, axis in zip(cols, ('x', 'y')):
    extends = []
    for extend in ('min', 'max'):
        draw.append(c + extend)
        col = getattr(df, c + extend)
        extends.append(getattr(col, extend)())
    getattr(ax, f'set_{axis}lim')(extends[0] * 1.05, extends[1] * 1.05)

cpubox = draw_boxes(ax, boxes[draw], 'gray')
slither = draw_boxes(ax, slices[draw], 'tomato')

coords = None
dirname = os.path.dirname(filename)
for fn in glob.glob(f'{dirname}/coords.*.txt'):
    data = pd.read_csv(fn, sep=' ', header=None, names='xyz')
    if coords is None:
        coords = data
    else:
        coords = coords.append(data)
if coords is not None:
    coords.plot.scatter(cols[0], cols[1], ax=ax, s=1)

ax.legend((slither, cpubox), ('Slice', 'CPU realm'))

for method, label in zip(('xlabel', 'ylabel'), cols):
    getattr(plt, method)(label)

plt.savefig(outname)
