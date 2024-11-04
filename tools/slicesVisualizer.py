##### TOOL TO DISPLAY THE SLICES OUTPUT AT output/Slices.N.txt
##### USAGE: mayavi2 -x ./slicesVisualizer.py output/Slices.N.txt

import sys
import numpy as np
from tvtk.api import tvtk
from mayavi.scripts import mayavi2
from numpy import array


# Uncomment the next two lines to save the dataset to a VTK XML file.
#w = tvtk.XMLPolyDataWriter(input=mesh, file_name='polydata.vtp')
#w.write()

file_gap = 9; #ie read 1 in every 'file_gap' rows

minCoords = []
maxCoords = []
slicePoints = []
sliceLines = []
sliceTriangles = []

def read_Histogram_Cuts():
    global sliceTriangles
    global slicePoints
   
    file = open(sys.argv[3],"r")

    #ignore line 1
    arr_text = file.readline()

    #read universe box from line 2
    arr_text = file.readline().split(' ')
    x0 = float(arr_text[0])
    y0 = float(arr_text[1])
    z0 = float(arr_text[2])
    x1 = float(arr_text[3])
    y1 = float(arr_text[4])
    z1 = float(arr_text[5])
    print "Universe: [" + str(x0) + ", " + str(y0) + ", " + str(z0)+ " ; " + str(x1) + ", " + str(y1) + ", " + str(z1)+ "]"

    #ignore line 3
    arr_text = file.readline()
    #get slices configuration from line 4
    arr_text = file.readline().split(' ')
    nb_subvolumes = int(arr_text[0])
    print "Slices: " + str(nb_subvolumes) + " " + arr_text[1]

    #ignore line 5
    arr_text = file.readline()
 
    #read all N slices
    index=0;
    for i in range (0,nb_subvolumes):
	arr_text = file.readline().split(' ')
	
        #skip rows as stated by file_gap
        if i % file_gap > 0:
            continue;

	#1st corner
	x0 = float(arr_text[0])
	y0 = float(arr_text[1])
	z0 = float(arr_text[2])
	
	#2nd corner
	x1 = float(arr_text[3])
	y1 = float(arr_text[4])
	z1 = float(arr_text[5])

	#8 points
	slicePoints.append([x0,y0,z0])
	slicePoints.append([x0,y0,z1])
	slicePoints.append([x1,y0,z1])
	slicePoints.append([x1,y0,z0])
	slicePoints.append([x0,y1,z0])
	slicePoints.append([x0,y1,z1])
	slicePoints.append([x1,y1,z1])
	slicePoints.append([x1,y1,z0])

	#y=0 plane
	sliceTriangles.append([index  , index+1, index+2])
	sliceTriangles.append([index  , index+2, index+3])

	#y=1 plane
	sliceTriangles.append([index+4, index+5, index+6])
	sliceTriangles.append([index+4, index+6, index+7])

	#x=0 plane
	sliceTriangles.append([index  , index+1, index+5])
	sliceTriangles.append([index  , index+5, index+4])

	#x=1 plane
	sliceTriangles.append([index+3, index+2, index+6])
	sliceTriangles.append([index+3, index+6, index+7])

	#z=0 plane
	sliceTriangles.append([index  , index+3, index+7])
	sliceTriangles.append([index  , index+7, index+4])

	#z=1 plane
	sliceTriangles.append([index+1, index+2, index+6])
	sliceTriangles.append([index+1, index+6, index+5])

	#slice lines
	sliceLines.append([index+0, index+1])
	sliceLines.append([index+1, index+2])
	sliceLines.append([index+2, index+3])
	sliceLines.append([index+3, index+0])

	sliceLines.append([index+4, index+5])
	sliceLines.append([index+5, index+6])
	sliceLines.append([index+6, index+7])
	sliceLines.append([index+7, index+4])

	sliceLines.append([index+0, index+4])
	sliceLines.append([index+1, index+5])
	sliceLines.append([index+2, index+6])
	sliceLines.append([index+3, index+7])

	index = index+8;
    file.close()
    #print slicePoints
    #print sliceTriangles

@mayavi2.standalone
def view():
    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi.modules.outline import Outline
    from mayavi.modules.surface import Surface
    from mayavi.modules.axes import Axes
    from mayavi.modules.orientation_axes import OrientationAxes
#    from mayavi.modules.hyper_streamline import HyperStreamline

    # OBJECT GROUP FOR BOUNDARIES -  The numpy array data.
    read_Histogram_Cuts()
    points = array (slicePoints)
    triangles = array (sliceTriangles)
    straightLines = array (sliceLines)
#   temperature = array([100., 120., 130., 140.])

    # The TVTK dataset.
    mesh = tvtk.PolyData(points=points, polys=triangles, lines=straightLines)
#    mesh.point_data.scalars = temperature
#    mesh.point_data.scalars.name = 'Temperature'
 
    src = VTKDataSource(data = mesh)
    mayavi.add_source(src)
    s = Surface()
    s.actor.property.opacity = 0.5
    s.actor.property.color = (0.4, 0.0, 0.8)
    s.actor.property.representation = 0 #points: 0, wireframe: 1, surface: 2

    a = Axes()
    a.axes.number_of_labels=11
    a.axes.font_factor=0.7
    mayavi.add_module(a)

    mayavi.add_module(s)
    mayavi.add_module(Outline())
    mayavi.add_module(OrientationAxes())

    print "Rendering shape in mayavi..."

#    mayavi.add_module(HyperStreamline())

if __name__ == '__main__':

    if (len(sys.argv) == 4):
        print "Program started, slices file: " + sys.argv[3]
        view()
    else:
        print "USAGE: mayavi2 -x ./slicesVisualizer.py output/Slices.N.txt"
        sys.exit()
