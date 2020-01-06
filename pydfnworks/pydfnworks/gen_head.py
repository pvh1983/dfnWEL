# Copyright 2019, Nevada System of Higher Education on Behalf of 
# the Desert Research Institute. All rights reserved.

# This is free software; you can redistribute it and/or modify 
# it under the terms of the GNU Lesser General Public License as 
# published by the Free Software Foundation; either version 3.0 
# of the License, or (at your option) any later version.

import math
import numpy as np
import vtk
import os
import sys
from array import *

# Calculate groundwater levels from pressures
#rho = 997.32
#gra = 9.81
#ele_head = 1094.8    # Bottom of the model domain
#write_head=0
#domain_thickness=100
def hhead(rho,gra,ele_head,domain_thickness,write_head):

	# Step 1: Set up the reader
	reader = vtk.vtkDataSetReader()

	runid=os.environ['runid']
	print 'runid = ', runid
	fname="parsed_vtk/iflow_"+runid+"-001.vtk"
	print "pressure file read: " + fname
	reader.SetFileName(fname)
	#reader.SetFileName("uGridEx.vtk")
	#reader.ReadAllVectorsOn()
	reader.ReadAllScalarsOn() # Activate the reading of all scalars
	reader.Update()

	# The reader is the top level object we use to access a VTK file.
	# For example, you can know the header of a VTK file by
	reader.GetHeader()

	# Step 2: Get data in DATASET block.
	data=reader.GetOutput()
	npoints=data.GetNumberOfPoints() # to know how may Points
	ncells = data.GetNumberOfCells()
	nscalars=reader.GetNumberOfScalarsInFile()
	nvectors=reader.GetNumberOfVectorsInFile()
	ntensors=reader.GetNumberOfTensorsInFile()
	#dim = data.GetDimensions()
	x = np.zeros(data.GetNumberOfPoints())
	y = np.zeros(data.GetNumberOfPoints())
	z = np.zeros(data.GetNumberOfPoints())

	for i in range(data.GetNumberOfPoints()):
			x[i],y[i],z[i] = data.GetPoint(i)

	# print "Number of points: " ,  npoints
	# print "Number of cells: " ,  ncells
	# print "Number of scalars: " ,  nscalars
	# print "Number of tensors: " ,  ntensors
	# print "Number of vectors: " ,  nvectors
	# for i in xrange(nscalars):
	#   Scalarnames=reader.GetScalarsNameInFile(i)
	#   print "Scalarnames",i,Scalarnames
	# for i in xrange(nvectors):
	#   Vectornames=reader.GetVectorsNameInFile(i)
	#   print "Vectornames",i,Vectornames

	# Step 3: get data in POINT_DATA/CELL_DATA block.
	d=data.GetPointData()
	array = d.GetArray('Liquid_Pressure')
	#print "Pressure at node 0 is", array.GetValue(0)
	#print array


	# Calculate groundwater levels from pressures
#	rho = 997.32
#	gra = 9.81
#	ele_head = 1094.8    # Bottom of the model domain

#	fid = open('xyzH.dat', 'w')
	head = np.zeros((npoints, 1), 'float')
	pressure = np.zeros((npoints, 1), 'float')
	for i in range(npoints):
		pres = array.GetValue(i)
		pressure[i,0]=pres
		head[i,0] = pres/rho/gra+z[i]+domain_thickness/2+ele_head
#		fid.write('%f %f %f %f %f \n' % (x[i], y[i], z[i],pres,head[i,0]))
#	fid.close()

	print('--> Finished writing simulated head to file')
	print('== Statistics of heads and pressures =====')
	print 'Min  = ', head.min(),  pressure.min()
	print 'Max  = ', head.max(),  pressure.max()
	print 'Mean = ', head.mean(), pressure.mean()
	print 'Std  = ', head.std(),  pressure.std()
	print 'Var  = ', head.var(),  pressure.var()

	print 'Location  Hmin | Hmax | Hmean | Hvar |  Pmin | Pmax | Pmean | Pvar'
	print '*', head.min(), head.max(), head.mean(), head.std(), head.var(), pressure.min(), pressure.max(), pressure.mean(), pressure.std(), pressure.var(), '    (*: entire_domain)' 

	# print head at observation locations
	fid2 = open('head_stas.out', 'w')
	fid2.write('obs,  Hmin,      Hmax,         Hmean,        H1std,     Hvar,     Pmin,      Pmax,         Pmean,        P1std,     Pvar \n')
	fid2.write('%d, %8.1f, %8.1f, %8.1f, %12.4f, %12.4f, %9.1f, %9.1f, %9.1f, %12.4f, %12.4f \n' % (0,head.min(), head.max(), head.mean(), head.std(), head.var(), pressure.min(), pressure.max(), pressure.mean(), pressure.std(), pressure.var()))
	for ifile in range(3):
		if ifile == 0:
			data=np.loadtxt('well_location_001.ex', delimiter=None, skiprows=1)
			ofile_H='H_ER2061.dat';
		elif ifile == 1:
			data=np.loadtxt('well_location_002.ex', delimiter=None, skiprows=1)          
			ofile_H='H_ER2062.dat';
		else:
			data=np.loadtxt('well_location_003.ex', delimiter=None, skiprows=1)
			ofile_H='H_ER2063.dat';
			#    print 'Size of vector = ', len(data[:,0]), type(data)
		node_id = data[:,0].astype(int) 
		Hsim=head[node_id,:]
		Psim=pressure[node_id,:]
		np.savetxt(ofile_H, Hsim,fmt='%.6e', delimiter=' ')  # 
		print ifile+1, Hsim.min(), Hsim.max(), Hsim.mean(), Hsim.std(), Hsim.var(), Psim.min(), Psim.max(), Psim.mean(), Psim.std(), Psim.var()		
		fid2.write('%d, %8.1f, %8.1f, %8.1f, %12.4f, %12.4f, %9.1f, %9.1f, %9.1f, %12.4f, %12.4f \n' % (ifile+1,Hsim.min(), Hsim.max(), Hsim.mean(), Hsim.std(), Hsim.var(), Psim.min(), Psim.max(), Psim.mean(), Psim.std(), Psim.var()))
		

	fid2.close()

	# Step xxx: Write a new vtk file

	if write_head==1:
		rmder = npoints%10 # get remainder
		head1 = head[:(npoints-rmder)]
		head1_rsh = np.reshape(head1,(npoints/10,10))
		head2 = head[(npoints-rmder):]
		head2_rsh = np.reshape(head2,(1,rmder))
		grid_file = "full_mesh.vtk"
		with open(grid_file, 'r') as f:
			grid = f.readlines()[3:]
		for line in grid:
			if 'POINTS' in line:
				num_cells = line.strip(' ').split()[1]
		header = ['# vtk DataFile Version 2.0 [hpham]\n',
				  'PFLOTRAN output\n',
				  'ASCII\n']
		#os.remove("head.vtk")
		#filename = "head.vtk"
		filename="head_run_"+runid+".vtk"
		with open(filename, 'w') as f:
			for line in header:
				f.write(line)
			for line in grid:
				f.write(line)
			f.write('\n')
			f.write('POINT_DATA ')
			f.write('%s \n' % (num_cells))
			f.write('SCALARS Hydraulic_Head_m      float 1 \n')
			f.write('LOOKUP_TABLE default \n')
		# Print numpy array to file
		filename2 = open(filename, "a")
		np.savetxt(filename2, head1_rsh,fmt='%.4e', delimiter=' ')  # X is an array
		np.savetxt(filename2, head2_rsh,fmt='%.4e', delimiter=' ')  # X is an array
		filename2.close()
		print "Done converting pressure to hydraulic head and writing data to head.vtk!"

	### Print aperture stats

# Copyright 2019, Nevada System of Higher Education on Behalf of 
# the Desert Research Institute. All rights reserved.

# This is free software; you can redistribute it and/or modify 
# it under the terms of the GNU Lesser General Public License as 
# published by the Free Software Foundation; either version 3.0 
# of the License, or (at your option) any later version.