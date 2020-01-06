# Copyright 2019, Nevada System of Higher Education on Behalf of 
# the Desert Research Institute. All rights reserved.

# This is free software; you can redistribute it and/or modify 
# it under the terms of the GNU Lesser General Public License as 
# published by the Free Software Foundation; either version 3.0 
# of the License, or (at your option) any later version.

#!/usr/bin/env python

# uniform particle distribution on fractures in sphere
# Version005 Last update: 03/20/2018 - Fix all locs at a same cell.
# Version006 Last update: 03/22/2018 - Read smaller file size full_mesh.vtk

import sys
import os
from time import time
import numpy as np


lagritpath = '/projects/DFN/apps/LAGRIT/lagrit_ulin3.2'

#sphere_ic.ic("full_mesh.inp",1000,1,(0,0,0))

main_time = time()

def ic(filename,n_part,sphere_r,sphere_center):
#filename = "full_mesh.vtk"
#n_part = 10000
#sphere_r = 50
#sphere_center = (0,0,0)

	main_elapsed0 = time() - main_time
	timing0 = 'Total time required for loading full_mesh.inp: %0.2f Minutes'%(main_elapsed0/60.0)
	print timing0

	#### load in node locations and connectivity table
	f = open(filename,'r') #file
	info_str = f.read() #string containing node locations
	f.close() #close node location file
	info_str = info_str.split('\n') #split into rows
	g = info_str[4]
	g = g.split()
	n_nodes = int(g[1])
	pos = np.array([[0 for j in range(3)] for i in range(1,n_nodes+1)], dtype=np.dtype(float))
		#x,y, and z positions of nodes
	for i in range(5,n_nodes+5):
		ln = info_str[i] #read node locations
		ln = ln.split() #split into columns
		for j in range(0,3):
			num = ln[j] #the i,jth element of pos
			pos[i-5][j] = float(num) #convert to float
	ect = np.array([[0 for j in range(3)] for i in range(n_nodes+6,len(info_str)-2)]) #ect
	for i in range(n_nodes+6,len(info_str)-2):
		ln = info_str[i] #read ect
		ln = ln.split()
		for j in range(1,4):
			num = ln[j]
			ect[i-n_nodes-6][j-1] = int(num)
	ect -= 1 #move to base 0
	del info_str #delete variable info_str (too big)

	#### place particles
	x_node = pos[0:,0] #x positions of nodes
	y_node = pos[0:,1] #y positions of nodes
	z_node = pos[0:,2] #z positions of nodes
	M = np.vstack((x_node-sphere_center[0],y_node-sphere_center[1],z_node-sphere_center[2])) #matrix of node locations relative to sphere center
	r = np.sqrt(np.sum(np.multiply(M,M),axis=0)) #distance between node and sphere center
	x_node[r>sphere_r] = np.nan #mark nodes that are outside of sphere
	x_el = x_node[ect] #x positions of element corners
	y_el = y_node[ect] #y positions of element corners
	z_el = z_node[ect] #z positions of element corners
	remove = np.isnan(x_el)
	remove = np.sum(remove.astype(int),axis=1)
	del_ind = remove.nonzero()
	x_el = np.delete(x_el,del_ind,axis=0)
	y_el = np.delete(y_el,del_ind,axis=0)
	z_el = np.delete(z_el,del_ind,axis=0)
	n_el = len(z_el) #number of elements
	A = np.zeros((n_el),dtype=np.dtype(float)) #area of each element
	np.seterr(invalid='ignore') #supress warning about invalid value in det
	for i in range(n_el):
		M = np.vstack([x_el[i,0:2]-x_el[i,2], y_el[i,0:2]-y_el[i,2], z_el[i,0:2]-z_el[i,2]]) #matrix used to compute area of element
		A[i] = .5*np.sqrt(np.linalg.det(np.matmul(M.transpose(),M))) #area of element i
		
	del_ind = (np.isnan(A)).nonzero() #indices where elements are outside of sphere
	A = np.delete(A,del_ind[0]) #only use elements in sphere
	x_el = np.delete(x_el,del_ind[0],0)
	y_el = np.delete(y_el,del_ind[0],0)
	z_el = np.delete(z_el,del_ind[0],0)
	n_el = len(z_el)
	cum_A = A.cumsum(axis=0) #cumulative sum of element areas
	cum_A = cum_A/cum_A[-1] #normalized
	x_part = np.zeros((n_part),dtype=np.dtype(float)) #x position of particles
	y_part = np.zeros((n_part),dtype=np.dtype(float)) #y position of particles
	z_part = np.zeros((n_part),dtype=np.dtype(float)) #z position of particles
	el_part = np.zeros((n_part),dtype=np.dtype(int)) #elements particles are in
	n_p_in_el = np.round(A/sum(A)*n_part) #approximate number of particles in each element
	cum_n_p_in_el = n_p_in_el.cumsum(axis=0) #cumulative minimum number of particles in elements
	cum_n_p_in_el = cum_n_p_in_el.astype(int) #convert to integer

	it = 0
	for i in range(n_el):
		el_part[it:cum_n_p_in_el[i]] = i
		it = cum_n_p_in_el[i]+1
	U = np.random.rand(n_part-it+1)
	for i in range(it,n_part):
		ind = (cum_A>=U[i-it]).nonzero() #indices of cum_A that are at least U
		ind = ind[0]
		el_part[i] = ind[0] #assign first element that has cum_A at least U

	U = np.random.rand(n_part,2)
	for i in range(n_part):
		r_1 = np.sqrt(U[i][0])
		r_2 = U[i][1]
		x = x_el[el_part[i],0:] #x position of nodes in element el
		y = y_el[el_part[i],0:] #y position of nodes in element el
		z = z_el[el_part[i],0:] #z position of nodes in element el
		x_part[i] = (1-r_1)*x[0] + r_1*(1-r_2)*x[1] + r_1*r_2*x[2] #convert particle x position to Cartesian coordinates
		y_part[i] = (1-r_1)*y[0] + r_1*(1-r_2)*y[1] + r_1*r_2*y[2] #convert particle y position to Cartesian coordinates
		z_part[i] = (1-r_1)*z[0] + r_1*(1-r_2)*z[1] + r_1*r_2*z[2] #convert particle z position to Cartesian coordinates



	## save
	f = open('ParticleInitCoordR.dat','w')
	#f.write('X,Y,Z,ele_ID,frac_ID' + '\n')
	f.write('npart ' + str(n_part) + '\n')
	for j in range(n_part):
	#    f.write(str(x_part[j]) + ' ' + str(y_part[j]) + ' ' + str(z_part[j]) + ' ' + str(el_part[j]) + ' ' + str(fract_part[j]) + '\n')
	    f.write(str(x_part[j]) + ' ' + str(y_part[j]) + ' ' + str(z_part[j]) + '\n')
	f.close()


	main_elapsed1 = time() - main_time
	timing1 = 'Total time required for placing particles: %0.2f Minutes'%(main_elapsed1/60.0)
	print timing1

	f=open('definedist.lgi','w')
	f.write("read / avs / "+str("full_mesh.inp")+" / mo2 \n")
	f.write("\n")
	f.write('cmo / create/ mo1  \n')
	f.write('cmo/readatt/ mo1/ npar / 1,0,0 / ParticleInitCoordR.dat  \n')
	f.write('cmo/readatt/ mo1/ xic, yic, zic/ 1,0,0 / ParticleInitCoordR.dat  \n')
	f.write('cmo / set_id/mo2 / node/ n_num \n')
	f.write('cmo / addatt/mo1 /idnum/ VINT/scalar / nnodes  \n')
	f.write('interpolate/voronoi/mo1,idnum/1 0 0/mo2,n_num  \n')
	f.write('dump / avs2/ClosestNodeR.inp / mo1/0 0 1 0  \n')

	f.write("finish \n")
	f.write("\n")
	f.flush()
	f.close()

	os.system(lagritpath + " <definedist.lgi >distance.out")


	main_elapsed2 = time() - main_time
	timing2 = 'Total time required for placing particles+finding nearby nodes: %0.2f Minutes'%(main_elapsed2/60.0)
	print timing2
	print 'End of file sphere_ic.py.'
	print ' '

# Copyright 2019, Nevada System of Higher Education on Behalf of 
# the Desert Research Institute. All rights reserved.

# This is free software; you can redistribute it and/or modify 
# it under the terms of the GNU Lesser General Public License as 
# published by the Free Software Foundation; either version 3.0 
# of the License, or (at your option) any later version.