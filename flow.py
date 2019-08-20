import os 
import subprocess
import sys
import helper
import glob
import shutil
from time import time 
import numpy as np
import h5py
import vtk     #hpham
def dfn_flow(self):
    ''' dfnFlow
    Run the dfnFlow portion of the workflow.
    ''' 

    print('='*80)
    print("\ndfnFlow Starting\n")
    print('='*80)

    tic_flow = time()

    tic = time()
    self.lagrit2pflotran()
    helper.dump_time(self._jobname, 'Function: lagrit2pflotran', time() - tic)   
    
    tic = time()    
    self.pflotran()
    helper.dump_time(self._jobname, 'Function: pflotran', time() - tic)  

    tic = time()    
    self.parse_pflotran_vtk_python()
    helper.dump_time(self._jobname, 'Function: parse_pflotran_vtk', time() - tic)    

    tic = time()    
    self.pflotran_cleanup()
    helper.dump_time(self._jobname, 'Function: parse_cleanup', time() - tic) 
    helper.dump_time(self._jobname,'Process: dfnFlow',time() - tic_flow)    

    print('='*80)
    print("\ndfnFlow Complete\n")
    print('='*80)
       
def lagrit2pflotran(self, inp_file='', mesh_type='', hex2tet=False):
    """  Takes output from LaGriT and processes it for use in PFLOTRAN.
    
    Kwargs:
        * inp_file (str): name of the inp (AVS) file produced by LaGriT 
        * mesh_type (str): the type of mesh
        * hex2tet (boolean): True if hex mesh elements should be converted to tet elements, False otherwise.
    """
    print ('='*80)
    print os.path.abspath(os.getcwd())     # hpham
    print("Starting conversion of files for PFLOTRAN ")
    print ('='*80)
    if inp_file:
        self._inp_file = inp_file
    else:
        inp_file = self._inp_file

    if inp_file == '':
        sys.exit('ERROR: Please provide inp filename!')

    if mesh_type:
        if mesh_type in mesh_types_allowed:
            self._mesh_type = mesh_type
        else:
            sys.exit('ERROR: Unknown mesh type. Select one of dfn, volume or mixed!')
    else:
        mesh_type = self._mesh_type

    if mesh_type == '':
        sys.exit('ERROR: Please provide mesh type!')

    self._uge_file = inp_file[:-4] + '.uge'
    # Check if UGE file was created by LaGriT, if it does not exists, exit
    failure = os.path.isfile(self._uge_file)
    if failure == False:
        sys.exit('Failed to run LaGrit to get initial .uge file')

    if mesh_type == 'dfn':
        self.write_perms_and_correct_volumes_areas() # Make sure perm and aper files are specified

    # Convert zone files to ex format
    #self.zone2ex(zone_file='boundary_back_n.zone',face='north')
    #self.zone2ex(zone_file='boundary_front_s.zone',face='south')
    #self.zone2ex(zone_file='boundary_left_w.zone',face='west')
    #self.zone2ex(zone_file='boundary_right_e.zone',face='east')
    #self.zone2ex(zone_file='boundary_top.zone',face='top')
    #self.zone2ex(zone_file='boundary_bottom.zone',face='bottom')
    self.zone2ex(zone_file='all')
    print ('='*80)
    print("Conversion of files for PFLOTRAN complete")
    print ('='*80)
    print("\n\n")

def zone2ex(self, uge_file='', zone_file='', face=''):
    '''zone2ex    
    Convert zone files from LaGriT into ex format for LaGriT
    inputs:
    uge_file: name of uge file
    zone_file: name of zone file
    face: face of the plane corresponding to the zone file

    zone_file='all' processes all directions, top, bottom, left, right, front, back
    '''

    print('--> Converting zone files to ex')    
    if self._uge_file:
        uge_file = self._uge_file
    else:
        self._uge_file = uge_file

    uge_file = self._uge_file
    if uge_file == '':
        sys.exit('ERROR: Please provide uge filename!')
    # Opening uge file
    print('\n--> Opening uge file')
    fuge = open(uge_file, 'r')

    # Reading cell ids, cells centers and cell volumes
    line = fuge.readline()
    line = line.split()
    NumCells = int(line[1])

    Cell_id = np.zeros(NumCells, 'int')
    Cell_coord = np.zeros((NumCells, 3), 'float')
    Cell_vol = np.zeros(NumCells, 'float')

    for cells in range(NumCells):
        line = fuge.readline()
        line = line.split()
        Cell_id[cells] = int(line.pop(0))
        line = [float(id) for id in line]
        Cell_vol[cells] = line.pop(3)
        Cell_coord[cells] = line
    fuge.close()

    print('--> Finished with uge file\n')

    # loop through zone files   hpham modified this.
    if zone_file is 'all':
            zone_files = ['pboundary_front_s.zone', 'pboundary_back_n.zone', 'pboundary_left_w.zone', \
                            'pboundary_right_e.zone', 'pboundary_top.zone', 'pboundary_bottom.zone', \
                            'pointset_well1.zone','pointset_well2.zone','pointset_well3.zone', \
                            'pointset_well12.zone','pointset_well22.zone','pointset_well32.zone']
            face_names = ['south', 'north', 'west', 'east', 'top', 'bottom','well1','well2','well3','well12','well22','well32']
    else: 
            if zone_file == '':
                sys.exit('ERROR: Please provide boundary zone filename!')
            if face == '':
                sys.exit('ERROR: Please provide face name among: top, bottom, north, south, east, west, wells !')
            zone_files = [zone_file]
            face_names = [face]
            
    for iface,zone_file in enumerate(zone_files):
            face = face_names[iface]
            # Ex filename
            ex_file = zone_file.strip('zone') + 'ex'

            # Opening the input file
            print '--> Opening zone file: ', zone_file
            fzone = open(zone_file, 'r')
            fzone.readline()
            fzone.readline()
            fzone.readline()

            # Read number of boundary nodes
            print('--> Calculating number of nodes')
            NumNodes = int(fzone.readline())
            Node_array = np.zeros(NumNodes, 'int')
            # Read the boundary node ids
            print('--> Reading boundary node ids')

            if (NumNodes < 10):
                g = fzone.readline()
                node_array = g.split()
                # Convert string to integer array
                node_array = [int(id) for id in node_array]
                Node_array = np.asarray(node_array)
            else:
                for i in range(NumNodes / 10 + 1):
                    g = fzone.readline()
                    node_array = g.split()
                    # Convert string to integer array
                    node_array = [int(id) for id in node_array]
                    if (NumNodes - 10 * i < 10):
                        for j in range(NumNodes % 10):
                            Node_array[i * 10 + j] = node_array[j]
                    else:
                        for j in range(10):
                            Node_array[i * 10 + j] = node_array[j]
            fzone.close()
            print('--> Finished with zone file')

            Boundary_cell_area = np.zeros(NumNodes, 'float')
            for i in range(NumNodes):
                Boundary_cell_area[i] = 1.e20  # Fix the area to a large number

            print('--> Finished calculating boundary connections')

            boundary_cell_coord = [Cell_coord[Cell_id[i - 1] - 1] for i in Node_array]
            epsilon = 1e-0  # Make distance really small
            if (face == 'top'):
                boundary_cell_coord = [[cell[0], cell[1], cell[2] + epsilon] for cell in boundary_cell_coord]
            elif (face == 'bottom'):
                boundary_cell_coord = [[cell[0], cell[1], cell[2] - epsilon] for cell in boundary_cell_coord]
            elif (face == 'north'):
                boundary_cell_coord = [[cell[0], cell[1] + epsilon, cell[2]] for cell in boundary_cell_coord]
            elif (face == 'south'):
                boundary_cell_coord = [[cell[0], cell[1] - epsilon, cell[2]] for cell in boundary_cell_coord]
            elif (face == 'east'):
                boundary_cell_coord = [[cell[0] + epsilon, cell[1], cell[2]] for cell in boundary_cell_coord]
            elif (face == 'west'):
                boundary_cell_coord = [[cell[0] - epsilon, cell[1], cell[2]] for cell in boundary_cell_coord]
            elif (face == 'none'):
                boundary_cell_coord = [[cell[0], cell[1], cell[2]] for cell in boundary_cell_coord]
            elif (face == 'well1'):
                boundary_cell_coord = [[cell[0], cell[1], cell[2]] for cell in boundary_cell_coord]
            elif (face == 'well2'):
                boundary_cell_coord = [[cell[0], cell[1], cell[2]] for cell in boundary_cell_coord]
            elif (face == 'well3'):
                boundary_cell_coord = [[cell[0], cell[1], cell[2]] for cell in boundary_cell_coord]
            elif (face == 'well12'):
                boundary_cell_coord = [[cell[0], cell[1], cell[2]] for cell in boundary_cell_coord]
            elif (face == 'well22'):
                boundary_cell_coord = [[cell[0], cell[1], cell[2]] for cell in boundary_cell_coord]
            elif (face == 'well32'):
                boundary_cell_coord = [[cell[0], cell[1], cell[2]] for cell in boundary_cell_coord]
            else:
                sys.exit('ERROR: unknown face. Select one of: top, bottom, east, west, north, south.')
            with open(ex_file, 'w') as f:
                f.write('CONNECTIONS\t%i\n' % Node_array.size)
                for idx, cell in enumerate(boundary_cell_coord):
                    f.write('%i\t%.6e\t%.6e\t%.6e\t%.6e\n' % (
                        Node_array[idx], cell[0], cell[1], cell[2], Boundary_cell_area[idx]))
            print('--> Finished writing ex file "' + ex_file + '" corresponding to the zone file: ' + zone_file+'\n')
    os.rename("pointset_well1.ex","well_location_001.ex")
    os.rename("pointset_well2.ex","well_location_002.ex")
    os.rename("pointset_well3.ex","well_location_003.ex")
    os.rename("pointset_well12.ex","well_location_0012.ex")
    os.rename("pointset_well22.ex","well_location_0022.ex")
    os.rename("pointset_well32.ex","well_location_0032.ex")
    print('--> Converting zone files to ex complete')    

def inp2gmv(self, inp_file=''):
    """ Convert inp file to gmv file, for general mesh viewer .
    
    Kwargs:
        inp_file (str): name of inp file
    """

    if inp_file:
        self._inp_file = inp_file
    else:
        inp_file = self._inp_file

    if inp_file == '':
        sys.exit('ERROR: inp file must be specified in inp2gmv!')

    gmv_file = inp_file[:-4] + '.gmv'

    with open('inp2gmv.lgi', 'w') as fid:
        fid.write('read / avs / ' + inp_file + ' / mo\n')
        fid.write('dump / gmv / ' + gmv_file + ' / mo\n')
        fid.write('finish \n\n')

    cmd = lagrit_path + ' <inp2gmv.lgi ' + '>lagrit_inp2gmv.txt'
    failure = os.system(cmd)
    if failure:
        sys.exit('ERROR: Failed to run LaGrit to get gmv from inp file!')
    print("--> Finished writing gmv format from avs format")


def write_perms_and_correct_volumes_areas(self, inp_file='', uge_file='', perm_file='', aper_file=''):
    """ Write permeability values to perm_file, write aperture values to aper_file, and correct volume areas in uge_file 
    """
    print("--> Writing Perms and Correct Volume Areas")
    if inp_file:
        self._inp_file = inp_file
    else:
        inp_file = self._inp_file
    
    if inp_file == '':
        sys.exit('ERROR: inp file must be specified!')

    if uge_file:
        self._uge_file = uge_file
    else:
        uge_file = self._uge_file

    if uge_file == '':
        sys.exit('ERROR: uge file must be specified!')

    if perm_file:
        self._perm_file = perm_file
    else:
        perm_file = self._perm_file

    if perm_file == '' and self._perm_cell_file == '':
        sys.exit('ERROR: perm file must be specified!')

    if aper_file:
        self._aper_file = aper_file
    else:
        aper_file = self._aper_file

    if aper_file == '' and self._aper_cell_file == '':
        sys.exit('ERROR: aperture file must be specified!')

    mat_file = 'materialid.dat'
    t = time()
    # Make input file for C UGE converter
    f = open("convert_uge_params.txt", "w")
    f.write("%s\n"%inp_file)
    f.write("%s\n"%mat_file)
    f.write("%s\n"%uge_file)
    f.write("%s"%(uge_file[:-4]+'_vol_area.uge\n'))
    if self._aper_cell_file:
            f.write("%s\n"%self._aper_cell_file)
            f.write("1\n")
    else:
            f.write("%s\n"%self._aper_file)
            f.write("-1\n")
    f.close()

    cmd = os.environ['correct_uge_PATH']+ 'correct_uge' + ' convert_uge_params.txt' 
    failure = os.system(cmd)
    if failure > 0:
            sys.exit('ERROR: UGE conversion failed\nExiting Program')
    elapsed = time() - t
    print '--> Time elapsed for UGE file conversion: %0.3f seconds\n'%elapsed

    # need number of nodes and mat ID file
    print('--> Writing HDF5 File')
    materialid = np.genfromtxt(mat_file, skip_header = 3).astype(int)
    materialid = -1 * materialid - 6
    NumIntNodes = len(materialid)

    if perm_file:
        filename = 'dfn_properties.h5'
        h5file = h5py.File(filename, mode='w')
        print('--> Beginning writing to HDF5 file')
        print('--> Allocating cell index array')
        iarray = np.zeros(NumIntNodes, '=i4')
        print('--> Writing cell indices')
        # add cell ids to file
        for i in range(NumIntNodes):
            iarray[i] = i + 1
        dataset_name = 'Cell Ids'
        h5dset = h5file.create_dataset(dataset_name, data=iarray)

        print ('--> Allocating permeability array')
        perm = np.zeros(NumIntNodes, '=f8')

        print('--> reading permeability data')
        print('--> Note: this script assumes isotropic permeability')
        perm_list = np.genfromtxt(perm_file,skip_header = 1)
        perm_list = np.delete(perm_list, np.s_[1:5], 1)

        matid_index = -1*materialid - 7
        for i in range(NumIntNodes):
            j = matid_index[i]
            if int(perm_list[j,0]) == materialid[i]:
                    perm[i] = perm_list[j, 1]
            else:
                    sys.exit('Indexing Error in Perm File')

        dataset_name = 'Permeability'
        h5dset = h5file.create_dataset(dataset_name, data=perm)

        h5file.close()
        print("--> Done writing permeability to h5 file")
        del perm_list

    if self._perm_cell_file:
        filename = 'dfn_properties.h5'
        h5file = h5py.File(filename, mode='w')

        print('--> Beginning writing to HDF5 file')
        print('--> Allocating cell index array')
        iarray = np.zeros(NumIntNodes, '=i4')
        print('--> Writing cell indices')
        # add cell ids to file
        for i in range(NumIntNodes):
            iarray[i] = i + 1
        dataset_name = 'Cell Ids'
        h5dset = h5file.create_dataset(dataset_name, data=iarray)
        print ('--> Allocating permeability array')
        perm = np.zeros(NumIntNodes, '=f8')
        print('--> reading permeability data')
        print('--> Note: this script assumes isotropic permeability')
        f = open(self._perm_cell_file, 'r')
        f.readline()
        perm_list = []
        while True:
            h = f.readline()
            h = h.split()
            if h == []:
                break
            h.pop(0)
            perm_list.append(h)

        perm_list = [float(perm[0]) for perm in perm_list]
        
        dataset_name = 'Permeability'
        h5dset = h5file.create_dataset(dataset_name, data=perm_list)
        f.close()

        h5file.close()
        print('--> Done writing permeability to h5 file')

# hpham: this is to run PFLOTRAN.
def pflotran2(self):
    ''' Run pflotran
    Copy PFLOTRAN run file into working directory and run with ncpus
    '''
    print os.path.abspath(os.getcwd())
    print self._local_dfnFlow_file
    print('hpham: This run pflotran2 where you have to manually copy dfn_explicit.in from the input folder to here.')
# hpham: comment these lines to manually copy the input file
#    try: 
#            shutil.copy(os.path.abspath(self._dfnFlow_file), os.path.abspath(os.getcwd()))
#    except:
#            print("-->ERROR copying PFLOTRAN input file")
#            exit()
    print("="*80)
    print("--> Running PFLOTRAN") 
    cmd = os.environ['PETSC_DIR']+'/'+os.environ['PETSC_ARCH']+'/bin/mpirun -np ' + str(self._ncpu) + \
          ' ' + os.environ['PFLOTRAN_DIR']+'/src/pflotran/pflotran -pflotranin ' + self._local_dfnFlow_file 
    os.system(cmd)    
    print('='*80)
    print("--> Running PFLOTRAN Complete")
    print('='*80)
    print("\n")

 
def pflotran(self):
    ''' Run pflotran
    Copy PFLOTRAN run file into working directory and run with ncpus
    '''
    print os.path.abspath(os.getcwd())
#    print('hpham: This run pflotran2 where you have to manually copy dfn_explicit.in from the input folder to here.')
    try: 
            shutil.copy(os.path.abspath(self._dfnFlow_file), os.path.abspath(os.getcwd()))
            shutil.copy('/projects/DFN/apps/dfnWorks-Version2.0/pydfnworks/bin/gen_head.py', os.path.abspath(os.getcwd()))
    except:
            print("-->ERROR copying PFLOTRAN input file or gen_head.py")
            exit()
    print("="*80)
    print("--> Running PFLOTRAN") 
    cmd = os.environ['PETSC_DIR']+'/'+os.environ['PETSC_ARCH']+'/bin/mpirun -np ' + str(self._ncpu) + \
          ' ' + os.environ['PFLOTRAN_DIR']+'/src/pflotran/pflotran -pflotranin ' + self._local_dfnFlow_file 
    os.system(cmd)    
    print('='*80)
    print("--> Running PFLOTRAN Complete")
    print('='*80)
    print("\n") 

def get_well_loc_old(self):      #hpham
    print("hpham: This run get_well_loc_old in flow.py")
    print os.path.abspath(os.getcwd())   #hpham
    
    # [1] Opening uge file
    fuge = open('full_mesh.uge', 'r')
    # Reading cell ids, cells centers and cell volumes
    line = fuge.readline()
    line = line.split()
    NumCells = int(line[1])
    Cell_id = np.zeros(NumCells, 'int')
    Cell_coord = np.zeros((NumCells, 3), 'float')
    Cell_vol = np.zeros(NumCells, 'float')
    for cells in range(NumCells):
        line = fuge.readline()
        line = line.split()
        Cell_id[cells] = int(line.pop(0))
        line = [float(id) for id in line]
        Cell_vol[cells] = line.pop(3)
        Cell_coord[cells] = line
    fuge.close()

    # [2] Call LaGriT to get 'pointset.zone'
    for ifile in range(3):
        if ifile == 0:
            os.system('cp -f $model_path/obsloc_001.dat obsloc.dat')
            os.system('$lagrit_dfn < obsloc.dat > out_LaGrit.dat')
            os.system('sed -i -- ''s/000001/000007/g'' pointset.zone')     # find and replace
            os.system('cp -f pointset.zone obs1.zone')
        elif ifile == 1:
            os.system('cp -f $model_path/obsloc_002.dat obsloc.dat')
            os.system('$lagrit_dfn < obsloc.dat > out_LaGrit.dat')
            os.system('sed -i -- ''s/000001/000008/g'' pointset.zone')     # find and replace
            os.system('cp -f pointset.zone obs2.zone')            
        else:
            os.system('cp -f $model_path/obsloc_003.dat obsloc.dat')
            os.system('$lagrit_dfn < obsloc.dat > out_LaGrit.dat')
            os.system('sed -i -- ''s/000001/000009/g'' pointset.zone')     # find and replace
            os.system('cp -f pointset.zone obs3.zone')
        inp_file = 'pointset.zone'
        f = open(inp_file, 'r')
        f.readline()
        f.readline()
        f.readline()
        line = f.readline()
        NumNodes = int(line.strip(' ').split()[0])
        #    print NumNodes
        Node_array = np.zeros(NumNodes, 'int')
        node_id = []
        nrows = np.ceil((NumNodes / 10.0))
        #    print 'Number of node rows = ', nrows
        for i in range(int(nrows)):
            #       print i
            line = f.readline()
            node_array = line.split()
            # Convert string to integer array:
            node_array = [int(id) for id in node_array]
            if (NumNodes - 10 * i < 10):
                for j in range(NumNodes % 10):
                    Node_array[i * 10 + j] = node_array[j]
            else:
                for j in range(10):
                    Node_array[i * 10 + j] = node_array[j]
        f.close()

        # [3] Writing to values again to files
        Boundary_cell_area = np.zeros(NumNodes, 'float')
        for i in range(NumNodes):
            Boundary_cell_area[i] = 1.e20  # Fix the area to a large number
        epsilon = 0  # hpham: Doesn't need this ... Make distance really small
        boundary_cell_coord = [Cell_coord[Cell_id[i - 1] - 1] for i in Node_array]
        boundary_cell_coord = [[cell[0], cell[1], cell[2] + epsilon] for cell in boundary_cell_coord]
        print 'NumNodes=', NumNodes
        if ifile == 0:
            fid = open('well_location_001.ex', 'w')
        elif ifile == 1:
            fid = open('well_location_002.ex', 'w')
        else:
            fid = open('well_location_003.ex', 'w')
        fid.write('CONNECTIONS\t%i\n' % Node_array.size)
        for idx, cell in enumerate(boundary_cell_coord):
            fid.write('%i\t%.6e\t%.6e\t%.6e\t%.6e\n' % (Node_array[idx], cell[0], cell[1], cell[2], Boundary_cell_area[idx]))
        fid.close()

    # copies well zone files to allboundaries.zone
    fall=open('allboundaries.zone','a')     # open file and write at the end of the file. 
    #copy all but last 2 lines of obs1.zone in allboundaries.zone
    fzone=open('obs1.zone','rb')
    lines=fzone.readlines()
    lines=lines[1:-2]
    fzone.close() 
    fall.writelines(lines)

    #copy all but last 2 lines of obs2.zone in allboundaries.zone
    fzone=open('obs2.zone','rb')
    lines=fzone.readlines()
    lines=lines[1:-2]
    fzone.close() 
    fall.writelines(lines)

    #copy all lines of obs3.zone in allboundaries.zone
    fzone=open('obs3.zone','rb')
    lines=fzone.readlines()
#    lines=lines[1:-2]
    lines=lines[1:]
    fzone.close() 
    fall.writelines(lines)
    fall.close()



def get_well_loc(self):      #hpham
    print("hpham: This runs get_well_loc funtion (new version) in flow.py ... ")
    print os.path.abspath(os.getcwd())   #hpham
    
   #[0] Read well locations
    os.system('ln -s $model_path/well_locations.dat .')
    fwell_id = open('well_locations.dat', 'r')
    line = fwell_id.readline()
    line = line.split()
    nwells = int(line[0])
    print "Number of well = ", nwells
    x_wel = np.zeros(nwells, 'float')
    y_wel = np.zeros(nwells, 'float')
    print "Well location: X Y"
    for i in range(nwells):
    #    print i
        line = fwell_id.readline()

        line = line.split()
    #    print line
        x_wel[i] = line.pop(0)
        y_wel[i] = line.pop(0)
        print x_wel[i], y_wel[i]
    print "Done reading well_locations.dat"


    # [...] Opening uge file
    fuge = open('full_mesh_vol_area.uge', 'r')
    #Reading cell ids, cells centers and cell volumes
        
    line = fuge.readline()  # Read the first line
    line = line.split()
    NumCells = int(line[1])  # get number of Delaunay tri nodes or vor cells
    Cell_id = np.zeros(NumCells, 'int')
    Cell_coord = np.zeros((NumCells, 3), 'float')
    Cell_vol = np.zeros(NumCells, 'float')
    for cells in range(NumCells):
        line = fuge.readline()
        line = line.split()
        Cell_id[cells] = int(line.pop(0))
        line = [float(id) for id in line]
        Cell_vol[cells] = line.pop(3)
        Cell_coord[cells] = line
    
    line = fuge.readline()  # Read the line with CONNECTIONS
    line = line.split()
    NumConns = int(line[1])  # get number of connections
    Node1 = np.zeros(NumConns, 'int')
    Node2 = np.zeros(NumConns, 'int')
    for cells in range(NumConns):
        line = fuge.readline()
        line = line.split()
        Node1[cells] = int(line.pop(0))
        Node2[cells] = int(line.pop(0))
    fuge.close()


    # Read materialid.dat
    mid_file = "materialid.dat"
    print('Reading materialid.dat')
    cell_mat = np.zeros((NumCells,1),'int')   # Cell material
    f = open(mid_file,'r')
    for i in range(3):
        g = f.readline()
#        print g
    for i in range(NumCells):
        g = f.readline()
        g = g.split()
        cell_mat[i] = float(g.pop(0))


    # # Opening avs file
    # avs_file = "full_mesh.inp"
    # print('--> Opening avs file')
    # f = open(avs_file,'r')
    # g = f.readline()
    # g = g.split()

    # Num_Nodes = int(g.pop(0))
    # Num_Elem = int(g.pop(0))

    # print "Num_Nodes=", Num_Nodes
    # print "Num_Elem=", Num_Elem

    # node_loc = np.zeros((Num_Nodes,3),'float')
    # #node_id_all = np.zeros((Num_Nodes,1),'int')
    # for i in range(Num_Nodes):
    #     g = f.readline()
    #     g = g.split()
    # #    node_id_all[i] = int(g.pop(0))
    #     node_loc[i,0] = float(g.pop(1))
    #     node_loc[i,1] = float(g.pop(1))
    #     node_loc[i,2] = float(g.pop(1))
    # #    print node_loc[i,0], node_loc[i,1], node_loc[i,2]

    # # elem_type
    # # tet = 4
    # # brick = 8

    # #Elements = np.zeros((Num_Elem,9),'int')
    # for i in range(Num_Elem):
    #     g = f.readline()
    # #    g = g.split()
    # #    elem_type = g.pop(2)
    # #    if (elem_type == 'tet'):
    # #        Elements[i,0] = 4
    # #    for j in range(Elements[i,0]):
    # #        Elements[i,j+1] = g.pop(2)

    # # Read and skip 7 lines
    # for i in range(7):
    #     g = f.readline()
    # #    print g


    # Vertices = np.zeros((Num_Nodes,2),'float')
    # for i in range(Num_Nodes):
    #     g = f.readline()
    #     g = g.split()
    #     Vertices[i,0] = float(g.pop(0))
    #     Vertices[i,1] = float(g.pop(0))
    # #    Vertices[i,2] = float(g.pop(1))
    # #    Vertices[i,3] = float(g.pop(1))
    # #    Vertices[i,4] = float(g.pop(1))
    # #    Vertices[i,5] = float(g.pop(1))
    # #    Vertices[i,6] = float(g.pop(1))
    # #    print i, Vertices[i,0], Vertices[i,1]
    # print('--> Done with avs file')

	
	

    # [2] Call LaGriT to get 'pointset.zone'
    for ifile in range(3):
        if ifile == 0:
            os.system('cp -f $model_path/obsloc_001.dat .')
            os.system('mv obsloc_001.dat obsloc.dat')
            os.system('$lagrit_dfn < obsloc.dat > out_LaGrit.dat')
            os.system('sed -i -- ''s/000001/000007/g'' pointset.zone')     # find and replace
            os.system('cp -f pointset.zone obs1.zone')
        elif ifile == 1:
            os.system('cp -f $model_path/obsloc_002.dat .')
            os.system('mv obsloc_002.dat obsloc.dat')
            os.system('$lagrit_dfn < obsloc.dat > out_LaGrit.dat')
            os.system('sed -i -- ''s/000001/000008/g'' pointset.zone')     # find and replace
            os.system('cp -f pointset.zone obs2.zone')            
        else:
            os.system('cp -f $model_path/obsloc_003.dat .')
            os.system('mv obsloc_003.dat obsloc.dat')
            os.system('$lagrit_dfn < obsloc.dat > out_LaGrit.dat')
            os.system('sed -i -- ''s/000001/000009/g'' pointset.zone')     # find and replace
            os.system('cp -f pointset.zone obs3.zone')


        inp_file = 'pointset.zone'
        f = open(inp_file, 'r')
        f.readline()
        f.readline()
        f.readline()
        line = f.readline()
        NumNodes = int(line.strip(' ').split()[0])
#        print "NumNodes before =", NumNodes
        Node_array = np.zeros(NumNodes, 'int')
        #node_id = []
        nrows = np.ceil((NumNodes / 10.0))
        #    print 'Number of node rows = ', nrows
        for i in range(int(nrows)):
            #       print i
            line = f.readline()
            node_array = line.split()
            # Convert string to integer array:
            node_array = [int(id) for id in node_array]
            if (NumNodes - 10 * i < 10):
                for j in range(NumNodes % 10):
                    Node_array[i * 10 + j] = node_array[j]
            else:
                for j in range(10):
                    Node_array[i * 10 + j] = node_array[j]
        f.close()

        #print Node_array
        id = np.zeros(NumNodes, 'int')
        for i in range(Node_array.size):
            id[i] = np.where(Cell_id == Node_array[i])[0]
        #    print id


        # [3] Writing to values again to files
        #Boundary_cell_area = np.zeros(NumNodes, 'float')
        epsilon = 0  # hpham: Make distance really small
        boundary_cell_coord = [Cell_coord[Cell_id[i - 1] - 1] for i in Node_array]
        boundary_cell_coord = [[cell[0], cell[1], cell[2] + epsilon] for cell in boundary_cell_coord]
        


        # with open(ex_file, 'w') as f:
        #     f.write('CONNECTIONS\t%i\n' % Node_array.size)
        #     for idx, cell in enumerate(boundary_cell_coord):
        #         f.write('%i\t%.6e\t%.6e\t%.6e\t%.6e\n' % (
        #             Node_array[idx], cell[0], cell[1], cell[2], Boundary_cell_area[idx]))
        # print('--> Finished writing ex file "' + ex_file + '" corresponding to the zone file: ' + zone_file+'\n')




        data = np.zeros((NumNodes, 6), 'float')
#        for i in range(NumNodes):
        for i, cell in enumerate(boundary_cell_coord):
        #    print i, Node_array[i]
            data[i, 0] = Node_array[i]
            data[i, 1] = cell[0]
            data[i, 2] = cell[1]
            data[i, 3] = cell[2]
            data[i, 4] = 1.e20
            data[i, 5] = cell_mat[id[i]]
        #    print i, data[i, 0], data[i, 5]

        #np.savetxt('data.dat',data)
        frac_id = np.unique(data[:, 5])
        #print "frac_id=", frac_id
        id_final = np.zeros(frac_id.size, 'int')
#        id_final = []
        cnt = 0
        for i in range(frac_id.size):
#            print "FracID = ", frac_id[i]
            id2 = np.where(data[:, 5] == frac_id[i])[0]

        #    print "id2=", id2
            # Calculate the distance from each node to well
            min_dist = 999
            id_min = -999

            for j in range(id2.size):
                x_curr = data[int(id2[j]), 1]
                y_curr = data[int(id2[j]), 2]
                nodeid_curr = data[int(id2[j]), 0]
                id3 = np.where(Node1 == nodeid_curr)[0]
                NumConns_i = id3.size
#                print NumConns_i
        #        print x_curr, y_curr, x_wel[k], y_wel[k], (x_curr-x_wel[k])**2, (y_curr-y_wel[k])**2
                dist2well = ((x_curr-x_wel[ifile])**2+(y_curr-y_wel[ifile])**2)**(0.5)
#                print "ifile | node_id | dist2well | NumConns_i | ", ifile, nodeid_curr, dist2well, NumConns_i
                if (dist2well <= 999) and (dist2well < min_dist) and (NumConns_i >= 1):      # Only get nodes that have >= xxx connections
                    min_dist = dist2well
                    id_min = int(id2[j])
            if id_min != -999:
                id_final[cnt] = int(id_min)
                cnt=cnt+1    # Only count the fractures that satisfy all conditions above.
        #print id_final
        id4 = np.where(id_final > 0)[0]
        data_final = np.zeros((id4.size, 5), 'int')
        data_final = data[id_final[id4], :]
        #print data_final[:, 0].size
#        np.savetxt('out.dat',data_final)


        if ifile == 0:
            fid = open('well_location_001.ex', 'w')
            fid_zone = open('obs1.zone', 'w')
            line2_ = '000007     obs1'
        elif ifile == 1:
            fid = open('well_location_002.ex', 'w')
            fid_zone = open('obs2.zone', 'w')
            line2_ = '000008     obs2'            
        else:
            fid = open('well_location_003.ex', 'w')
            fid_zone = open('obs3.zone', 'w')
            line2_ = '000009     obs3'            
        NumNodes_final = data_final[:, 0].size    
        fid.write('CONNECTIONS\t%i\n' % NumNodes_final)
        for idx in range(data_final[:, 0].size):
            fid.write('%i\t%.6e\t%.6e\t%.6e\t%.6e\n' % (data_final[idx,0], data_final[idx,1], data_final[idx,2], data_final[idx,3], data_final[idx,4]))
        fid.close()

        # Re-write obs*.zone
        fid_zone.write('zone \n')
        fid_zone.write('%s \n' % line2_)
        fid_zone.write('nnum \n')
        fid_zone.write('%d \n' % NumNodes_final)

        nrows = np.ceil((NumNodes_final / 10.0))
        for i in range(int(nrows)):
            if (NumNodes_final - 10 * i < 10):
                for j in range(NumNodes_final % 10):
                    fid_zone.write('%11d' % data_final[i * 10 + j, 0])
            else:
                for j in range(10):
                    fid_zone.write('%11d' % data_final[i * 10 + j, 0])
            fid_zone.write('\n')
        fid_zone.write('\n')
        fid_zone.write('stop')
        f.close()
        fid_zone.close()







    # copies well zone files to allboundaries.zone
    fall=open('allboundaries.zone','a')     # open file and write at the end of the file. 
    #copy all but last 2 lines of obs1.zone in allboundaries.zone
    fzone=open('obs1.zone','rb')
    lines=fzone.readlines()
    lines=lines[1:-2]
    fzone.close() 
    fall.writelines(lines)

    #copy all but last 2 lines of obs2.zone in allboundaries.zone
    fzone=open('obs2.zone','rb')
    lines=fzone.readlines()
    lines=lines[1:-2]
    fzone.close() 
    fall.writelines(lines)

    #copy all lines of obs3.zone in allboundaries.zone
    fzone=open('obs3.zone','rb')
    lines=fzone.readlines()
#    lines=lines[1:-2]
    lines=lines[1:]
    fzone.close() 
    fall.writelines(lines)
    fall.close()


def gen_new_boundary_files(self):
    input_files = {'pboundary_back_n.ex', 'pboundary_left_w.ex', 'pboundary_right_e.ex', 'pboundary_front_s.ex'}

    # Load datum data
    os.system('cp -f $model_path/datum.dat .')
    datum_tmp=np.loadtxt('datum.dat',delimiter=',', skiprows=1)
    print 'Running gen_new_boundary_files ... make sure you check file datum.dat for correctness!'

    for ifile in input_files:
        print ifile
        if ifile == 'pboundary_back_n.ex':
            datum = datum_tmp[0, :]
            dir = 0   # 0: x direction
            ofname1 = ifile[:-3] + '_right.ex2'
            ofname2 = ifile[:-3] + '_left.ex2'
        elif ifile == 'pboundary_front_s.ex':
            datum = datum_tmp[1, :]
            dir = 0
            ofname1 = ifile[:-3] + '_right.ex2'
            ofname2 = ifile[:-3] + '_left.ex2'
        elif ifile == 'pboundary_right_e.ex':
            datum = datum_tmp[2, :]
            dir = 1  # 1: y direction
            ofname1 = ifile[:-3] + '_top.ex2'
            ofname2 = ifile[:-3] + '_bot.ex2'
        else:    #pboundary_left_w.ex
            datum = datum_tmp[3, :]
            dir = 1
            ofname1 = ifile[:-3] + '_top.ex2'
            ofname2 = ifile[:-3] + '_bot.ex2'
        # Load data
        data=np.loadtxt(ifile, skiprows=1)


        # Split matrix by x or y or z
        #print datum
        data1 = data[data[:, dir+1] > datum[dir]]
        data2 = data[data[:, dir+1] <= datum[dir]]
        if os.path.exists(ofname1) == 'true':
            os.remove(ofname1)
        if os.path.exists(ofname2) == 'true':
            os.remove(ofname2)

        fid1 = open(ofname1, 'a')
        fid1.write('CONNECTIONS\t%i\n' % data1.shape[0])
        np.savetxt(fid1, data1, fmt='%i\t%.6e\t%.6e\t%.6e\t%.6e', newline='\n')
        fid1.close()


        fid2 = open(ofname2, 'a')
        fid2.write('CONNECTIONS\t%i\n' % data2.shape[0])
        np.savetxt(fid2, data2, fmt='%i\t%.6e\t%.6e\t%.6e\t%.6e', newline='\n')
        fid2.close()








def pflotran_cleanup(self):
    '''pflotran_cleanup
    Concatenate PFLOTRAN output files and then delete them 
    '''
    print '--> Processing PFLOTRAN output' 
    
    cmd = 'cat '+self._local_dfnFlow_file[:-3]+'-cellinfo-001-rank*.dat > cellinfo.dat'
    os.system(cmd)

    cmd = 'cat '+self._local_dfnFlow_file[:-3]+'-darcyvel-001-rank*.dat > darcyvel.dat'
    os.system(cmd)

    for fl in glob.glob(self._local_dfnFlow_file[:-3]+'-cellinfo*.dat'):
            os.remove(fl)    
    for fl in glob.glob(self._local_dfnFlow_file[:-3]+'-darcyvel*.dat'):
            os.remove(fl)    

def create_dfn_flow_links():
    os.symlink('../full_mesh.uge', 'full_mesh.uge')
    os.symlink('../full_mesh_vol_area.uge', 'full_mesh_vol_area.uge')
    os.symlink('../full_mesh.inp', 'full_mesh.inp')
    os.symlink('../pboundary_back_n.zone', 'pboundary_back_n.zone')
    os.symlink('../pboundary_front_s.zone', 'pboundary_front_s.zone')
    os.symlink('../pboundary_left_w.zone', 'pboundary_left_w.zone')
    os.symlink('../pboundary_right_e.zone', 'pboundary_right_e.zone')
    os.symlink('../pboundary_top.zone', 'pboundary_top.zone')
    os.symlink('../pboundary_bottom.zone', 'pboundary_bottom.zone')
    os.symlink('../materialid.dat', 'materialid.dat')
    
def uncorrelated(sigma):
    print '--> Creating Uncorrelated Transmissivity Fields'
    print 'Variance: ', sigma
    print 'Running un-correlated'
    x = np.genfromtxt('../aperture.dat', skip_header = 1)[:,-1]
    k = np.genfromtxt('../perm.dat', skip_header = 1)[0,-1]
    n = len(x)

    print np.mean(x)

    perm = np.log(k)*np.ones(n) 
    perturbation = np.random.normal(0.0, 1.0, n)
    perm = np.exp(perm + np.sqrt(sigma)*perturbation) 

    aper = np.sqrt((12.0*perm))
    aper -= np.mean(aper)
    aper += np.mean(x)

    print '\nPerm Stats'
    print '\tMean:', np.mean(perm)
    print '\tMean:', np.mean(np.log(perm))
    print '\tVariance:',np.var(np.log(perm))
    print '\tMinimum:',min(perm)
    print '\tMaximum:',max(perm)
    print '\tMinimum:',min(np.log(perm))
    print '\tMaximum:',max(np.log(perm))

    print '\nAperture Stats'
    print '\tMean:', np.mean(aper)
    print '\tVariance:',np.var(aper)
    print '\tMinimum:',min(aper)
    print '\tMaximum:',max(aper)


    output_filename = 'aperture_' + str(sigma) + '.dat'
    f = open(output_filename,'w+')
    f.write('aperture\n')
    for i in range(n):
        f.write('-%d 0 0 %0.5e\n'%(i + 7, aper[i]))
    f.close()

    cmd = 'ln -s ' + output_filename + ' aperture.dat '
    os.system(cmd)

    output_filename = 'perm_' + str(sigma) + '.dat'
    f = open(output_filename,'w+')
    f.write('permeability\n')
    for i in range(n):
        f.write('-%d 0 0 %0.5e %0.5e %0.5e\n'%(i+7, perm[i], perm[i], perm[i]))
    f.close()

    cmd = 'ln -s ' + output_filename + ' perm.dat '
    os.system(cmd) 
        

def parse_pflotran_vtk(self, grid_vtk_file=''): 
    """ Using C++ VTK library, convert inp file to VTK file, then change name of CELL_DATA to POINT_DATA.
    """
    print '--> Parsing PFLOTRAN output using C++'
    files = glob.glob('*-[0-9][0-9][0-9].vtk')
    out_dir = 'parsed_vtk'
    vtk_filename_list = []
    replacements = {'CELL_DATA':'POINT_DATA'} 
    header = ['# vtk DataFile Version 2.0\n',
              'PFLOTRAN output\n',
              'ASCII\n']
   
    inp_file = self._inp_file
    inp_file_copy = self._inp_file[:-4] + '_copy.inp'
    subprocess.call('cp ' + inp_file + ' ' + inp_file_copy, shell=True)
    jobname = self._jobname + '/'

    for fle in files:

        if os.stat(fle).st_size == 0:
            print 'ERROR: opening an empty pflotran output file'
            exit()
        
        temp_file = fle[:-4] + '_temp.vtk'
        with open(fle, 'r') as infile, open(temp_file, 'w') as outfile:
            ct = 0 
            for line in infile:
                if 'CELL_DATA' in line:
                    num_cells = line.strip(' ').split()[1]
                    outfile.write('POINT_DATA\t ' + num_cells + '\n')
                else: 
                    outfile.write(line)
        infile.close()
        outfile.close()
        vtk_filename = out_dir + '/' + fle.split('/')[-1]
        if not os.path.exists(os.path.dirname(vtk_filename)):
            os.makedirs(os.path.dirname(vtk_filename))
        arg_string = os.environ['VTK_PATH'] + ' '  +  jobname + inp_file + ' ' + jobname + vtk_filename  
        subprocess.call(arg_string, shell=True)
        arg_string = 'tail -n +6 ' + jobname + temp_file + ' > ' + jobname + temp_file + '.tmp && mv ' + jobname + temp_file +  '.tmp ' + jobname + temp_file  
        subprocess.call(arg_string, shell=True)
        arg_string = 'cat ' +  jobname + temp_file + ' >> ' + jobname + vtk_filename
        subprocess.call(arg_string, shell=True) 

    print '--> Parsing PFLOTRAN output complete'

def inp2vtk_python(self, inp_file=''):
    import pyvtk as pv
    """ Using Python VTK library, convert inp file to VTK file.  then change name of CELL_DATA to POINT_DATA.
    """
    print("--> Using Python to convert inp files to VTK files")
    if self._inp_file:
        inp_file = self._inp_file
    else:
        self._inp_file = inp_file

    if inp_file == '':
        sys.exit('ERROR: Please provide inp filename!')

    if self._vtk_file:
        vtk_file = self._vtk_file
    else:
        vtk_file = inp_file[:-4]
        self._vtk_file = vtk_file + '.vtk'

    print("--> Reading inp data")

    with open(inp_file, 'r') as f:
        line = f.readline()
        num_nodes = int(line.strip(' ').split()[0])
        num_elems = int(line.strip(' ').split()[1])

        coord = np.zeros((num_nodes, 3), 'float')
        elem_list_tri = []
        elem_list_tetra = []

        for i in range(num_nodes):
            line = f.readline()
            coord[i, 0] = float(line.strip(' ').split()[1])
            coord[i, 1] = float(line.strip(' ').split()[2])
            coord[i, 2] = float(line.strip(' ').split()[3])

        for i in range(num_elems):
            line = f.readline().strip(' ').split()
            line.pop(0)
            line.pop(0)
            elem_type = line.pop(0)
            if elem_type == 'tri':
                elem_list_tri.append([int(i) - 1 for i in line])
            if elem_type == 'tet':
                elem_list_tetra.append([int(i) - 1 for i in line])

    print('--> Writing inp data to vtk format')

    vtk = pv.VtkData(pv.UnstructuredGrid(coord, tetra=elem_list_tetra, triangle=elem_list_tri),
                     'Unstructured pflotran grid')
    vtk.tofile(vtk_file)


def parse_pflotran_vtk_python(self, grid_vtk_file=''):
    """ Replace CELL_DATA with POINT_DATA in the VTK output."""
    print '--> Parsing PFLOTRAN output with Python'
    if grid_vtk_file:
        self._vtk_file = grid_vtk_file
    else:
        self.inp2vtk_python()

    grid_file = self._vtk_file
    
    files = glob.glob('*-[0-9][0-9][0-9].vtk')
    with open(grid_file, 'r') as f:
        grid = f.readlines()[3:]

    out_dir = 'parsed_vtk'
    for line in grid:
        if 'POINTS' in line:
            num_cells = line.strip(' ').split()[1]

    for file in files:
        with open(file, 'r') as f:
            pflotran_out = f.readlines()[4:]
        pflotran_out = [w.replace('CELL_DATA', 'POINT_DATA ') for w in pflotran_out]
        header = ['# vtk DataFile Version 2.0\n',
                  'PFLOTRAN output\n',
                  'ASCII\n']
        filename = out_dir + '/' + file
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        with open(filename, 'w') as f:
            for line in header:
                f.write(line)
            for line in grid:
                f.write(line)
            f.write('\n')
            f.write('\n')
            if 'vel' in file:
                f.write('POINT_DATA\t ' + num_cells + '\n')
            for line in pflotran_out:
                f.write(line)
    print '--> Parsing PFLOTRAN output complete'
