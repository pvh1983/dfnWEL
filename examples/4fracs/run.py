#"""
#   :synopsis: run file for dfnworks 
#   :version: 1.0
#   :maintainer: Jeffrey Hyman, Carl Gable, Nathaniel Knapp
#.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
#"""

import os, sys
from time import time
from pydfnworks import * 
import subprocess
#import cylinder_ic # hpham
#import sphere_ic # hpham
#import gen_head # hpham



define_paths()
main_time = time()
DFN = dfnworks.create_dfn()

DFN.make_working_directory()
DFN.check_input()
DFN.create_network()       # no need if run flow only
#DFN.output_report()        # no need
DFN.mesh_network()

DFN.lagrit2pflotran()
#DFN.get_well_loc()       # hpham no need
#DFN.get_well_loc_old()   # hpham no need 

DFN.pflotran()
DFN.parse_pflotran_vtk_python()       
DFN.pflotran_cleanup()   
#gen_head.hhead(997.32,9.81,1094.8,100,1)   # hhead(rho,gra,ele_head,domain_thickness,write_output) 

DFN.copy_dfn_trans_files()
#sphere_ic.ic("full_mesh.vtk",50000,48,(0.0,0.0,0.0))
DFN.run_dfn_trans()
#DFN.run_dfn_trans_hpham()   # include sphere_ic.py

#DFN.run_dfn_trans()
#DFN.run_dfn_trans_sphere_ic()   # hpham: To run sphere_ic.py
#DFN.run_dfn_trans_cylinder_ic()   # hpham: To run cylinder_ic.py

#DFN.cleanup_files_at_end()

main_elapsed = time() - main_time
timing = 'Time Required: %0.2f Minutes'%(main_elapsed/60.0)
print("*"*80)
print(DFN._jobname+' complete')
print("Thank you for using dfnWorks")
print("*"*80)

'''
if __name__ == "__main__":

    define_paths()
    main_time = time()
    print 'Compiling executables'
    subprocess.call(os.environ['python_dfn'] + ' compile.py', shell=True)  
    
    DFN = create_dfn()
    if type(DFN) is ' NoneType':
        print 'ERROR: DFN object not created correctly'
        exit()
    # General Work Flow
    DFN.dfn_gen() # Comment out if this is the first run
    DFN.dfn_flow()
#    DFN.dfn_trans() # Org code
    DFN.copy_dfn_trans_files() # hpham: copy input files for DFNTrans
    DFN.run_dfn_trans() # New hpham
    DFN.cleanup_files_at_end()

    #hpham: Speficy a full path to tje WEL parameter file
    #os.system('export ifile_well="/home/hpham/apps/dfnWorks-Version2.0/tests/4frac.wel"')
    os.environ['ifile_well'] = '/home/hpham/apps/dfnWorks-Version2.0/tests/4frac.wel'
    #os.environ['PETSC_DIR']='/home/hpham/apps/petsc'
    p = os.environ['ifile_well']
    print "Path to ifile_well:", p
    
    #sys.exit()


    main_elapsed = time() - main_time
    timing = 'Time Required: %0.2f Minutes'%(main_elapsed/60.0)
    print timing
    dump_time(DFN._local_jobname, DFN._jobname,main_elapsed) 
    #dfn.print_run_time()	
    print("*"*80)
    print(DFN._jobname+' complete')
    print("Thank you for using dfnWorks")
    print("*"*80)

'''