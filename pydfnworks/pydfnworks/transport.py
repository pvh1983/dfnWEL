# Copyright 2019, Nevada System of Higher Education on Behalf of the Desert Research Institute, and
# Copyright 2018, Triad National Secuirty, LLC. All rights reserved.

# This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation; either version 3.0 of the License,
# or (at your option) any later version.  This software requires you separately obtain dfnWorks under an
# appropriate license from the Los Alamos National Laboratory (LANL), which is operated by
# Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.

import os
import sys
import shutil
import helper
from time import time


def dfn_trans(self):
    '''dfnTrans
    Copy input files for dfnTrans into working directory and run DFNTrans
    '''
    print('='*80)
    print("\ndfnTrans Starting\n")
    print('='*80)

    self.copy_dfn_trans_files()
    tic = time()
    self.run_dfn_trans()
    # self.cleanup_files_at_end()
    #helper.dump_time(self._jobname, 'Process: dfnTrans', time() - tic)


def copy_dfn_trans_files(self):
    '''create link to DFNTRANS and copy input file into local directory
    '''
#    import cylinder_ic ### hpham

    # Create Path to DFNTrans
    try:
        os.symlink(os.environ['DFNTRANS_PATH']+'DFNTrans', './DFNTrans')
    except OSError:
        os.remove('DFNTrans')
        os.symlink(os.environ['DFNTRANS_PATH']+'DFNTrans', './DFNTrans')
    except:
        sys.exit("Cannot create link to DFNTrans. Exiting Program")

    # Copy DFNTrans input file
    print(os.getcwd())

    print("Attempting to Copy %s\n" % self._dfnTrans_file)
    try:
        shutil.copy(self._dfnTrans_file, os.path.abspath(os.getcwd()))
        # shutil.copy('/projects/DFN/apps/dfnWorks-Version2.0/pydfnworks/bin/cylinder_ic.py', os.path.abspath(os.getcwd())) #hpham
        # shutil.copy('/projects/DFN/apps/dfnWorks-Version2.0/pydfnworks/bin/sphere_ic.py', os.path.abspath(os.getcwd())) #hpham
    except OSError:
        print("--> Problem copying %s file" % self._local_dfnTrans_file)
        print("--> Trying to delete and recopy")
        os.remove(self._local_dfnTrans_file)
        shutil.copy(self._dfnTrans_file, os.path.abspath(os.getcwd()))
    except:
        print("--> ERROR: Problem copying %s file" % self._dfnTrans_file)
        sys.exit("Unable to replace. Exiting Program")


def run_dfn_trans(self):
    '''run dfnTrans simulation'''
    print("hpham: run_dfn_trans in transport.py")
    failure = os.system('./DFNTrans '+self._local_dfnTrans_file)
    if failure == 0:
        print('='*80)
        print("\ndfnTrans Complete\n")
        print('='*80)
    else:
        sys.exit("--> ERROR: dfnTrans did not complete\n")


def run_dfn_trans_sphere_ic(self):
    '''run dfnTrans simulation'''
#    import RandomPositGener ### hpham
    import sphere_ic  # hpham
    failure = os.system('./DFNTrans '+self._local_dfnTrans_file)
    if failure == 0:
        print('='*80)
        print("\ndfnTrans Complete\n")
        print('='*80)
    else:
        sys.exit("--> ERROR: dfnTrans did not complete\n")


def run_dfn_trans_cylinder_ic(self):
    '''run dfnTrans simulation'''
#    import RandomPositGener ### hpham
#    import cylinder_ic ### hpham
    failure = os.system('./DFNTrans '+self._local_dfnTrans_file)
    if failure == 0:
        print('='*80)
        print("\ndfnTrans Complete\n")
        print('='*80)
    else:
        sys.exit("--> ERROR: dfnTrans did not complete\n")


def create_dfn_trans_links():
    os.symlink('../params.txt', 'params.txt')
    os.symlink('../allboundaries.zone', 'allboundaries.zone')
    os.symlink('../tri_fracture.stor', 'tri_fracture.stor')
    os.symlink('../poly_info.dat', 'poly_info.dat')

# Copyright 2019, Nevada System of Higher Education on Behalf of the Desert Research Institute, and
# Copyright 2018, Triad National Secuirty, LLC. All rights reserved.

# This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation; either version 3.0 of the License,
# or (at your option) any later version.  This software requires you separately obtain dfnWorks under an
# appropriate license from the Los Alamos National Laboratory (LANL), which is operated by
# Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
