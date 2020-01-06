# Copyright 2019, Nevada System of Higher Education on Behalf of the Desert Research Institute, and
# Copyright 2018, Triad National Secuirty, LLC. All rights reserved.

# This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation; either version 3.0 of the License,
# or (at your option) any later version.  This software requires you separately obtain dfnWorks under an
# appropriate license from the Los Alamos National Laboratory (LANL), which is operated by
# Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.

from tempfile import mkstemp
from shutil import move
import os
import subprocess


def valid(name):
    if not (os.path.isfile(os.path.abspath(os.environ[name])) or os.path.isdir(os.path.abspath(os.environ[name]))):
        error_msg = "ERROR: " + name + \
            " has an invalid path name: " + os.environ[name]
        print error_msg
        exit()


def define_paths():

    # ================================================
    # THESE PATHS MUST BE SET BY THE USER.
    # ================================================

    # the dfnWorks-Version2.0  repository
    #    os.environ['DFNWORKS_PATH'] = '/home/hpham/dfnWorks-Version2.0/'
    os.environ['DFNWORKS_PATH'] = '/home/hpham/apps/dfnWorks-Version2.0/'
#    os.environ['DFNWORKS_PATH'] = '/projects/DFN/apps/orgdfnWorks/dfnWorks-Version2.0/'
    valid('DFNWORKS_PATH')
    if not (os.path.isdir(os.path.abspath(os.environ['DFNWORKS_PATH'] + 'tests/'))):
        print "INVALID VERSION OF DFNWORKS - does not have tests folder of official release 2.0"
        exit()

    # PETSC paths
    os.environ['PETSC_DIR'] = '/home/hpham/apps/petsc'  # Updated 07/05/18
#    os.environ['PETSC_DIR']='/home/hpham/apps/petsc'
    os.environ['PETSC_ARCH'] = 'arch-linux2-c-debug'
    valid('PETSC_DIR')
#    valid('PETSC_ARCH')

    # PFLOTRAN path
    # Updated 07/05/18
    os.environ['PFLOTRAN_DIR'] = '/home/hpham/apps/pflotran'
#    os.environ['PFLOTRAN_DIR']='/projects/DFN/apps/pflotran-dev'
#    os.environ['PFLOTRAN_DIR']='/projects/DFN/apps/new_pflotran'
    valid('PFLOTRAN_DIR')

    # Python executable
    os.environ['python_dfn'] = '/home/hpham/anaconda3/envs/b2dfn/bin/python'
    valid('python_dfn')

    # LaGriT executable
#    os.environ['lagrit_dfn'] = '/n/swqa/LAGRIT/lagrit.lanl.gov/downloads/lagrit_ulin3.2'
    os.environ['lagrit_dfn'] = '/home/hpham/apps/LAGRIT/lagrit_ulin3.2'
    valid('lagrit_dfn')

    # ===================================================
    # THESE PATHS ARE AUTOMATICALLY SET. DO NOT CHANGE.
    # ====================================================

    # Directories
    os.environ['DFNGEN_PATH'] = os.environ['DFNWORKS_PATH']+'DFNGen/'
    os.environ['DFNTRANS_PATH'] = os.environ['DFNWORKS_PATH'] + \
        'ParticleTracking/'
    os.environ['PYDFNWORKS_PATH'] = os.environ['DFNWORKS_PATH'] + 'pydfnworks/'
    os.environ['connect_test'] = os.environ['DFNWORKS_PATH'] + \
        'DFN_Mesh_Connectivity_Test/'
    os.environ['correct_uge_PATH'] = os.environ['DFNWORKS_PATH'] + \
        'C_uge_correct/'
    os.environ['VTK_PATH'] = os.environ['DFNWORKS_PATH'] + 'inp_2_vtk/'

# Copyright 2019, Nevada System of Higher Education on Behalf of the Desert Research Institute, and
# Copyright 2018, Triad National Secuirty, LLC. All rights reserved.

# This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation; either version 3.0 of the License,
# or (at your option) any later version.  This software requires you separately obtain dfnWorks under an
# appropriate license from the Los Alamos National Laboratory (LANL), which is operated by
# Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
