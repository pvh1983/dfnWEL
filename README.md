dfnWEL - A Method to Represent a Well in a Three-Dimensional Discrete Fracture Network Model for dfnWorks 2.0 (https://github.com/lanl/dfnWorks)

To use dfnWEL:

1. Download the dfnWorks 2.0 from the LANL GitHub at https://github.com/lanl/dfnWorks 
2. Download all dfnWEL files and overwrite the origin dfnWorks 2.0 files
3. Specifile the path to:
	(1) LAGRIT for InitialParPositions.c, line 452 (e.g., lagritpath = '/home/hpham/apps/LAGRIT/lagrit_ulin3.2')
	(2) the well paramter files *.wel at line 698, file /pydfnworks/pydfnworks/lagrit_scripts.py (will generalize this later).
4. Re-compile dfnWorks 2.0 code (go to folder pydfnworks and run python setup.py install).
5. Re-compile DFNTrans if needed (go to folder ParticleTracking and type make from your terminal. Make sure you have a suitable compiler. For me, I am using module load gnu/4.9.4). 

To run a test example:
1. Go to examples/4fracs folder and run the script run.sh
2. Change the paths for input files
3. Specify the path to the DFN mesh input file inside 4fracs/gen_4_user_rectangles.dat, at line 424


Other notes: 

1. This dfnWEL use Python 2.x
2. Search keywords hpham to see new modifications both in the input files and the source codes. 



