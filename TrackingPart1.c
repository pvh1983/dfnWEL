#include <stdio.h>
#include <search.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "FuncDef.h" 
#include <unistd.h>
#include <time.h>
#include<sys/stat.h>
#include <stdio.h>

struct inpfile {
  char filename[120];
  long int flag;
  double param;
}; 

struct lagrangian{
  double tau;
  double betta;
  double initx;
  double inity;
  double initz;

};

struct lagrangian lagvariable;
unsigned int FLAG_OUT=0; 
unsigned int np, t, nodeID=0, avs_o=0, traj_o=0, curv_o=0, no_out=0;
unsigned int marfa=0, plumec=0, disp_o=0, timecounter=0, frac_o=0;
struct intcoef {
  double weights[3];
};

struct posit3d {
  double cord3[3];
};

struct tempout{  // strucuture of the particle data, saved temporary for output perposes
  unsigned int times; // count of the time steps
  double position2d[2]; // 2d position of particles on a fracture plane
  double position3d[3]; // 3d position of particle in the 3d space
  double velocity3d[3]; //particles velocity, 3d
  unsigned int cellp; // number of cell in current particle position
  unsigned int fracturep; // number of fracture in current particle position
  double timep; // particle's current travel time
  double betap; // current beta parameter
  double length_t; //particle's trajectory length
  double pressure; // fluid pressure at particle's position
};

struct tempout *tempdata;
    

static FILE *wpt;
static FILE *p2well2;   // hpham
static FILE *fp_hpham;   // hpham
static FILE *wpt_att;
static FILE *wv;
static FILE *wint;
double x0w, y0w, rwell,top_well, bot_well, opt_mwell_on;  //hpham
int counthp1=0;    


//////////////////////////////////////////////////////////////////////////////
void ParticleTrack ()
{

  /*********Functions holds the loops on particles and loop on time steps for 
	    every partcile. Controls data output. ******************************/
    
  /*** read output options *****/    
  int res, par_ts_freq, par_np_freq, par_prt_cutoff;
  struct inpfile inputfile;
  char filename[125];
  unsigned int tort_o=0;
  int hp_count=1 ;     //hpham_var
  double rpw_sum=0, rpw=1000, rpw0=0;
  // Output of tortuosity file 
  inputfile=Control_File_Optional("out_tort:",9);
  if (inputfile.flag<0)
    tort_o=0;
  else
    {
      res=strncmp(inputfile.filename,"yes",3);
      if (res==0)
	tort_o=1;
    }

   // Output of fractures ID 
  inputfile=Control_File_Optional("out_fract:",10);
  if (inputfile.flag<0)
    frac_o=0;
  else
    {
      res=strncmp(inputfile.filename,"yes",3);
      if (res==0)
        frac_o=1;
    }

  // Output according to trajectory curvature (not every time step)  
  inputfile=Control_File("out_curv:",9);
  res=strncmp(inputfile.filename,"yes",3);
  if (res==0)
    curv_o=1;
   
  // output AVS file for each trajectory       
  inputfile=Control_File("out_avs:",8);
  res=strncmp(inputfile.filename,"yes",3);
  if (res==0)
    avs_o=1;

  // Output ASCII trajectory file    
  inputfile=Control_File("out_traj:",9);
  res=strncmp(inputfile.filename,"yes",3);
  if (res==0)
    traj_o=1;

  // Output MARFA input file
  inputfile=Control_File_Optional("out_marfa:",10);
  if (inputfile.flag<0)
    marfa=0;
  else
    {
      res=strncmp(inputfile.filename,"yes",3);
      if (res==0)
	{
	  marfa=1;
	  traj_o=1;
	  curv_o=1;
	} 
    }

  // Output PLUMECALC input file
  inputfile=Control_File_Optional("out_plumecalc:",14);
  if (inputfile.flag<0)
    plumec=0;
  else 
    {
      res=strncmp(inputfile.filename,"yes",3);
      if (res==0)
	{
	  plumec=1;
	  traj_o=1;
	  curv_o=1;
	}
    }

  // Output Time control planes
  inputfile=Control_File_Optional("out_disp:",9);
  if (inputfile.flag<0)
    disp_o=0;
  else
    {
      res=strncmp(inputfile.filename,"yes",3);
      if (res==0 )
	disp_o=1;
    }
      
  //if there is an additional output - hpham commented these two lines. 
  if ((avs_o+traj_o)==0)
    no_out=1;   
  
  // open file with total results         
  inputfile = Control_File("out_time:",9 );
  sprintf(filename,"%s/%s",maindir,inputfile.filename);
  FILE *tp = OpenFile (filename,"w");
  fprintf(tp,"# of time steps, flux weights, total travel time, x-, y-, z-final pos, beta, total length[m] \n");

  // hpham read input parameters for a monitoring well
  inputfile=Control_Param("mwell_on:",8);
  opt_mwell_on=inputfile.param;
//  printf("opt_mwell_on = %f \n",opt_mwell_on);
  inputfile=Control_Param("wloc_x:",7);
  x0w=inputfile.param;
  inputfile=Control_Param("wloc_y:",7);
  y0w=inputfile.param;
  inputfile=Control_Param("wrad:",5);
  rwell=inputfile.param;  
  inputfile=Control_Param("wtop:",5);
  top_well=inputfile.param;  
  inputfile=Control_Param("wbot:",5);
  bot_well=inputfile.param; 
  // Read input frequency to print out particle's info
  inputfile=Control_Param("out_par_freq:",13);
  par_np_freq=inputfile.param;     // print every par_np_freq particles
  inputfile=Control_Param("out_ts_freq:",12);
  par_ts_freq=inputfile.param;     // print out every par_ts_freq timesteps
  inputfile=Control_Param("out_prt_cutoff:",15);
  par_prt_cutoff=inputfile.param;     // print out every par_ts_freq timesteps
  
  
  
  /* initial positions of particle, output file****/
  FILE *inp;
  int outinit=0;
  inputfile=Control_File("out_init:",9);
  res=strncmp(inputfile.filename,"yes",3);
  
  if (res==0)
    {
      sprintf(filename,"%s/initpos.txt",maindir);      //hpham added .txt here
      inp = OpenFile (filename,"w");
      outinit=1;
      fprintf(inp,"NoofParticle, NoofCell, NoofFracture, x-pos, y-pos, z-pos, fluxweight \n"); 
    }

  // open tortuosity file
  char path[125];
  FILE *tort;
  if (tort_o>0)
    {
      sprintf(filename,"%s/torts.dat",maindir);
      tort =  OpenFile (filename,"w");
      fprintf(tort,"Data for trortuosity calculation: total length of trajectory, x-,y-, z- of initial pos., x-, y-, z- final pos., number of intersections. \n");
    } 
 
 // open FractureID file
   FILE *frac;
    
   if (frac_o>0)
     {
      sprintf(filename,"%s/FractureID",maindir);
      frac =  OpenFile (filename,"w");
      fprintf(frac,"Fractures ID  \n");

     }
  // Create path for trajectory outputs
  if (no_out!=1)
    {
      //out_path:
 
      inputfile = Control_File("out_path:",9 );
      sprintf(path,"%s/%s",maindir,inputfile.filename );
   
  
      mkdir(path, 0777);
      printf("\n All output trajectory files will be written in %s/ \n", path);
      printf(" Note: if the directory exists all the files in it will be replaced. \n \n");

          
    }
    
    
  /**** settings for Control Plane/Cylinder Output *****/
    
  int out_control=0, out_plane=0, out_cylinder=0; 
    
  inputfile=Control_File("ControlPlane:",13);
  res=String_Compare(inputfile.filename,"yes");
  if (res==0)
    out_plane=1;
  else
    {  
      inputfile=Control_File_Optional("ControlCylinder:",16);
      if (inputfile.flag<0)
	out_cylinder=0;
      else
	{
	  res=strncmp(inputfile.filename,"yes",3);
	  if (res==0)
	    out_cylinder=1;
	}
    }
 
       
  if ((out_cylinder+out_plane)>0)
    out_control=1;
 
  // reading variables for dispersivity (time-control) calculation

  unsigned int time_d=1, kd=1;
  double dtime=0.0;
  double epsl=0.0;
  if (disp_o==1)
    {
      inputfile = Control_Param("out_dtimest:",12 );
      time_d=(int)inputfile.param;

   
      inputfile = Control_Param("out_dtime:",10 );
      dtime=inputfile.param;
    
      epsl=0.05*dtime;
    }

  double sumsquares[time_d][2];//keeps sum of squares, like (x0-x)^2 for all particles 
  FILE * dis;
  unsigned int nd[time_d];
   
    
  int itime;
  char filetime[15];
  double inflowcoord=0.0; //define automatically
  double outflowcoord=0.0;
  double controllength=0.0;
  double wellthick=0.0;
  char pathcontrol[125];
  double deltaCP=0, current_CP=0.0;
  int ic, icl=0, idist=0, flowd=0, welld=0;
	   
  if (disp_o==1)
    {
 
      for (ic=0; ic<time_d-1; ic++)
	{
          sumsquares[ic][0]=0.0;
          sumsquares[ic][1]=0.0;
          nd[ic]=0;
          
          
	  sprintf(filename,"%s/ControlTime_t%d", maindir, ic+1);
	  dis = OpenFile(filename,"w");

          // fprintf(dis,"positions for longitudinal dispersivity calculation at time %f \n", (ic+1)*dtime);
	  fclose(dis);
          
	}

      inputfile=Control_Param("dflowdir:",8);
      flowd=inputfile.param;
    }
	  
		
  // settings for Control Plane / Control Cylinder options	  
  if (out_control==1)
    {	 
     
      inputfile = Control_File("control_out:",12 );
      sprintf(pathcontrol,"%s/%s",maindir,inputfile.filename );     
      
      mkdir(pathcontrol,0777);

    
      inputfile=Control_Param("delta_Control:",14);
      deltaCP=inputfile.param;
      if (out_plane==1)
	{  
	  inputfile=Control_Param("flowdir:",8);
	  flowd=inputfile.param;
    
	  node[nodezonein[0]-1].fracture[0];
    
	  inflowcoord=node[nodezonein[0]-1].coord[flowd];
	  outflowcoord=node[nodezoneout[0]-1].coord[flowd];
     
	  controllength=fabs(outflowcoord)+fabs(inflowcoord);
	  icl=controllength/deltaCP +1;
	}
     
      if (out_cylinder==1)
	{
	  inputfile=Control_Param("welldir:",8);
	  welld=inputfile.param;
     
	  inputfile=Control_Param("lengthtowell:",13);
	  controllength=inputfile.param;
      
	  inputfile=Control_Param("wellthickness:",14);
	  wellthick=inputfile.param;
      
	  if (deltaCP<wellthick)
	    icl=(controllength-wellthick/2)/deltaCP+1; 
	  else
	    icl=controllength/deltaCP +1;
	}
    
     
 
    }
       
  int cross[icl];
  
  
  
    
  /** define particle's initial positions **/
  int numbpart, initweight=0;
  
  /*** set up initial positions of particles ***/

  numbpart=InitPos();     //init_par_line

//  numbpart=InitPos2();     //init_par_line hpham

  //Question: How to get all values of particle array? 
//  int ii;
//  for (ii=0; ii<30; ii++)     //hpham
//  {
//	  printf(" [TrackingPart1.c] np = %4d, fractures = %6d \n", ii, particle[ii].fracture);
//  }


  printf("\n [hpham] %ld time steps in the current run \n", timesteps);
  printf(" [hpham] %d particles actually placed (after flux weighting if selected) \n", numbpart);
  
 
  
  
  
  
  
  
  
  
  
  inputfile=Control_File("flux_weight:",12);
  res=strncmp(inputfile.filename,"yes",3);
  if ((res==0)&&(flag_w==1))
    initweight=1;

//  printf("\n =============== Hello, I am here!!!! ============\n"); 

  
  int i=1, ins=0, intersm=0,curr_n=1;;
  int  t_end=0, fracthit=0; 
  struct posit3d particle3dposit, particle3dvelocity;
 
  int counttimestep=0, prevcell=0, prevfract=0; 
  double xcop, ycop, zcop, currentlength=0, totallength=0;
  double xx=0, yy=0,zz=0, xinit=0, yinit=0, zinit=0;
  
  /**** calculate initial flux weight of particles *****/
  /**** works for first 3 cases of particles seed ******/
         
 
  if (initweight==1)
    {
   
      FlowInWeight(numbpart);
   
    }
  

	
//	  sprintf(filename,"%s/part2well2_%d.out",path,curr_n);
	 sprintf(filename,"part2well2.out");     // hpham
	 p2well2 = OpenFile(filename,"w");	    // hpham
	 fprintf(p2well2,"Par_No, Dist2Well, # of time steps, flux weights, total travel time, x-, y-, z-final pos, beta, total length[m] \n");
   fp_hpham = OpenFile("all_par_locs.txt", "w");
//   fprintf(fp_hpham,"    Par_ID, time,          xloc,     yloc,     zloc,     vel_x,     vel_y,     vel_z  \n");
//   fprintf(fp_hpham,"Par_ID,timestep,partraveltime,xloc,yloc,zloc, rpw \n");
	 fprintf(fp_hpham,"Par_ID,partraveltime,xloc,yloc,zloc,totallength \n");
	// fprintf(fp_hpham,"Par_ID,timestep,partraveltime,xloc,yloc,zloc \n");

  
  /************ LOOP ON PARTICLES  **********/   // hpham
  for (np=0; np<numbpart; np++) 
    {

    if ((numbpart<100)&&(((np-1) % 10)==0)&&(np>2))
		printf("Done %.1f %% of particles. \n", ((float)(np-1)/numbpart)*100);
    else
	if ((numbpart<10000)&&(((np-1) % 100)==0)&&(np>2))
		printf("Done %.1f %% of particles. \n", ((float)(np-1)/numbpart)*100);
	else
	if ((numbpart<100000)&&(((np-1) % 1000)==0)&&(np>2))
		printf("Done %.1f %% of particles. \n", ((float)(np-1)/numbpart)*100);
	else
	if ((((np-1)%10000)==0)&&(np>2))
		printf("Done %.1f %% of particles. \n", ((float)(np-1)/numbpart)*100);
         
    if (avs_o==1)
	{
	  // AVS output (should be optional) 
	  sprintf(filename,"%s/part3D_%d.inp",path,curr_n);
	  wpt = OpenFile(filename,"w");
	  fprintf(wpt,"%10d    %10d    %10d    %10d    %10d\n", timesteps, 0,0,0,0);  
          //open a separate file for attributes, will be attached to the original AVS later
	  sprintf(filename,"%s/part3D_%d.att",path,curr_n);
	  wpt_att = OpenFile(filename,"w");  
	  fprintf(wpt_att,"0008   1    1     1    1    1     1    1\n");
	  fprintf(wpt_att,"fracture, integer\n");
	  fprintf(wpt_att,"time, real\n");
	  fprintf(wpt_att,"velocity, real\n");
	  fprintf(wpt_att,"vel_x, real\n");
	  fprintf(wpt_att,"vel_y, real\n");
	  fprintf(wpt_att,"vel_z, real\n");
	  fprintf(wpt_att,"aperture, real\n");
          fprintf(wpt_att,"pressure, real\n"); 
	} 
       
      if (traj_o==1)
	{   
	  // ascii output of: 3d positions, 3d velocities, cell, fracture, time and beta     
	  sprintf(filename,"%s/traject_%d",path,curr_n);
	  wv = OpenFile(filename,"w");
	  fprintf(wv,"    Current time step, x-, y-, z- pos., Vx, Vy, Vz at ths positions, # of cell, #of fracture, travel time, aperture, beta nters. fracture ID, pressure at particle position \n");

	   // output data on intersections only
	  sprintf(filename,"%s/inters_%d",path,curr_n);
	  wint = OpenFile(filename,"w");
	  fprintf(wint,"     Current traj. length, travel time, x-, y-, z- pos., fracture ID, beta, pressure at particle's position \n");
    

	}
 
      // define capacity for temp data used for outputs
      int capacity= (int) timesteps/10;
      
       if (disp_o=!1)
            time_d=1;

      
      double part_squares[time_d][3];
   
        if (disp_o=!1)
		{
		  for (ic=0; ic<time_d; ic++)
			{
				part_squares[ic][0]=0.0;
				part_squares[ic][1]=0.0;
				part_squares[ic][2]=0.0;
			}
		}
		
		
      // control plane/cylinder output 
      t_end=0;
      intersm=0; 
      /* define an initial cell  */
      ins=0;
//	  printf("\n =============== Hello, I am here!!!! ============ \n"); 
//	  printf(" [TrackingPart1.c] np = %4d, particle[np].fracture = %6d. Now calling InitCell ... \n", np, particle[np].fracture);

      ins=InitCell();


//      printf(" initial cell ins = %d \n", ins);
      
      if (ins==0)
		{
//			printf("Initial cell is not found for particle %4d %8.4f %8.4f %8.4f in fract %4d \n", np+1, particle[np].position[0], particle[np].position[1], particle[np].position[2], particle[np].fracture);
			counthp1++;
		}
      else      // ins=!0 (long loop, upto line 1091)
		{
//	    printf("init pos %5.8e %5.8e \n", particle[np].position[0], particle[np].position[1]); 
	  // set up initial values for Lagrangian variables
        lagvariable.tau=0.0;
        lagvariable.betta=0.0;
          
        prevcell=particle[np].cell;
		prevfract=particle[np].fracture;
        unsigned int fract_id[nfract], id=0;            
		for (id=0; id<nfract; id++)
			fract_id[id]=0;
			
	   fract_id[0]=particle[np].fracture;
	
	  /************LOOP ON TIME *********/   // hpham
	  FLAG_OUT=0;
	  int stuck=0, stuckcell=0, cur_ind=0, cur_node=0;
	  // double delat=0;
	  t=0;
	  nodeID=0;
	  counttimestep=0;
	  xinit=0.0;
	  yinit=0.0;
	  zinit=0.0;
	  fracthit=0;
	  // output particles initial positions	
	  particle3dposit=CalculatePosition3D();      //hpham's mark: This is to get particle positions
	  xcop=particle3dposit.cord3[0];
	  ycop=particle3dposit.cord3[1];
	  zcop=particle3dposit.cord3[2];
           
	  xinit=xcop;
	  yinit=ycop;
	  zinit=zcop;

//	  printf("init pos %8.4f %8.4f %8.4f \n", xcop, ycop, zcop);     //hpham


	  lagvariable.initx=xinit;
	  lagvariable.inity=yinit;
	  lagvariable.initz=zinit;
	  kd=1;
          
	  double tautau=0.0;
	  double beta=0.0;
	  if  (outinit>0)
	    {
//		  printf("Start printing outinit ... \n");     //hpham: delete later. 
	      fprintf(inp,"%05d, %05d, %05d, %5.12E, %5.12E, %5.12E, %5.12E \n",np+1, particle[np].cell, particle[np].fracture, particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2], particle[np].fl_weight); 
	    }
	  totallength=0.0;
	  currentlength=0.0;
	         
	  //counts for control plane/time output
	  
	  FILE* tmp2;
	  
	  if (out_control==1)
	    {	  
	   
	      for  (ic=0; ic<icl; ic++)
		{
	     
		  cross[ic]=0;
		}  
	      idist=0;
	      sprintf(filename,"%s/part_control_%d",pathcontrol, curr_n);
	  
	      tmp2=OpenFile(filename,"w");
	      fprintf(tmp2," travel time, x-, y-, z- position, Vx, Vy, Vz, trajectory length, #  of current fracture, aperture \n");
	     
	      if (out_plane==1)
		{    
		  fprintf(tmp2,  "%5.12E  %5.12E   %5.12E  %5.12E   %5.12E   %5.12E   %5.12E  %5.12E  %05d  %5.12E\n",particle[np].time, xcop, ycop,zcop, 0.0, 0.0,0.0, totallength, particle[np].fracture, node[cell[particle[np].cell-1].node_ind[0]-1].aperture);
	    
	  
		  // in case of inflow is negative 
		  if (inflowcoord<0) 
		    current_CP=inflowcoord+deltaCP;
		  else
		    current_CP=inflowcoord-deltaCP;
	     
		}  
	   
	      if (out_cylinder==1)
		current_CP=controllength;
	     
	    }

          timecounter=0;
	  /////////////////////////  TIME LOOP /////////////////////////////////// 
	         
	  for (t=0; t<timesteps; t++)     // hpham time loop 
	    {
//		  printf("hpham: np = %d | current_ts = %d | timecounter =  %d \n", np, t, timecounter); 
//		  printf("hpham: np = %d | current_ts = %d | particle[np].time =  %d \n", np, t, particle[np].time);
	      particle3dposit=CalculatePosition3D();
           
	      xx=particle3dposit.cord3[0]-xcop;
	      yy=particle3dposit.cord3[1]-ycop;
	      zz=particle3dposit.cord3[2]-zcop;
           
	      currentlength=sqrt(xx*xx+yy*yy+zz*zz);
           
	      totallength=totallength+currentlength;
           
	      xcop=particle3dposit.cord3[0];
	      ycop=particle3dposit.cord3[1];
	      zcop=particle3dposit.cord3[2];    
//          printf("Line 642 \n");   

	      if (no_out!=1)
		{
//          printf("Line 646 \n");      
		  if (timecounter==0)
		    tempdata=(struct tempout*) malloc (capacity*sizeof(struct tempout));

		  if (tempdata==NULL)
		    printf("Allocation memory problem - tempdata\n");    //hpham 
    
//	      printf("Line 653 \n");   
		  particle3dvelocity=CalculateVelocity3D();
//		  printf("timecounter %d time %d capacity %d \n", timecounter, t, capacity);
		  tempdata[timecounter].times=t;   
		  tempdata[timecounter].position2d[0]= particle[np].position[0]; 
		  tempdata[timecounter].position2d[1]= particle[np].position[1]; 
		  tempdata[timecounter].position3d[0]= particle3dposit.cord3[0]; 
		  tempdata[timecounter].position3d[1]= particle3dposit.cord3[1];
		  tempdata[timecounter].position3d[2]= particle3dposit.cord3[2];
		  tempdata[timecounter].velocity3d[0]= particle3dvelocity.cord3[0]; 
		  tempdata[timecounter].velocity3d[1]= particle3dvelocity.cord3[1];
		  tempdata[timecounter].velocity3d[2]= particle3dvelocity.cord3[2];
		  tempdata[timecounter].cellp= particle[np].cell; 
		  tempdata[timecounter].fracturep = particle[np].fracture; 
		  tempdata[timecounter].timep= particle[np].time; 
		  tempdata[timecounter].betap= beta; 
                  tempdata[timecounter].length_t=totallength;
                  tempdata[timecounter].pressure=particle[np].pressure;
                  timecounter++;
//          printf("hpham_TP.c: ipar, xp, yp, ele_ID, frac_ID = %d %f %f %d %d \n",np, particle[np].position[0], particle[np].position[1], particle[np].cell, particle[np].fracture);

		  // if memory should be reallocated
//		  printf("hpham: timecounter_ = %d \n", timecounter);
		  if (timecounter==capacity)
            {
		      capacity=2*capacity;
		      tempdata=(struct tempout*)realloc(tempdata, sizeof(struct tempout)*capacity);
		      if (tempdata==NULL)
			printf("REAllocation memory problem - tempdata\n");
		    }
    
                
		 
		}


	      //calculations for dispersivity: square of (xo-x) is calculated for transverse disersivity only
	      // for longitudinal dispersivity we save the actual coordination of the particle (commented out)
	      double ctime=0.0;
            if (disp_o==1)
              {
	      ctime=kd*dtime;
	      if (((ctime-epsl)<=particle[np].time) && ((ctime+epsl)>=particle[np].time))
		{
		  //                   if (flowd==0)
		  //                     {
		  //                     part_squares[kd-1][0]=(ycop-yinit)*(ycop-yinit);
		  //                     part_squares[kd-1][1]=(zcop-zinit)*(zcop-zinit);
		  //                     part_squares[kd-1][2]=xcop;
		  //                     }
		  //                      if (flowd==1)
		  //                     {
		  //                     part_squares[kd-1][0]=(xcop-xinit)*(xcop-xinit);
		  //                     part_squares[kd-1][1]=(zcop-zinit)*(zcop-zinit);
		  //                     part_squares[kd-1][2]=ycop;
		  //                     }
		  //                      if (flowd==2)
		  //                     {
		  //                    part_squares[kd-1][0]=(xcop-xinit)*(xcop-xinit);
		  //                    part_squares[kd-1][1]=(ycop-yinit)*(ycop-yinit);
		  //                    part_squares[kd-1][2]=zcop;
		  //                    }
		  // 3d position of particle
		  part_squares[kd-1][0]=xcop;
		  part_squares[kd-1][1]=ycop;
		  part_squares[kd-1][2]=zcop;

		  if (kd<time_d) kd++;
		}
	       }
		
	      /***** output data at each control plane ********/		
	      double welldist=0.0; //shortest distance from particle to well
	     
	     
	      if (out_control==1)
		{
		  if (out_plane==1)
	            {   
		      if (( (inflowcoord<0) && (cross[idist]==0) && (particle3dposit.cord3[flowd] >=current_CP)) || ( (inflowcoord>0) && (cross[idist]==0) && (particle3dposit.cord3[flowd] <=current_CP)))
			{
	        
			  cross[idist]=1;
			  //	         printf(" %d cross %d current_CP %lf at z %lf \n",np, idist, current_CP, zcop);  
			  if (no_out==1)
			    particle3dvelocity=CalculateVelocity3D();
	         
			  fprintf(tmp2,"%5.12E  %5.12E   %5.12E  %5.12E   %5.12E   %5.12E   %5.12E  %5.12E  %05d  %5.12E\n",particle[np].time, xcop, ycop,zcop, particle3dvelocity.cord3[0], particle3dvelocity.cord3[1],particle3dvelocity.cord3[2], totallength, particle[np].fracture, node[cell[particle[np].cell-1].node_ind[0]-1].aperture);
	         
			  if (inflowcoord<0)
			    current_CP=current_CP+deltaCP;
			  else
			    current_CP=current_CP-deltaCP;
			  idist=idist+1;
	          
			}
		    }
	         
		  if ((out_cylinder==1)&&(current_CP>=wellthick/2))
	            {
	         
	         
		      if (welld==0)
	         	welldist=sqrt(ycop*ycop+zcop*zcop);
		      if (welld==1)
	         	welldist=sqrt(xcop*xcop+zcop*zcop);
		      if (welld==2)
	         	welldist=sqrt(xcop*xcop+ycop*ycop);	
	         
	   
		      if ((cross[idist]==0) && (welldist<=current_CP)) 
			{
	         
			  //	printf(" %d cross %d current_CP %lf at z %lf \n",np, idist, current_CP, welldist);         
			  cross[idist]=1;
			  if (no_out==1)
			    particle3dvelocity=CalculateVelocity3D();
	         
			  fprintf(tmp2,"%5.12E  %5.12E   %5.12E  %5.12E   %5.12E   %5.12E   %5.12E  %5.12E  %05d  %5.12E\n",particle[np].time, xcop, ycop,zcop, particle3dvelocity.cord3[0], particle3dvelocity.cord3[1],particle3dvelocity.cord3[2], totallength, particle[np].fracture, node[cell[particle[np].cell-1].node_ind[0]-1].aperture);
	         
			  current_CP=current_CP-deltaCP;
			  idist=idist+1;
			}
	            }
	         
		}

		
		
		
	      /*** if particle is on intersection cell **/
	      /*** check distance to intersection ***/ 
	      
	      if ((particle[np].intcell==1) || (particle[np].intcell==3)) 
		{
		  intersm=CheckDistance (); 
		  prevcell=particle[np].cell;
		 
		}

	      /*** if particle's new cell was not found move to the next particle ***/ 
	      if (particle[np].cell==0)
		{
	      
		  break;
                }
    
	      /** Get new particle's velocity ***/
	      
	     
              PredictorStep();
	      //                   printf("%5.8e % 5.8e %d \n", particle[np].velocity[0], particle[np].velocity[1], particle[np].cell);
	      //             printf(" %5.8e %5.8e %5.8e %5.8e %5.8e %5.8e \n",node[5270-1].velocity[0][0],node[5270-1].velocity[0][1], node[5271-1].velocity[0][0],node[5271-1].velocity[0][1],node[5286-1].velocity[0][0],node[5268-1].velocity[0][1]); 
	      /** Calculate new weights and check if particle is in new cell ***/
 
	      CheckNewCell();
	      if ((FLAG_OUT==1))
		break;

	      if (particle[np].cell!=0) 
		{
		  if ((particle[np].intcell!=1)&&(particle[np].intcell!=3))
		    /** Get new particle's position ***/
		    //delat=CalculateCurrentDT();
              
		    CorrectorStep();
		  lagvariable = CalculateLagrangian(xcop,ycop,zcop,lagvariable.initx, lagvariable.inity, lagvariable.initz);
		  tautau=tautau+lagvariable.tau;
		  beta=beta+lagvariable.betta;
		  particle[np].time=particle[np].time+lagvariable.tau;
		}
	      else
		{
		  if ((FLAG_OUT!=1)&&(FLAG_OUT!=3))
		    {
		      //	    printf("%5.8e % 5.8e %d \n", particle[np].velocity[0], particle[np].velocity[1], particle[np].cell);
		      //	      		    printf("cell=0 %d for particle %d at time %d after predictor step \n",particle[np].cell, np+1, t );
		      
		    }
		  else
		    {
		      if (FLAG_OUT==1)
			break;
		      else
			if ((particle[np].cell==0)&& (FLAG_OUT==3))
			  {
			    Moving2NextCellBound(prevcell);
			    FLAG_OUT=0;
			  }
		    }
		}
	      /** Calculate new weights and check if particle is in new cell ***/
              if (particle[np].cell!=0)        
		CheckNewCell();
	   
	      if ((particle[np].cell==0)&& (FLAG_OUT==3))
                {
		  FLAG_OUT=0;
		  Moving2NextCellBound(prevcell);
		}
       
	      if ((particle[np].cell!=prevcell) && (particle[np].cell!=0))
		{
		  /*** calculate Lagrangian variables tau and beta *****/	
		  /*** the sum accumulates every time particle crosses triangular edge***/
	
		  //		lagvariable = CalculateLagrangian(xcop,ycop,zcop,lagvariable.initx, lagvariable.inity, lagvariable.initz);
		  //		tautau=tautau+lagvariable.tau;
		  //		beta=beta+lagvariable.betta;
		
		
		  //		printf("current cell  %d  init pos %lf %lf %lf \n",particle[np].cell, lagvariable.initx, lagvariable.inity, lagvariable.initz);
		  //		printf ("tau= %lf betta = %lf \n",tautau, beta);
		
		  counttimestep=0;
	          prevcell=particle[np].cell;
		  if (prevfract!=particle[np].fracture)
		    {  
		      stuck=0;
		      stuckcell=0;
		      cur_ind=0;
		      cur_node=0;
		      prevfract=particle[np].fracture;
		      fracthit=fracthit+1;
                      fract_id[fracthit]=particle[np].fracture;

		      //		   printf(" %d %d %15.8e \n",np+1, prevfract, particle[np].time); 
		    }
		}
        
	      if (particle[np].cell==prevcell)
		{
		  counttimestep++;
		  if (counttimestep > (int)timesteps/3.0)
		    {
		      //	                         printf("stuck %05d %5.12E %5.12E %5.12E %5.12E %5.12E %5.12E %5.12E %5.12E %05d %05d %5.12E %d %d\n",t+1,particle[np].position[0], particle[np].position[1], particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2],particle3dvelocity.cord3[0], particle3dvelocity.cord3[1],particle3dvelocity.cord3[2],particle[np].cell, particle[np].fracture, particle[np].time, counttimestep, (int)timesteps/3.0);
                  
		      break;
		      if (particle[np].cell==stuckcell)
			stuck=cur_ind;
		      else
			{
                          stuck=0; 
                          stuckcell=particle[np].cell;
			} 
		      if (cur_ind<node[cell[particle[np].cell-1].node_ind[cur_node]-1].numneighb)        
			cur_ind=Moving2NextCell(stuck, cell[particle[np].cell-1].node_ind[cur_node]);
		      else
			{
			  if (cur_node<2)
			    {
			      cur_ind=0;
			      cur_node++;
			      cur_ind=Moving2NextCell(stuck, cell[particle[np].cell-1].node_ind[cur_node]);
			    }
			  else
			    {
			      //		      			      printf(" got stuck \n");
			      break;
			    }
			}

		      /** time correction step ***/    
		      double delt=CalculateCurrentDT();
		      particle[np].time=particle[np].time-counttimestep*delt;
		      //		      printf(" cot %15.8e %d\n",counttimestep*delt, counttimestep); 
		      t=t-counttimestep;
		    }
		}
	      t_end=t;

	      /*** if particle's new cell was not found move to the next particle ***/ 
	      if (particle[np].cell==0)
		{    
		  break;
                }

		 // ================================================================================== 
		// == hpham added these lines to print out particles to wells ======================
		 if (opt_mwell_on==1)
		   {
         double xpw, ypw;
         int curr_time;
		//	  int obs_well_o;
		//	  obs_well_o=1;  // generalize this later. 
		//      if (obs_well_o==1)
	    //	 { 
		  xpw=fabs(x0w-particle3dposit.cord3[0]);  // x_distance of particle to well
		  ypw=fabs(y0w-particle3dposit.cord3[1]);  // y_distance of particle to well
		  rpw=sqrt(xpw*xpw+ypw*ypw);

//	      printf("%8d %8d %12.3E %12.3E %12.3E %12.3E %8.3f \n", np+1, pfrac, time, posit[0], posit[1], posit[2], rpw); 
//            printf("%d %6.2f %d  %5.12E  %5.12E  %5.12E  %5.12E  %5.12E  %5.12E  %5.12E \n",np, rpw, t_end, particle[np].fl_weight, particle[np].time, particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2], beta, totallength);     //hpham
//              fprintf(fp_hpham,"%d %6.2f %d  %5.12E  %5.12E  %5.12E  %5.12E  %5.12E  %5.12E  %5.12E \n",np, rpw, t_end, particle[np].fl_weight, particle[np].time, particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2], beta, totallength);     //hpham
//              printf("fmod = %8.6f \n", fmod(round(totallength*10000)/10000,0.01));

		if ( particle[np].time <= par_prt_cutoff )  // hpham: Only consider particles travel in < 0.3425 years (~3000 hours)
		  {
//        if ( (t % par_ts_freq == 0) && (np % par_np_freq == 0) && (rpw > rwell) )     // print particle location every 100 timesteps; do not print time > xxx years
        curr_time = ceil(particle[np].time*3600);
        fprintf("Current time = %d\n", curr_time);
				if ( (curr_time % par_ts_freq == 0) && (np % par_np_freq == 0) && (rpw > rwell) )     // print particle location every 100 timesteps; do not print time > xxx years
						{
//              fprintf(fp_hpham,"%8d, %8d, %12.6f, %6.4f, %6.4f, %6.4f, %6.4f \n", np+1, t, particle[np].time, particle3dposit.cord3[0], particle3dposit.cord3[1], particle3dposit.cord3[2], rpw);     //hpham
							fprintf(fp_hpham,"%8d, %12.6f, %6.4f, %6.4f, %6.4f, %12.6f  \n", np+1, particle[np].time, particle3dposit.cord3[0], particle3dposit.cord3[1], particle3dposit.cord3[2], totallength);     //hpham
		//				    printf("%8d, %8d, %5.12E, %6.4f, %6.4f, %6.4f, %6.4f \n", np+1, t, particle[np].time, particle3dposit.cord3[0], particle3dposit.cord3[1], particle3dposit.cord3[2], rpw);     //hpham

				  }
					  if ( (rpw <= rwell) && (bot_well<=particle3dposit.cord3[2]) && (particle3dposit.cord3[2]<=top_well) )   // rwell is the radius of the monitoring well
						{
						 fprintf(p2well2,"%8d, %6.2f, %10d, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f \n",np+1, rpw, t_end, particle[np].fl_weight, particle[np].time, particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2], beta, totallength);
		//                 printf("np = %d | rpw = %6.2f | t_end = %d | weight = %5.12E | par_time = %5.12E | %5.12E  %5.12E  %5.12E  %5.12E  %5.12E \n",np, rpw, t_end, particle[np].fl_weight, particle[np].time, particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2], beta, totallength);     //hpham
		//	             printf("%d %6.2f %6.2f %6.2f %6.2f \n", np, rpw, particle3dposit.cord3[0], particle3dposit.cord3[1], particle3dposit.cord3[2]);
						 break;    // hpham: This ends the time loop and moves to a next particle.
		//				  printf("%8d %8d %12.3E %12.3E %12.3E %12.3E %8.3f \n", np+1, pfrac, time, posit[0], posit[1], posit[2], rpw); 
		//				  printf("%8d %8d %12.3E %12.3E %12.3E %12.3E %8.3f \n", np+1, particle[np].fracture, tempdata[i].timep, particle3dposit.cord3[0], particle3dposit.cord3[1], particle3dposit.cord3[2], rpw); 
						}
					//	 } 
		  }
		  else
		  {
			break;
		  }



		// == end hpham_added ================================================================ 
		// ===================================================================================				
		   }
//		rpw_sum = rpw_sum + fabs(rpw-rpw0);
//		rpw0=rpw;
//		hp_count=hp_count+1;  
////		printf("particle_id=%d, hp_count=%d, rpw=%8.5f, rpw_sum=%8.5f. \n", np,hp_count,rpw,rpw_sum);
//		if (rpw_sum < 1E-12)
//		{
//			printf("particle id %d is strapped at distantce = %8.5f (m). \n", np,rpw);
//			break;
//		}
//		if (hp_count>1000)
//		{
//			hp_count=0;
//			rpw_sum=0;
//		}
		   
		   
	     } //loop on time hpham
 
//			  printf("hpham: np = %d | t = %d \n", np+1, t);


	  /** if particle did not go out through flow-out zone ****/ 
	  if (FLAG_OUT!=1)
	    {
	 
	
	      //	  	  printf("Done for particle %d in fracture %d at time %d, located at %f, %f; cell %d, Flag %d intc %d.  \n\n", np+1, particle[np].fracture, t_end,particle[np].position[0], particle[np].position[1], prevcell, FLAG_OUT, particle[np].intcell);  
	  
	    
	      if (out_control==1)
		{	
		  fclose(tmp2);  
		}		  
	    }
	  else
	    {  
	      /** if particle went out through flow out zone ****/
	      if (FLAG_OUT==1)
		{
		  if (no_out!=1)
	  
		    {
		      if (particle[np].cell !=0)
			{
		    
			  FinalPosition();
			}
		      particle3dposit=CalculatePosition3D();
		      particle3dvelocity=CalculateVelocity3D();

		      tempdata[timecounter].times=t;
		      tempdata[timecounter].position2d[0]= particle[np].position[0];
		      tempdata[timecounter].position2d[1]= particle[np].position[1];
		      tempdata[timecounter].position3d[0]= particle3dposit.cord3[0];
		      tempdata[timecounter].position3d[1]= particle3dposit.cord3[1];
		      tempdata[timecounter].position3d[2]= particle3dposit.cord3[2];
		      tempdata[timecounter].velocity3d[0]= particle3dvelocity.cord3[0];
		      tempdata[timecounter].velocity3d[1]= particle3dvelocity.cord3[1];
		      tempdata[timecounter].velocity3d[2]= particle3dvelocity.cord3[2];
		      tempdata[timecounter].cellp= particle[np].cell;
		      tempdata[timecounter].fracturep = particle[np].fracture;
		      tempdata[timecounter].timep= particle[np].time;
		      tempdata[timecounter].betap= beta;
		      tempdata[timecounter].length_t=totallength;
		      tempdata[timecounter].pressure=particle[np].pressure;
		    
			  
		      ParticleOutput(t, 0);
            
	  
		    }
		  else
		    {
		      if (particle[np].cell !=0)
			{
		    
			  FinalPosition(); 
			}         
		      particle3dposit=CalculatePosition3D();
		  
		  
		    }
		  //adding data to dispersivity


		  for (ic=0; ic<kd-1; ic++)
		    {
		      nd[ic]++;

		      //            sumsquares[ic][0]=sumsquares[ic][0]+part_squares[ic][0];
		      //            sumsquares[ic][1]=sumsquares[ic][1]+part_squares[ic][1];
		      if (disp_o==1)
			{
			  sprintf(filename,"%s/dispers_t%d", maindir, ic+1);
			  dis = OpenFile(filename,"r+");
			  fseek(dis, 0, SEEK_END);
			  fprintf(dis,"%f %f %f \n", part_squares[ic][0], part_squares[ic][1], part_squares[ic][2]);
			  //        printf("%f %d part %d time_d %d\n",part_squares[ic][2], ic, np, time_d);
			  fclose(dis);
			}
		    }
	  
		  if (particle[np].cell!=0)
		    {
		      lagvariable = CalculateLagrangian(xcop,ycop,zcop,lagvariable.initx, lagvariable.inity, lagvariable.initz);
		      tautau=tautau+lagvariable.tau;
		      beta=beta+lagvariable.betta;
		    }
           
		  xx=particle3dposit.cord3[0]-xcop;
		  yy=particle3dposit.cord3[1]-ycop;
		  zz=particle3dposit.cord3[2]-zcop;
           
		  currentlength=sqrt(xx*xx+yy*yy+zz*zz);
           
		  totallength=totallength+currentlength;
	  
	  
	   
		  /**  output travel time *****/ 
		  fprintf(tp,"%d  %5.12E  %5.12E  %5.12E  %5.12E  %5.12E  %5.12E  %5.12E \n",t_end, particle[np].fl_weight, particle[np].time, particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2], beta, totallength);
		  if (tort_o>0)
		    {
		      fprintf(tort,"%5.12E %5.12E %5.12E %5.12E %5.12E %5.12E %5.12E  %d\n",totallength, xinit, yinit, zinit, particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2], fracthit);
		    }   


                   if (frac_o>0)
                   { 
                   id=0;
                     do
                       { 
                      fprintf(frac," %5d ", fract_id[id]);         
                      id++;
                        }
                     while (fract_id[id]!=0);
                       fprintf(frac,"\n"); 
            
                   }
		  //      printf("%d %5.12E %5.12E %5.12E\n", np, tautau, particle[np].time, beta );
		  //	  printf("%d %5.12E %5.12E %5.12E %5.12E %5.12E \n",t_end, particle[np].fl_weight, particle[np].time, particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2]);
	      
	     
		  // output of last control plane - outflow plane    
          
		  if (out_control==1)
		    {
           
		      if (particle[np].cell==0)
			particle[np].cell=prevcell;
           
		      particle3dvelocity=CalculateVelocity3D();
          
            
		      fprintf(tmp2,"%5.12E  %5.12E   %5.12E  %5.12E   %5.12E   %5.12E   %5.12E  %5.12E  %05d  %5.12E\n",particle[np].time, particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2], particle3dvelocity.cord3[0], particle3dvelocity.cord3[1],particle3dvelocity.cord3[2], totallength, particle[np].fracture, node[cell[particle[np].cell-1].node_ind[0]-1].aperture);
           
           
		      fclose(tmp2);
           
	        
		    }    
	        
	              
		  curr_n++;
	  
		  if (avs_o==1)
		    {
		      /*** write a connectivity list in inp files ***/
		      for (i=0; i<nodeID-1; i++)
			{
			  fprintf(wpt, "%10d %5d line %10d %10d\n", i+1, 1, i+1,i+2);  
			}
		      rewind(wpt);
		      fprintf(wpt,"%10d    %10d    %10d    %10d    %10d\n", nodeID, nodeID-1 ,8,0,0); 
		      fclose(wpt);
		      fclose(wpt_att);
		      char filename1[125]={0}, filename2[125]={0}, filename3[125]={0}, buffer[500]={0};
                        char wspace=' ';
		      sprintf(filename1,"%s/part3D_%d.inp",path,curr_n-1);
   		      sprintf(filename2,"%s/part3D_%d.att",path,curr_n-1);
		      sprintf(filename3,"%s/part_%d.inp",path,curr_n-1);

		      sprintf(buffer,"cat%c%s%c%s%c>%c%s",wspace,filename1, wspace,filename2,wspace,wspace,filename3); 
		      
			system (buffer);
		      sprintf(buffer,"rm%c-f%c%s%c%s  ",wspace,wspace,filename1,wspace,filename2); 
		      system (buffer);
		
		    }
	   
		}
	    }
	} //end if ins!=0 (the initial cell was found)

	
      if (traj_o==1)
	{
          rewind(wv);
          fprintf(wv,"%d \n",nodeID);
	  fclose(wv);

          fclose(wint); 
	
	}

     
    } //end of particle loop  hpham



  printf("\n  NOT FOUND initial cells for %d particle. \n", counthp1); 
  printf("\n  Number of particles that went out through flow-out boundary: %d \n", curr_n-1); 
  fclose(p2well2);    // hpham    
  fclose(tp);
  fclose(fp_hpham);
  if (tort_o>0)
    fclose(tort);
  if (outinit>0)
    fclose(inp);
  if (frac_o>0)
     fclose(frac); 
    
  sprintf(filename,"%s/TotalNumberP",maindir);
  FILE *tn = OpenFile (filename,"w");
  fprintf(tn," %10d \n", curr_n-1);
  fclose(tn);  
  
  if (marfa==1 || plumec==1 )
    {
      printf("\n    Working on additional outputs    \n");
      OutputMarPlumDisp (curr_n-1, path);
    }

    
  return;
}
/////////////////////////////////////////////////////////////////////////////

void CheckNewCell()
{
  /******Function looks for the current / new cell for particle************/
  struct intcoef lambda;
  int n1=0, n2=0, n3=0, pcell; 
  double delta_t;
  int cb=0;
  n1=cell[particle[np].cell-1].node_ind[0];
  n2=cell[particle[np].cell-1].node_ind[1];
  n3=cell[particle[np].cell-1].node_ind[2];
     
  delta_t=CalculateCurrentDT(); 

   
  /** first, calculate weights in the current cell****/

  lambda=CalculateWeights(n1, n2, n3);

  if ((lambda.weights[0]<=1.)&&(lambda.weights[0]>=0.)&&(lambda.weights[1]<=1.)&&(lambda.weights[1]>=0.)&&(lambda.weights[2]<=1.)&&(lambda.weights[2]>=0.))
    {
      /** particle is in the current cell ***/
      particle[np].weight[0]=lambda.weights[0];
      particle[np].weight[1]=lambda.weights[1];
      particle[np].weight[2]=lambda.weights[2];
    }
  else
    {
      /* particle is not in the current cell ***/
      pcell=particle[np].cell;
      int pfract;
      pfract=particle[np].fracture;
      particle[np].cell=0;

      SearchNeighborCells(n1,n2,n3);
      cb=0;
      if ((node[n1-1].typeN==210)||(node[n1-1].typeN==212)||(node[n1-1].typeN==200)||(node[n1-1].typeN==202))
	cb=cb+1;
      if ((node[n2-1].typeN==210)||(node[n2-1].typeN==212)||(node[n2-1].typeN==200)||(node[n2-1].typeN==202))
	cb=cb+1;
      if ((node[n3-1].typeN==210)||(node[n3-1].typeN==212)||(node[n3-1].typeN==200)||(node[n3-1].typeN==202))
	cb=cb+1;
      if (cb!=0)
	{
	  //   printf(" Particle %d IS OUT of flow out zone. \n", np+1);
	  FLAG_OUT=1;
	} 
      if (particle[np].cell==0)
	{
	  cb=0;
	  
	  if (((node[n1-1].typeN==10)||(node[n1-1].typeN==12)))
	    cb=cb+1;
	  if (((node[n2-1].typeN==10)||(node[n2-1].typeN==12)))
	    cb=cb+1;
	  if (((node[n3-1].typeN==10)||(node[n3-1].typeN==12)))
	    cb=cb+1;
	  if ((cb>0)&&(node[n1-1].typeN!=210)&&(node[n1-1].typeN!=212)&&(node[n2-1].typeN!=210)&&(node[n2-1].typeN!=212)&&(node[n3-1].typeN!=210)&&(node[n3-1].typeN!=212))
	    {
	      //       printf("Particle %d is out of fracture %f  %f %d %d\n", np+1, particle[np].position[0], particle[np].position[1], particle[np].fracture, pcell);

	      /** if out of fracture -   make a flip in x direction ***/ 

	      double  cx, cy;
	      cy=1.0;
	      cx=-1.0;
	      particle[np].position[0]=particle[np].prev_pos[0]+delta_t*particle[np].velocity[0]*cx;
	      particle[np].position[1]=particle[np].prev_pos[1]+delta_t*particle[np].velocity[1]*cy;
	      SearchNeighborCells(n1,n2,n3);
            
	      if (particle[np].cell==0)
		{
		  /** if still  out of fracture -   make a flip in y direction ***/ 
		  cy=-1.0;
		  cx=1.0;
		  particle[np].position[0]=particle[np].prev_pos[0]+delta_t*particle[np].velocity[0]*cx;
		  particle[np].position[1]=particle[np].prev_pos[1]+delta_t*particle[np].velocity[1]*cy;
		  SearchNeighborCells(n1,n2,n3);
		}
		
	      if (particle[np].cell==0)
		{
		  FLAG_OUT=3;
		   
		  //		     particle[np].cell=pcell;
		  //             printf(" before %f %f \n", particle[np].position[0], particle[np].position[1]); 
                     
                     
             
		  //		     Moving2Center(np, pcell);
		  //		      printf(" afteer %f %f \n", particle[np].position[0], particle[np].position[1]); 
		  // 		   printf("particle %d was moved to cell center %d  time \n", np+1, particle[np].cell);
		} 
	    }
	  cb=0;
	  if ((node[n1-1].typeN==210)||(node[n1-1].typeN==212)||(node[n1-1].typeN==200))
	    cb=cb+1;
	  if ((node[n2-1].typeN==210)||(node[n2-1].typeN==212)||(node[n2-1].typeN==200))
	    cb=cb+1;
	  if ((node[n3-1].typeN==210)||(node[n3-1].typeN==212)||(node[n3-1].typeN==200))
	    cb=cb+1;
	  if (cb!=0)
	    {
	      //   printf(" Particle %d IS OUT of flow out zone. \n", np+1);
	      FLAG_OUT=1;
	    } 
	  else
	    {
	      int ii=0;
	      for (ii=0; ii<node[n1-1].numneighb; ii++)
		{
		  if ((node[node[n1-1].indnodes[ii]-1].typeN==212)||(node[node[n1-1].indnodes[ii]-1].typeN==210))
		    {
		      cb++;
		      FLAG_OUT=1;
		      break;
		    }
		}
	      if (cb==0)
		{
		  for (ii=0; ii<node[n2-1].numneighb; ii++)
		    {
		      if ((node[node[n2-1].indnodes[ii]-1].typeN==212)||(node[node[n2-1].indnodes[ii]-1].typeN==210))
			{
			  cb++;
			  FLAG_OUT=1;
			  break;
			}
		    }
		} 
	      if (cb==0)
		{
		  for (ii=0; ii<node[n3-1].numneighb; ii++)
		    {
		      if ((node[node[n3-1].indnodes[ii]-1].typeN==212)||(node[node[n3-1].indnodes[ii]-1].typeN==210))
			{
			  cb++;
			  FLAG_OUT=1;
			  break;
			}
		    }
		} 
	     
	    } 
	}
    }
  return;
}
/////////////////////////////////////////////////////////////////////////////

struct intcoef  CalculateWeights(int nn1, int nn2, int nn3)
{
  /****** function calculates weights for velocity interpolation *****************/
  double n1x,n1y,n2x,n2y,n3x,n3y, deter;
  struct intcoef lambda;
   double eps=10e-5;


  n1x=node[nn1-1].coord_xy[Xindex(nn1,np)];
  n1y=node[nn1-1].coord_xy[Yindex(nn1,np)];
     
  n2x=node[nn2-1].coord_xy[Xindex(nn2,np)];
  n2y=node[nn2-1].coord_xy[Yindex(nn2,np)];
 
  n3x=node[nn3-1].coord_xy[Xindex(nn3,np)];
  n3y=node[nn3-1].coord_xy[Yindex(nn3,np)];
 
  deter=(n2y-n3y)*(n1x-n3x)+(n3x-n2x)*(n1y-n3y);

  lambda.weights[0]=((n2y-n3y)*(particle[np].position[0]-n3x)+(n3x-n2x)*(particle[np].position[1]-n3y))/deter;
  lambda.weights[1]=((n3y-n1y)*(particle[np].position[0]-n3x)+(n1x-n3x)*(particle[np].position[1]-n3y))/deter;
  
  if ((lambda.weights[0]>=-eps) && (lambda.weights[0]<=eps))
    lambda.weights[0]=0.0;
  if ((lambda.weights[1]>=-eps) && (lambda.weights[1]<=eps))
    lambda.weights[1]=0.0;
  lambda.weights[2]=1-lambda.weights[0]-lambda.weights[1];
  if ((lambda.weights[2]>=-eps) && (lambda.weights[2]<=eps))
    lambda.weights[2]=0.0;

  return lambda;
}
////////////////////////////////////////////////////////////////////////////////
/************** Searching particle in Neighboring Cells ***********************/
void SearchNeighborCells(int nn1, int nn2, int nn3)
{
 

  int k=0;
   
  /*search for the node with highest weight - highest probability that particle is in one of it's neighbors*/
  if ((particle[np].weight[0]>=particle[np].weight[1]) && (particle[np].weight[0]>=particle[np].weight[2]))
    k=nn1;

  if ((particle[np].weight[1]>=particle[np].weight[0]) && (particle[np].weight[1]>=particle[np].weight[2]))
    k=nn2;

  if ((particle[np].weight[2]>=particle[np].weight[1]) && (particle[np].weight[2]>=particle[np].weight[0]))
    k=nn3;

  if (k!=0)

    NeighborCells (k);

  /* if not found, search between neighbors of other two nodes */ 
  if ((particle[np].cell==0)&&(k!=nn1))
    NeighborCells (nn1);  

        
  if ((particle[np].cell==0)&&(k!=nn2))
    NeighborCells (nn2);  


  if ((particle[np].cell==0)&&(k!=nn3))
    NeighborCells (nn3);  

  return;
}
//////////////////////////////////////////////////////////////////////////////
/**** function returns 1 if particle is inside cell ***************/

int InsideCell (int numc)
{ 
  struct intcoef lambda;
  int nk_1, nk_2, nk_3, inside, nb=0;
  nk_1=cell[numc-1].node_ind[0];
  nk_2=cell[numc-1].node_ind[1];
  nk_3=cell[numc-1].node_ind[2];
  double eps=1e-7;
     
  lambda=CalculateWeights(nk_1, nk_2, nk_3);
  int intc=0;
  if (particle[np].intcell==4)
    intc=1;
  
       
  if((lambda.weights[0]<=1.+eps)&&(lambda.weights[0]>=-eps)&&(lambda.weights[1]<=1.+eps)&&(lambda.weights[1]>=-eps)&&(lambda.weights[2]<=1.+eps)&&(lambda.weights[2]>=-eps))
    {
      inside=1;
      particle[np].weight[0]=fabs(lambda.weights[0]);
      particle[np].weight[1]=fabs(lambda.weights[1]);
      particle[np].weight[2]=fabs(lambda.weights[2]);
      particle[np].cell=numc;

      /* particle.intcell is a flag of particle being in intersection (=1) or boundary(=2) cell */

      if (node[nk_1-1].typeN==10) 
	nb++;
      if (node[nk_2-1].typeN==10) 
	nb++;
      if (node[nk_3-1].typeN==10) 
	nb++;
	

      if (nb>1)
	particle[np].intcell=2;
      else
      	particle[np].intcell=0;
 
      if (((node[nk_1-1].typeN==2) || (node[nk_2-1].typeN==2)||(node[nk_3-1].typeN==2) ||(node[nk_1-1].typeN==12) || (node[nk_2-1].typeN==12)||(node[nk_3-1].typeN==12))&& (intc==0))

	particle[np].intcell=1;
 
      if ((nb>0) && (particle[np].intcell==1))
	particle[np].intcell=3;
            
    }
  else
    inside=0;
  return inside;
}
//////////////////////////////////////////////////////////////////////////////
/**** function calculates new velocities and new particle positions ***********/
void PredictorStep()

{
  int n1=0, n2=0, n3=0, v1=0,v2=0,v3=0; 
  double delta_t;
  n1=cell[particle[np].cell-1].node_ind[0];
  n2=cell[particle[np].cell-1].node_ind[1];
  n3=cell[particle[np].cell-1].node_ind[2];

  v1=cell[particle[np].cell-1].veloc_ind[0];
  v2=cell[particle[np].cell-1].veloc_ind[1];
  v3=cell[particle[np].cell-1].veloc_ind[2];

  delta_t=CalculateCurrentDT();
  
  // particle[np].time=particle[np].time+delta_t; 

  /** velocity interpolation ***/
  particle[np].velocity[0]=particle[np].weight[0]*node[n1-1].velocity[v1][0]+particle[np].weight[1]*node[n2-1].velocity[v2][0]+particle[np].weight[2]*node[n3-1].velocity[v3][0];
  particle[np].velocity[1]=particle[np].weight[0]*node[n1-1].velocity[v1][1]+particle[np].weight[1]*node[n2-1].velocity[v2][1]+particle[np].weight[2]*node[n3-1].velocity[v3][1];
   
 particle[np].pressure=particle[np].weight[0]*node[n1-1].pressure+particle[np].weight[1]*node[n2-1].pressure+particle[np].weight[2]*node[n3-1].pressure;
       
  particle[np].prev_pos[0]=particle[np].position[0];
  particle[np].prev_pos[1]=particle[np].position[1];

  particle[np].position[0]=particle[np].position[0]+delta_t*particle[np].velocity[0];
  particle[np].position[1]=particle[np].position[1]+delta_t*particle[np].velocity[1];

  return;
}

//////////////////////////////////////////////////////////////////////////////
/**** function calculates new velocities and new particle positions ***********/
void CorrectorStep()

{
  int n1=0, n2=0, n3=0, v1=0,v2=0,v3=0; 
  double delta_t;
  n1=cell[particle[np].cell-1].node_ind[0];
  n2=cell[particle[np].cell-1].node_ind[1];
  n3=cell[particle[np].cell-1].node_ind[2];

  v1=cell[particle[np].cell-1].veloc_ind[0];
  v2=cell[particle[np].cell-1].veloc_ind[1];
  v3=cell[particle[np].cell-1].veloc_ind[2];

 
  delta_t=CalculateCurrentDT();

  //   particle[np].time=particle[np].time+delta_t; 

  particle[np].position[0]=particle[np].prev_pos[0]+delta_t*particle[np].velocity[0];
  particle[np].position[1]=particle[np].prev_pos[1]+delta_t*particle[np].velocity[1];
  return;
}
///////////////////////////////////////////////////////////////////////////////
/**** function check neighboring cells to find a particle *************/
void NeighborCells (int k)
{
  int i=0,j,inscell,nc;
  do
    {
      for (j=0; j<4; j++)
	{
          if ((node[k-1].fracts[i][j]==particle[np].fracture))
	    {
	      nc=node[k-1].cells[i][j];  
	      inscell=InsideCell(nc);     // function returns 1 if particle is inside cell
	    }
	}
      i++;
    }
  while ((inscell==0)&&(i<node[k-1].numneighb));
  return;
}
//////////////////////////////////////////////////////////////////////////////



void DefineTimeStep()
{
  /**** function defines time step as function of polygon volume and velocity ****/
  int i; 
  double  dotvel1, epsd=1e-10;

  //  FILE *wtt = OpenFile ("timestep","w");

  for (i=0; i<nnodes; i++)
    {
      short int  j=0;
      for (j=0; j<4; j++)
	{
	  
	  dotvel1=node[i].velocity[j][0]*node[i].velocity[j][0]+node[i].velocity[j][1]*node[i].velocity[j][1];
	  if (dotvel1>epsd)  
	    {
	    
	      if (node[i].typeN==10)
	         
		node[i].timestep[j]=0.005*sqrt(((node[i].pvolume)/node[i].aperture)/dotvel1);
	      else 
		node[i].timestep[j]=0.005*sqrt(((node[i].pvolume)/node[i].aperture)/dotvel1);
	    }
	  else
	    node[i].timestep[j]=0.005*sqrt(((node[i].pvolume)/node[i].aperture)/epsd);
 
	  //	  fprintf(wtt,"%d  %d %5.8e %5.8e %5.8e\n", i+1, j+1, node[i].timestep[j],node[i].velocity[j][0],node[i].velocity[j][1]); 
	}
	
    }
  //   fclose(wtt);
  return;
}     // hpham: end of DefineTimeStep
//////////////////////////////////////////////////////////////////////////////
/************ Function checks the distance from particle to intersection line *****************/

int CheckDistance()
{
  /*** define the intersection segment/line ***/
  int dn1, dn2, dn3, int1=0, int2=0, int3=0, i, intm=0;
  int ind_int2, fract_p;
  double px, py, dist, delta_t;
  double cx1=0,cy1=0,cx2=0,cy2=0;
  struct intcoef lambda;

  dn1=cell[particle[np].cell-1].node_ind[0];
  dn2=cell[particle[np].cell-1].node_ind[1];
  dn3=cell[particle[np].cell-1].node_ind[2];
  delta_t=CalculateCurrentDT();

  int prevcell=0;
  prevcell=particle[np].cell;

  px=particle[np].position[0];
  py=particle[np].position[1];

  /* distance that particle will make during next step */

  dist=sqrt(pow((particle[np].velocity[0]*delta_t),2)+pow((particle[np].velocity[1]*delta_t),2));

  /* check: one edge of cell belongs to intersection line */
  if ((node[dn1-1].typeN==2)||(node[dn1-1].typeN==12))
    int1=dn1;
  if ((node[dn2-1].typeN==2)||(node[dn2-1].typeN==12))
    { 
      if (int1!=0)
	int2=dn2;
      else
	int1=dn2;           
    }
  if ((node[dn3-1].typeN==2)||(node[dn3-1].typeN==12))
    {
      if (int1!=0)
	int2=dn3;
      else
	int1=dn3;
    }
  if ((int1!=0) && (int2!=0))
    int3=-1;

  /* check: only one node of cell belongs to intersection - find the second node in neighboring list */
  if (int2==0) 
    {
      for(i=0; i<node[int1-1].numneighb; i++)
	{
          if (((node[int1-1].type[i]==2)|| (node[int1-1].type[i]==12)) && ((node[node[int1-1].indnodes[i]-1].fracture[0]==particle[np].fracture)||(node[node[int1-1].indnodes[i]-1].fracture[1]==particle[np].fracture)))
	    {
	      if (int2==0)
		int2=node[int1-1].indnodes[i];
              else
		int3=node[int1-1].indnodes[i];
            }
          if ((int3!=0) &&(int3!=-1))
	    break;
	}
    }
    
  //       printf("int1 %d int2 %d int3 %d cell %d frac %d \n", int1, int2, int3, particle[np].cell, particle[np].fracture);

  /*** if two nodes on intersection are found. the intersection edge will be defined  *****/
  if ((int1!=0) && (int2!=0))
    {
      cx1=node[int1-1].coord_xy[Xindex(int1,np)];
      cy1=node[int1-1].coord_xy[Yindex(int1,np)];

      cx2=node[int2-1].coord_xy[Xindex(int2,np)];
      cy2=node[int2-1].coord_xy[Yindex(int2,np)];

   
      if ((int1!=0) && (int2!=0) && (int3!=0))
	{
	
	  /* define height of triangle, where the base is intersection segment */
	  double base, side1, side2, perim, area, height;   
	  base=sqrt(((cx2-cx1)*(cx2-cx1))+((cy2-cy1)*(cy2-cy1)));
	  side2=sqrt(((cx2-px)*(cx2-px))+((cy2-py)*(cy2-py)));
	  side1=sqrt(((cx1-px)*(cx1-px))+((cy1-py)*(cy1-py)));
	  perim=(side1+side2+base)/2.0;
	  area=sqrt(perim*(perim-base)*(perim-side1)*(perim-side2));
	  height=area/(base*0.5);

	  /* if distance between particle and intersection line is small enough, **/
	  /**  particle makes an additional step and we check did it cross the intersection line or not**/
	  /**  if particle crossed intersection,  find the intersection point and move particle there**/
	 
	  if (height<=dist)
	    {
	      double  pr1, pr2,pr3,pr4, px1, py1,px2,py2;
 
	      PredictorStep();
 
	      px1=particle[np].position[0];
	      py1=particle[np].position[1];
	      px2=particle[np].prev_pos[0];
	      py2=particle[np].prev_pos[1];

	      pr1=(px1-cx1)*(py2-cy1)-(py1-cy1)*(px2-cx1);
	      pr2=(cx1-px1)*(cy2-py1)-(cy1-py1)*(cx2-px1);
	      pr3=(px1-cx2)*(py2-cy2)-(py1-cy2)*(px2-cx2);
	      pr4=(cx1-px2)*(cy2-py2)-(cy1-py2)*(cx2-px2);
   
	      if ((pr1*pr3<0)&&(pr2*pr4<0))
		{
		

		
		  /* particle crossed intersection */        
		  pr1=cx1*cy2-cy1*cx2;
		  pr2=px1*py2-py1*px2;
		  pr3=(cx1-cx2)*(py1-py2)-(cy1-cy2)*(px1-px2);
		  px=((px1-px2)*pr1-(cx1-cx2)*pr2)/pr3;
		  py=((py1-py2)*pr1-(cy1-cy2)*pr2)/pr3;
   
		  particle[np].position[0]=px; 
		  particle[np].position[1]=py;
                  
		  CheckNewCell();
		  intm=1;
               
		  if (particle[np].cell!=0)
		    {
		    if (no_out!=1)
                      {
                     if (node[int1-1].fracture[0]!=particle[np].fracture)
                        fract_p=node[int1-1].fracture[0];
                      else
                         fract_p=node[int1-1].fracture[1];   


			ParticleOutput(t, fract_p);
                        }
		      AcrossIntersection (prevcell, int1, int2);
		    }
		  else
		    {
		      //		    printf("Particle is lost on intersection. \n");
		    }
         
		}
	    }//if height
	} 
   
      /* particle is in a cell at the ending point of intersection */ 
      if ((int1!=0)&&(int2!=0)&&(int3==0))
	{
	
	  double nposx=0, nposy=0, dist_init=0, dist_fin=0, inout=0;
	  int coutf=0;
	  
	  nposx=(cx1+cx2)/2;
	  nposy=(cy1+cy2)/2;

	  dist_init=pow((particle[np].position[0]-nposx),2)+pow((particle[np].position[1]-nposy),2);
     
	  PredictorStep();

	  dist_fin=pow((particle[np].position[0]-nposx),2)+pow((particle[np].position[1]-nposy),2);
         
         
 
	  if (dist_init>=dist_fin)
	    {
	    
	      /* particle moves toward the intersection */    
	      
	      particle[np].position[0]=nposx; 
	      particle[np].position[1]=nposy;
	      
	      CheckNewCell();
	       
	      intm=1;
	      if (particle[np].cell!=0)
		{
		  prevcell=particle[np].cell;
		  
		   
		  if (no_out!=1)
                     { 
                      if (node[int1-1].fracture[0]!=particle[np].fracture)
                        fract_p=node[int1-1].fracture[0];
                      else
                         fract_p=node[int1-1].fracture[1];

		    ParticleOutput(t, fract_p);
                     }
		  AcrossIntersection (prevcell, int1, int2);
		   
		}
	      else
        
		printf("Particle is lost on end of intersection. \n");
	    }
	  else
	    {
	      ind_int2=0;
	      for (i=0; i<node[int1-1].numneighb; i++)
		{
		  if (node[int1-1].indnodes[i]==int2)
	            {
		      ind_int2=i;
		      break;
		    }
		}      
	      if (ind_int2==0)
		printf("ind=0 \n");
	      for (i=0; i<4; i++)
		{
		  if (node[int1-1].fracts[ind_int2][i]==particle[np].fracture)
		    {
		      int indc=0;
		      indc=node[int1-1].cells[ind_int2][i];
		      inout=InOutFlowCell(indc, int1, nposx,nposy);
		      if (inout>0)
			coutf++;

		    } 
		}
	    
	      if (coutf>1)
		{
	   
		  particle[np].position[0]=nposx; 
		  particle[np].position[1]=nposy;
	      
		  CheckNewCell();
	       
		  intm=1;
		  if (particle[np].cell!=0)
		    {
		      prevcell=particle[np].cell;
		      
		       
		      if (no_out!=1)
                       {
                       if (node[int1-1].fracture[0]!=particle[np].fracture)
                        fract_p=node[int1-1].fracture[0];
                      else
                         fract_p=node[int1-1].fracture[1];
			ParticleOutput(t, fract_p);
                        }
		      AcrossIntersection (prevcell, int1, int2);
		 
		    }
		  else
		    {
		      //		    printf("Particle is lost on end of intersection. \n");
		    }
		}
	    }
	
	}
    }
  if ((int1==0)||(int2==0))
    {
      /* there is no node of cell belongs to intersection */    
      //      printf("Check if %d cell is on intersection %d %d %d %d %d dn %d %d %d  %d. \n", particle[np].cell, node[dn1-1].typeN, node[dn2-1].typeN, node[dn3-1].typeN, int1, int2, dn1, dn2, dn3,  particle[np].intcell);
    }
  return intm;
}     // hpham: end of CheckDistance


//////////////////////////////////////////////////////////////////////////////
double InOutFlowCell(int indcell, int int1, double nposx, double nposy)
{
  /********* function defines is cell is in flow or out flow onto intersection***/
  double inoutf=0;
  int n1n,n2n,n3n,v1v,v2v,v3v;
  struct intcoef lambda;

  double prevpos0=particle[np].position[0], prevpos1=particle[np].position[1];
  int prevfract=particle[np].fracture, previouscell=particle[np].cell;
  double products=0, product=0;
  particle[np].position[0]=nposx;
  particle[np].position[1]=nposy;

  n1n=cell[indcell-1].node_ind[0];
  n2n=cell[indcell-1].node_ind[1];
  n3n=cell[indcell-1].node_ind[2];

  v1v=cell[indcell-1].veloc_ind[0];
  v2v=cell[indcell-1].veloc_ind[1];
  v3v=cell[indcell-1].veloc_ind[2];
	      
  int thirdnode=0;
  double tnx=0,tny=0;

  if ((node[n1n-1].typeN!=12)&&(node[n1n-1].typeN!=2))
    thirdnode=n1n;
  if ((node[n2n-1].typeN!=12)&&(node[n2n-1].typeN!=2))
    thirdnode=n2n;
  if ((node[n3n-1].typeN!=12)&&(node[n3n-1].typeN!=2))
    thirdnode=n3n;

  tnx=node[thirdnode-1].coord_xy[0];
  tny=node[thirdnode-1].coord_xy[1];

  double vintx=0, vinty=0, velocx=0, velocy=0;

  vintx=node[int1-1].coord_xy[XindexC(int1,indcell-1)];
  vinty=node[int1-1].coord_xy[YindexC(int1,indcell-1)];
		  
  ChangeFracture(indcell);
	      
  lambda=CalculateWeights(n1n, n2n, n3n);    
 
  velocx=lambda.weights[0]*node[n1n-1].velocity[v1v][0]+lambda.weights[1]*node[n2n-1].velocity[v2v][0]+lambda.weights[2]*node[n3n-1].velocity[v3v][0];
  velocy=lambda.weights[0]*node[n1n-1].velocity[v1v][1]+lambda.weights[1]*node[n2n-1].velocity[v2v][1]+lambda.weights[2]*node[n3n-1].velocity[v3v][1];  

  /* calculate vector's cross product to define outgoing and incoming flow cells */
  product=((particle[np].position[0]-tnx)*(particle[np].position[1]-vinty))-((particle[np].position[0]-vintx)*(particle[np].position[1]-tny));

  products=(velocx*(particle[np].position[1]-vinty))-((particle[np].position[0]-vintx)*velocy);
  inoutf=products*product;
	      
	      
  particle[np].position[0]=prevpos0;
  particle[np].position[1]=prevpos1;
  particle[np].cell=previouscell;
  particle[np].fracture=prevfract;
	      
  return inoutf;
}
//////////////////////////////////////////////////////////////////////////////
/***** Function moves particle through  intersection ********/
void AcrossIntersection (int prevcell, int int1, int int2)
{
  struct intcoef lambda;
  int  indj=-1, k,  indcell, cell_win=0, indk=0;
  int n1n,n2n,n3n,v1v,v2v,v3v;
  double speedsq[4]={0.0,0.0, 0.0, 0.0},  velocx, velocy;
  double  products[4]={0.0, 0.0, 0.0, 0.0};
 

  if ((int1!=0)&&(int2!=0))
    {
  
      /* defines indj - index of current cell in int1 node list */
      k=0;
      do
	{
          if (node[int1-1].indnodes[k]==int2)
	    indj=k; 
	  k++; 
	}
      while ((indj<0) && (k<node[int1-1].numneighb));            
           
      if (indj<0) 
	printf(" Current cell not found: NODES %d %d pw %f %f %f \n", int1, int2, particle[np].weight[0], particle[np].weight[1], particle[np].weight[2]); 


      /* the loop on 4 neighboring cells with common edge: int1 - int2 */
      for (k=0; k<4; k++)
	{
	  if (node[int1-1].cells[indj][k]!=0)
	    {
	      indcell=node[int1-1].cells[indj][k];

	      n1n=cell[indcell-1].node_ind[0];
	      n2n=cell[indcell-1].node_ind[1];
	      n3n=cell[indcell-1].node_ind[2];

	      v1v=cell[indcell-1].veloc_ind[0];
	      v2v=cell[indcell-1].veloc_ind[1];
	      v3v=cell[indcell-1].veloc_ind[2];

	      int thirdnode=0;
	      double tnx=0,tny=0;

	      if ((node[n1n-1].typeN!=12)&&(node[n1n-1].typeN!=2))
		thirdnode=n1n;
	      if ((node[n2n-1].typeN!=12)&&(node[n2n-1].typeN!=2))
		thirdnode=n2n;
	      if ((node[n3n-1].typeN!=12)&&(node[n3n-1].typeN!=2))
		thirdnode=n3n;

	      tnx=node[thirdnode-1].coord_xy[0];
	      tny=node[thirdnode-1].coord_xy[1];

	      double product, vintx=0, vinty=0;

	      vintx=node[int1-1].coord_xy[XindexC(int1,indcell-1)];
	      vinty=node[int1-1].coord_xy[YindexC(int1,indcell-1)];
 
	      if ((indcell)==prevcell)
                indk=k;
	      /** move to the intersecting  fracture and recalculate coordinations  ***/
	      ChangeFracture(indcell);

	      lambda=CalculateWeights(n1n, n2n, n3n);    
 
	      velocx=lambda.weights[0]*node[n1n-1].velocity[v1v][0]+lambda.weights[1]*node[n2n-1].velocity[v2v][0]+lambda.weights[2]*node[n3n-1].velocity[v3v][0];
	      velocy=lambda.weights[0]*node[n1n-1].velocity[v1v][1]+lambda.weights[1]*node[n2n-1].velocity[v2v][1]+lambda.weights[2]*node[n3n-1].velocity[v3v][1];  

	      /* calculate vector's cross product to define outgoing and incoming flow cells */
	      product=((particle[np].position[0]-tnx)*(particle[np].position[1]-vinty))-((particle[np].position[0]-vintx)*(particle[np].position[1]-tny));

	      products[k]=(velocx*(particle[np].position[1]-vinty))-((particle[np].position[0]-vintx)*velocy);
	      products[k]=products[k]*product;
	      speedsq[k]=velocx*velocx+velocy*velocy;
	    }    
        }//loop on k
 
      cell_win=RandomSampling(products, speedsq, indj, int1, indk);  
 
      ChangeFracture(cell_win);
      particle[np].intcell=4;
      particle[np].prev_pos[0]=particle[np].position[0];
      particle[np].prev_pos[1]=particle[np].position[1];

    }
  return;
}
////////////////////////////////////////////////////////////////////////////// 
int RandomSampling(double products[4], double speedsq[4], int indj, int int1, int indk)
{    
  /*********** Weighted Random Sampling ****************/ 
  int win_cell=0, k;
  int count=0, outc[4]={0,0,0,0};
  double random_number, totalspeed=0, minsp=0.0;

  /* find outgoing flow cells */  
  for (k=0; k<4; k++)
    {
      if (products[k]<0)
	{
	  outc[count]=k;
	  count++;
	  totalspeed=totalspeed+speedsq[k];
 
	}

    }
  /* if no outgoing flow cells found - move to cell with largest velocity magnitude 
     (this should not happen, it will mean we have a physical flow sink)*/
  if (count==0)
    {
    
      minsp=0.0;
      for (k=0; k<4; k++)
	if ((speedsq[k]>minsp)&&(k!=indk)) 
	  {
	    minsp=speedsq[k];
	    win_cell=node[int1-1].cells[indj][k];
	  }
    }
  /* if only one cell found */
  if (count==1)
    {
      win_cell=node[int1-1].cells[indj][outc[0]];
    }


  if (count==2)
    {
      random_number=drand48();
      if (random_number<=(speedsq[outc[0]]/totalspeed))
	//  if (random_number<0.5)
	win_cell=node[int1-1].cells[indj][outc[0]];
      else
	win_cell=node[int1-1].cells[indj][outc[1]];
    }


  if (count==3)
    {
      random_number=drand48();
      if (random_number>((speedsq[outc[0]]+speedsq[outc[1]])/totalspeed))
	//      if (random_number<0.3)
	win_cell=node[int1-1].cells[indj][outc[2]];
      else
	{
          if (random_number<=(speedsq[outc[0]]/totalspeed))
	    //         if (random_number<0.6)
            win_cell=node[int1-1].cells[indj][outc[0]];
	  else
            win_cell=node[int1-1].cells[indj][outc[1]];
	}
    }


  if (count==4)
    {
      random_number=drand48();
      if (random_number>((speedsq[outc[0]]+speedsq[outc[1]]+speedsq[outc[2]])/totalspeed))
	win_cell=node[int1-1].cells[indj][outc[3]];
      else
	{
          if (random_number>((speedsq[outc[0]]+speedsq[outc[1]])/totalspeed))
            win_cell=node[int1-1].cells[indj][outc[2]];
	  else
	    {
	      if (random_number<=(speedsq[outc[0]]/totalspeed))
		win_cell=node[int1-1].cells[indj][outc[0]];
	      else
		win_cell=node[int1-1].cells[indj][outc[1]];
	    }
	}
    }
  return win_cell;
}
//////////////////////////////////////////////////////////////////////////////
void Moving2Center (int nnp, int cellnumber)
{
  /*** function moves particle to the center of the same cell ******/
  double centx=0, centy=0, n1x=0, n2x=0, n3x=0, n1y=0, n2y=0,n3y=0;
  int n1=0, n2=0, n3=0; 

//  printf("nnp = %d | cellnumber = %d \n",nnp, cellnumber);     //hpham

  n1=cell[cellnumber-1].node_ind[0];
  n2=cell[cellnumber-1].node_ind[1];
  n3=cell[cellnumber-1].node_ind[2];

  n1x=node[n1-1].coord_xy[Xindex(n1,nnp)];
  n1y=node[n1-1].coord_xy[Yindex(n1,nnp)];

  n2x=node[n2-1].coord_xy[Xindex(n2,nnp)];
  n2y=node[n2-1].coord_xy[Yindex(n2,nnp)];

  n3x=node[n3-1].coord_xy[Xindex(n3,nnp)];
  n3y=node[n3-1].coord_xy[Yindex(n3,nnp)];

  centx=n1x+n2x+n3x;
  centy=n1y+n2y+n3y;

  particle[nnp].position[0]=centx/3;
  particle[nnp].position[1]=centy/3;

  int in;
  in=InsideCell (particle[nnp].cell);     //function returns 1 if particle is inside cell

  return;
}
/////////////////////////////////////////////////////////////////////////////


int Moving2NextCell (int stuck, int k)
{
  /* in case when particle got stuck in one cell, we move it to the center of neighboring cell */ 
 
  int current_index=0;
  int  j, nc=0, i=0; 
 
 
  if (stuck<node[k-1].numneighb)
    {
      i=stuck;
   
      do   
	{     
	  for (j=0; j<4; j++)
	    {
	      if (node[k-1].fracts[i][j]==particle[np].fracture)
		{
	      
	          if (node[k-1].cells[i][j]!=particle[np].cell)
		    {
		      nc=node[k-1].cells[i][j];
	      
		      current_index=i+1;
		      //	      printf("current_index %d i %d j %d \n", current_index, i,j);
		      break;
	            }
	        }
	       
	    }



	  if (nc!=0)
	    {
	      //	      printf("particle %d moved from cell %d to cell %d \n", np+1, particle[np].cell, nc); 
	      particle[np].cell=nc;
	      Moving2Center (np, nc);
	
	      break;
	    }
	  i++;
	}
      while ((nc==0)&&(i<node[k-1].numneighb));
    }
  else
    {
      //    printf(" 'stuck' is too big  \n");
    }
  if (nc==0)
    {
      //   printf("moving cell was not found part %d fract %d cell %d k %d\n", np+1, particle[np].fracture, particle[np].cell, k);
    }
 
  return current_index;

}
/////////////////////////////////////////////////////////////////////////////
double CalculateCurrentDT()
{
  /*** functions returns current time step delta t ***********/
  int dn1=0, dn2=0, dn3=0,  dv1=0, dv2=0, dv3=0; 
  double current_delta_t;


  dn1=cell[particle[np].cell-1].node_ind[0];
  dn2=cell[particle[np].cell-1].node_ind[1];
  dn3=cell[particle[np].cell-1].node_ind[2];


  dv1=cell[particle[np].cell-1].veloc_ind[0];
  dv2=cell[particle[np].cell-1].veloc_ind[1];
  dv3=cell[particle[np].cell-1].veloc_ind[2];

 
  current_delta_t=node[dn1-1].timestep[dv1]*particle[np].weight[0]+node[dn2-1].timestep[dv2]*particle[np].weight[1]+node[dn3-1].timestep[dv3]*particle[np].weight[2];

  return current_delta_t;
}
////////////////////////////////////////////////////////////////////////////
int Xindex(int nodenum, int nnp)
{
  /*** functions returns the correct index of node's coordination ***********/
  int xind=0;
  if (node[nodenum-1].fracture[0]==particle[nnp].fracture)
    xind=0;
  if (node[nodenum-1].fracture[1]==particle[nnp].fracture)
    xind=3;
  return xind;
} 
//////////////////////////////////////////////////////////////////////////////  
int Yindex(int nodenum, int nnp)
{
  /*** functions returns the correct index of node's coordination ***********/
  int yind=0;
  if (node[nodenum-1].fracture[0]==particle[nnp].fracture)
    yind=1;
  if (node[nodenum-1].fracture[1]==particle[nnp].fracture)
    yind=4;
  return yind;
}  
/////////////////////////////////////////////////////////////////////////////// 
void Moving2NextCellBound(int prevcell)
{
  /**** function moves particle to internal cell if it gets out of fracture****/
  int n1=0, n2=0, n3=0; 
  n1=cell[prevcell-1].node_ind[0];
  n2=cell[prevcell-1].node_ind[1];
  n3=cell[prevcell-1].node_ind[2];
  
  int k, j, nc=0, i=0, nc0=0, nc10=0, nc12=0, ncb=0; 
  k=0;
  
  if (node[n1-1].typeN==0)
    k=n1;
  else
    {  
      if (node[n2-1].typeN==0)
	k=n2; 
      else 
	k=n3;
    }  
 
  
  do   
    {     
      for (j=0; j<4; j++)
	{
	  if (node[k-1].fracts[i][j]==particle[np].fracture)
	    {
	      
	      if (node[k-1].cells[i][j]!=prevcell)
		{
		  nc=node[k-1].cells[i][j];
	      
	      
		  if (node[cell[nc-1].node_ind[0]-1].typeN+node[cell[nc-1].node_ind[1]-1].typeN+node[cell[nc-1].node_ind[2]-1].typeN==0)
		    {
		      nc0=nc;
		      break;
		    }
		  if (node[cell[nc-1].node_ind[0]-1].typeN+node[cell[nc-1].node_ind[1]-1].typeN+node[cell[nc-1].node_ind[2]-1].typeN==10)
		    nc10=nc;
		  if (node[cell[nc-1].node_ind[0]-1].typeN+node[cell[nc-1].node_ind[1]-1].typeN+node[cell[nc-1].node_ind[2]-1].typeN==12)
		    nc12=nc;
		  if (node[cell[nc-1].node_ind[0]-1].typeN+node[cell[nc-1].node_ind[1]-1].typeN+node[cell[nc-1].node_ind[2]-1].typeN >15)
		    ncb=nc;
		}
	    }
	       
	}
	     	
      if (nc0!=0)
	{
	  //	  printf("particle %d moved from boundary cell %d to cell %d \n", np+1, prevcell, nc0); 
	  particle[np].cell=nc0;
	  Moving2Center (np, nc0);
	
	  break;
	}
      i++;
    }
  while (i<node[k-1].numneighb);
  
  
  if ((nc0==0) && (nc10!=0))
    {
      //      printf("particle %d moved from boundary cell %d to cell %d \n", np+1, prevcell, nc10); 
      particle[np].cell=nc10;
      Moving2Center (np, nc10);
    }
  else
    {
      if ((nc0==0) && (nc10==0) && (nc12!=0))
  	{
	  //	  printf("particle %d moved from boundary cell %d to cell %d \n", np+1, prevcell, nc12); 
	  particle[np].cell=nc12;
	  Moving2Center (np, nc12);
   	}
      else
  	{
	  if ((nc0==0) && (nc10==0) && (nc12==0) && (ncb!=0))
	    {
	      //	      printf("particle %d moved from boundary cell %d to cell %d \n", np+1, prevcell, ncb); 
	      particle[np].cell=ncb;
	      Moving2Center (np, ncb);
	    }
	  else
	    {
	      if (nc0 +nc10+nc12+ncb ==0)
		{
		  //		printf("moving cell from bound was not found part %d fract %d cell %d k %d\n", np+1, particle[np].fracture, prevcell, k);
		}
	    }
   	}
    }
  return;
}
//////////////////////////////////////////////////////////////////////////////
void ParticleOutput (int currentt, int fract_p)
{
  /****** function outputs particles data int external file *******************/


  double posit[3]={0.0, 0.0, 0.0}, veloc[3]={0.0, 0.0, 0.0};
  double startx, starty, endx, endy, midx, midy, time, velocity_t=0.0, obeta=0.0, length_t=0.0;
  int i,tstart=-1, pcell=0, pfrac=0, tend=0, tmid; 
  int time_l, kdiv=2;
  double eps=0.05, pressure=0.0;
  struct posit3d particle3dp, particle3dv;
 
  startx= tempdata[0].position2d[0];
  starty= tempdata[0].position2d[1];
  posit[0]= tempdata[0].position3d[0];
  posit[1]= tempdata[0].position3d[1];
  posit[2]= tempdata[0].position3d[2];
  veloc[0]= tempdata[0].velocity3d[0];
  veloc[1]= tempdata[0].velocity3d[1];
  veloc[2]= tempdata[0].velocity3d[2];
  pcell= tempdata[0].cellp;
  pfrac= tempdata[0].fracturep;
  time= tempdata[0].timep;
  obeta= tempdata[0].betap;
  length_t=tempdata[0].length_t;
  pressure=tempdata[0].pressure;

  
 //             if (traj_o==1)
 //                   fprintf(wint,"%5.12E %5.12E  %5.12E   %5.12E  %5.12E %d %5.12E %5.12E\n", length_t, time, posit[0], posit[1], posit[2], pfrac, obeta, pressure);

  tstart=tempdata[0].times;

  if (tstart>=0)
    {
      if (curv_o==1)
	{
	  /* output according to trajectory's curvature */ 
	  endx=particle[np].position[0];
	  endy=particle[np].position[1];
	  tend=currentt;
 
	  if (tstart!=tend) 
	    {
      
	      if (traj_o==1)
		{
  
  
  
		  fprintf(wv,"%05d  %5.12E %5.12E %5.12E %5.12E %5.12E %5.12E %05d %05d %5.12E %5.12E %5.12E %d %5.12E \n", tstart, posit[0], posit[1], posit[2], veloc[0], veloc[1], veloc[2], pcell, pfrac, time,node[cell[pcell-1].node_ind[0]-1].aperture, obeta, 0, pressure);
	
		}
	      nodeID++;
	      if (avs_o==1)
		{
		  
		  fprintf(wpt,"%05d %5.12E %5.12E %5.12E \n",nodeID, posit[0], posit[1], posit[2]); 
		  velocity_t=sqrt(pow(veloc[0],2)+pow(veloc[1],2)+pow(veloc[2],2));
		  fprintf(wpt_att,"%010d  %06d  %5.12E  %5.12E  %5.12E %5.12E  %5.12E  %5.12E %5.12E\n", nodeID, pfrac, time, velocity_t, veloc[0], veloc[1], veloc[2], node[cell[pcell-1].node_ind[0]-1].aperture, pressure);
		}
	      int tstep, flag=0, isch=0; 
	      double angle_m;

	      do 
		{

		  time_l=tend-tstart;
		  kdiv=kdiv*2;
		  tstep=(int) (time_l/2.0);
		  
	  
		  for (isch=0; isch<time_l; isch++)
		    {
		      tmid= tempdata[isch].times;
		      midx= tempdata[isch].position2d[0];
		      midy= tempdata[isch].position2d[1];
		      posit[0]= tempdata[isch].position3d[0];
		      posit[1]= tempdata[isch].position3d[1];
		      posit[2]= tempdata[isch].position3d[2];
		      veloc[0]= tempdata[isch].velocity3d[0];
		      veloc[1]= tempdata[isch].velocity3d[1];
		      veloc[2]= tempdata[isch].velocity3d[2];
		      pcell= tempdata[isch].cellp;
		      pfrac= tempdata[isch].fracturep;
		      time= tempdata[isch].timep;
		      obeta= tempdata[isch].betap;
                      pressure=tempdata[isch].pressure;
                      length_t=tempdata[isch].length_t;
		      if (tmid==tstart+tstep)
			break;
 
		    }
		  

		  angle_m=DefineAngle(startx-midx, starty-midy, endx-midx, endy-midy);
 
		  if (((angle_m>-eps)&& (angle_m<eps)) || ((angle_m>pi-eps)&& (angle_m<pi+eps)))
		    {
 
		      flag=1;
		      break;
		    }
		  else
		    {
		      tend=tmid;
		      endx=midx;
		      endy=midy;
 
		    }
		  if (tend==tstart)
		    {
		      flag=1;
		      break;
		    }

		}

	      while (flag==0);

	      if (kdiv<10)
		kdiv=kdiv*4;
    
	      time_l=t-tstart;
	      if (kdiv!=0)
		tstep=(int)(time_l/kdiv);
      
	      if (kdiv>time_l)
		{
		  tstep=2;
		  kdiv=(int)time_l/2.0;
		}
	      if (kdiv==0)
		{
		  kdiv=2;
		  tstep=1;
		}
	      for(i=0; i<kdiv-1; i++)
		{
 
		  for (isch=0; isch<time_l; isch++)
		    
		    {

		      tmid= tempdata[isch].times;
		      midx= tempdata[isch].position2d[0];
		      midy= tempdata[isch].position2d[1];
		      posit[0]= tempdata[isch].position3d[0];
		      posit[1]= tempdata[isch].position3d[1];
		      posit[2]= tempdata[isch].position3d[2];
		      veloc[0]= tempdata[isch].velocity3d[0];
		      veloc[1]= tempdata[isch].velocity3d[1];
		      veloc[2]= tempdata[isch].velocity3d[2];
		      pcell= tempdata[isch].cellp;
		      pfrac= tempdata[isch].fracturep;
		      time= tempdata[isch].timep;
		      obeta= tempdata[isch].betap;
                      pressure=tempdata[isch].pressure;
			length_t=tempdata[isch].length_t;
		      if (tmid==tstart+tstep*(i+1))
			break;

		    }
		

		  if (traj_o==1)
		    {

		      fprintf(wv,"%05d  %5.12E %5.12E %5.12E %5.12E %5.12E %5.12E %05d %05d %5.12E %5.12E %5.12E %d %5.12E\n", tmid, posit[0], posit[1], posit[2], veloc[0], veloc[1], veloc[2], pcell, pfrac, time, node[cell[pcell-1].node_ind[0]-1].aperture, obeta, 0, pressure);

		    }	    
		  nodeID++;
		  if (avs_o==1)
		    {
		      
		      fprintf(wpt,"%05d %5.12E %5.12E %5.12E \n",nodeID, posit[0], posit[1], posit[2]); 
		      velocity_t=sqrt(pow(veloc[0],2)+pow(veloc[1],2)+pow(veloc[2],2));
		      fprintf(wpt_att,"%010d  %06d  %5.12E  %5.12E  %5.12E %5.12E  %5.12E  %5.12E  %5.12E\n", nodeID, pfrac, time, velocity_t, veloc[0], veloc[1], veloc[2], node[cell[pcell-1].node_ind[0]-1].aperture, pressure);

		    }
		}


	    }
    
	}
      else
	{
	  /* in case of every time step output */
	  if (curv_o!=1)
	    {

	      tstart=tempdata[0].times;


	      time_l=currentt-tstart;
	      for (i=0; i<time_l; i++)
		{

                  tstart= tempdata[i].times;
                  startx= tempdata[i].position2d[0];
                  starty= tempdata[i].position2d[1];
                  posit[0]= tempdata[i].position3d[0];
                  posit[1]= tempdata[i].position3d[1];
                  posit[2]= tempdata[i].position3d[2];
                  veloc[0]= tempdata[i].velocity3d[0];
                  veloc[1]= tempdata[i].velocity3d[1];
                  veloc[2]= tempdata[i].velocity3d[2];
                  pcell= tempdata[i].cellp;
                  pfrac= tempdata[i].fracturep;
                  time= tempdata[i].timep;
                  obeta= tempdata[i].betap;
                  pressure=tempdata[i].pressure;
                  length_t=tempdata[i].length_t;
		  nodeID++;
		  if (avs_o==1)
		    {
		      
		      fprintf(wpt,"%05d %5.12E %5.12E %5.12E \n",nodeID, posit[0], posit[1], posit[2]); 
		      velocity_t=sqrt(pow(veloc[0],2)+pow(veloc[1],2)+pow(veloc[2],2));
		      fprintf(wpt_att,"%010d  %06d  %5.12E  %5.12E  %5.12E %5.12E  %5.12E  %5.12E %5.12E\n", nodeID, pfrac, time, velocity_t, veloc[0], veloc[1], veloc[2], node[cell[pcell-1].node_ind[0]-1].aperture, pressure);

		    }
		  if (traj_o==1)
		    fprintf(wv,"%05d  %5.12E %5.12E %5.12E %5.12E %5.12E %5.12E %05d %05d %5.12E %5.12E %5.12E %d %5.12E\n", tstart, posit[0], posit[1], posit[2], veloc[0], veloc[1], veloc[2], pcell, pfrac, time, node[cell[pcell-1].node_ind[0]-1].aperture, obeta, 0, pressure);
		}

	    } 
	}
      particle3dp=CalculatePosition3D();
      particle3dv=CalculateVelocity3D();
  
      if (traj_o==1)
	{
       fprintf(wint,"%5.12E %5.12E  %5.12E   %5.12E  %5.12E %d %5.12E %5.12E\n", length_t,particle[np].time, particle3dp.cord3[0], particle3dp.cord3[1],particle3dp.cord3[2], particle[np].fracture, obeta, particle[np].pressure);

	  fprintf(wv,"%05d  %5.12E %5.12E %5.12E %5.12E %5.12E %5.12E %05d %05d %5.12E %5.12E %5.12E %d %5.12E\n", t, particle3dp.cord3[0], particle3dp.cord3[1],particle3dp.cord3[2],particle3dv.cord3[0], particle3dv.cord3[1],particle3dv.cord3[2],particle[np].cell, particle[np].fracture, particle[np].time, node[cell[pcell-1].node_ind[0]-1].aperture, obeta, fract_p, particle[np].pressure);
    
	}
      nodeID++;
      if (avs_o==1)
	{
	  
	  fprintf(wpt,"%05d %5.12E %5.12E %5.12E \n",nodeID, particle3dp.cord3[0], particle3dp.cord3[1],particle3dp.cord3[2]); 
	  
	  velocity_t=sqrt(pow(particle3dv.cord3[0],2)+pow(particle3dv.cord3[1],2)+pow(particle3dv.cord3[2],2));
	  fprintf(wpt_att,"%010d  %06d  %5.12E  %5.12E  %5.12E %5.12E  %5.12E  %5.12E %5.12E\n", nodeID, particle[np].fracture, particle[np].time, velocity_t, particle3dv.cord3[0], particle3dv.cord3[1],particle3dv.cord3[2], node[cell[pcell-1].node_ind[0]-1].aperture, particle[np].pressure);
	 
	}
    }     // hpham: end of if (tstart>=0) line 2247
  
  timecounter=0;
  free(tempdata);

  return;
}
/////////////////////////////////////////////////////////////////////////////
void FinalPosition()
{
  /* function calculates particles final position on out-flow surface */
  int n1, n2, n3;
  int n1out=0, n2out=0;
  n1=cell[particle[np].cell-1].node_ind[0];
  n2=cell[particle[np].cell-1].node_ind[1];
  n3=cell[particle[np].cell-1].node_ind[2];
 
  double cx1,cx2,cy1,cy2, px1,px2,py1,py2, ap,bp,cp,as,bs,cs, deter,xint,yint;
 
      
           
  if ((node[n1-1].typeN==210)||(node[n1-1].typeN==212)||(node[n1-1].typeN==200)||(node[n1-1].typeN==202))
    n1out=n1;
	    
  if ((node[n2-1].typeN==210)||(node[n2-1].typeN==212)||(node[n2-1].typeN==200)||(node[n2-1].typeN==202))
    {
      if (n1out==0)
	n2out=n2;
      else
	n1out=n2;
    }
  if ((node[n3-1].typeN==210)||(node[n3-1].typeN==212)||(node[n3-1].typeN==200)||(node[n3-1].typeN==202))
    {
      if (n2out==0)
	n2out=n3;
    }
	    
  if ((n1out!=0) && (n2out!=0))
    {
      //   two vertexes of particles cell are on out flow boundary
      double cx1,cx2,cy1,cy2, px1,px2,py1,py2, ap,bp,cp,as,bs,cs, deter,xint,yint; 
      cx1=node[n1out-1].coord_xy[Xindex(n1out,np)];
      cy1=node[n1out-1].coord_xy[Yindex(n1out,np)];

      cx2=node[n2out-1].coord_xy[Xindex(n2out,np)];
      cy2=node[n2out-1].coord_xy[Yindex(n2out,np)];
 
      as=cy2-cy1;
      bs=cx1-cx2;
      cs=as*cx1+bs*cy1;
      px1=particle[np].position[0];
      py1=particle[np].position[1];
      
      px2=particle[np].prev_pos[0];
      py2=particle[np].prev_pos[1];
      
      ap=py2-py1;
      bp=px1-px2;
      cp=ap*px1+bp*py1;
      deter=ap*bs-as*bp;
      xint=(bs*cp-bp*cs)/deter;
      yint=(ap*cs-as*cp)/deter;
      
      double distance=0, tfinal;
      distance=sqrt((xint-px1)*(xint-px1)+(yint-py1)*(yint-py1));
      tfinal=distance/sqrt(particle[np].velocity[0]*particle[np].velocity[0]+particle[np].velocity[1]*particle[np].velocity[1]);
     
      particle[np].prev_pos[0]=particle[np].position[0];
      particle[np].prev_pos[1]=particle[np].position[1];
      particle[np].position[0]=xint;
      particle[np].position[1]=yint;
      particle[np].time=tfinal+particle[np].time;
      
      
      //     printf("x %12.5e y %12.5e z %12.5e %d time %5.12e\n", particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2], particle[np].cell, particle[np].time);
 
    } 
  else
    {
      int ncent=0, nnext=0;
      ncent=n1out+n2out;
	  
      if (ncent!=0)
	{
	    
	  int ii=0;
	  for (ii=0; ii<node[ncent-1].numneighb; ii++)
	    {
		                                                                                                                         
	      if ((node[node[ncent-1].indnodes[ii]-1].typeN==212)||(node[node[ncent-1].indnodes[ii]-1].typeN==210)||(node[node[ncent-1].indnodes[ii]-1].typeN==202))
		{
		   
		  nnext=node[ncent-1].indnodes[ii];
		  int newcell=0;
		  newcell=node[ncent-1].cells[ii][0]; 
		    	       
		  cx1=node[ncent-1].coord_xy[Xindex(ncent,np)];
		  cy1=node[ncent-1].coord_xy[Yindex(ncent,np)];

		  cx2=node[nnext-1].coord_xy[Xindex(nnext,np)];
		  cy2=node[nnext-1].coord_xy[Yindex(nnext,np)];
 
     
		  as=cy2-cy1;
		  bs=cx1-cx2;
		  cs=as*cx1+bs*cy1;
		  px1=particle[np].position[0];
		  py1=particle[np].position[1];
      
		  px2=particle[np].prev_pos[0];
		  py2=particle[np].prev_pos[1];
      
		  ap=py2-py1;
		  bp=px1-px2;
		  cp=ap*px1+bp*py1;
		  deter=ap*bs-as*bp;
		  xint=(bs*cp-bp*cs)/deter;
		  yint=(ap*cs-as*cp)/deter;
      
		  double distance=0, tfinal;
		  distance=sqrt((xint-px1)*(xint-px1)+(yint-py1)*(yint-py1));
		  tfinal=distance/sqrt(particle[np].velocity[0]*particle[np].velocity[0]+particle[np].velocity[1]*particle[np].velocity[1]);
      
		  particle[np].prev_pos[0]=particle[np].position[0];
		  particle[np].prev_pos[1]=particle[np].position[1];
		  particle[np].position[0]=xint;
		  particle[np].position[1]=yint;
		  particle[np].time=tfinal+particle[np].time;
		  particle[np].cell=newcell;
           
		  //      printf("INSIDE! x %12.5e y %12.5e z %12.5e %d time %5.12e\n", particle3dposit.cord3[0], particle3dposit.cord3[1],particle3dposit.cord3[2], newcell, particle[np].time);
		    
		     
		}
	      if (nnext!=0)
		break;
	    }
	      
	    
	    
	}
      //	    else
      //	  printf(" %d not in the cell, %d\n", np, particle[np].cell);  
	    

    }
    
    

  return;
}
///////////////////////////////////////////////////////////////////////////
struct lagrangian CalculateLagrangian(double xcurrent, double ycurrent, double zcurrent, double xprev, double yprev, double zprev)
{
  /**** calculates Lagrangian variables tau and beta;
	we take one triangular cell as a segment;
  *****/ 
  struct lagrangian lagvariable;
  struct posit3d particle3dv;
  double currentdistance=0.0, deltatau=0.0, deltabeta=0.0, velsquare=0.0;
  particle3dv=CalculateVelocity3D();
  currentdistance=pow((xcurrent-xprev),2)+pow((ycurrent-yprev),2)+pow((zcurrent-zprev),2);
 
  velsquare=(pow(particle3dv.cord3[0],2)+pow(particle3dv.cord3[1],2)+pow(particle3dv.cord3[2],2));
 
  deltatau=currentdistance/velsquare;
 
 
  if (currentdistance>0.0)
    {
      deltabeta=sqrt(currentdistance)/(sqrt(velsquare)*(node[cell[particle[np].cell-1].node_ind[0]-1].aperture*0.5));
      lagvariable.tau=sqrt(deltatau);
    }
  else
    {
      deltabeta=0.0;
      lagvariable.tau=0.0;
    }
  lagvariable.betta=deltabeta;
  lagvariable.initx=xcurrent;
  lagvariable.inity=ycurrent;
  lagvariable.initz=zcurrent;

   

  return lagvariable;
}
