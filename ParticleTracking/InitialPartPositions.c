// Copyright 2019, Nevada System of Higher Education on Behalf of the Desert Research Institute, and 
// Copyright 2018, Triad National Secuirty, LLC. All rights reserved.

// This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser 
// General Public License as published by the Free Software Foundation; either version 3.0 of the License, 
// or (at your option) any later version.  This software requires you separately obtain dfnWorks under an 
// appropriate license from the Los Alamos National Laboratory (LANL), which is operated by 
// Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.

#include <stdio.h> 
#include <search.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "FuncDef.h" 
#include <unistd.h>
#include <time.h>

int nnode_iwel;
unsigned int *nodeinj_wel;

struct inpfile {
  char filename[120];
  long int flag;
  double param;
}; 
  

int InitPos2()     //hpham marked this line for tracking purpose.
{
  /***** Function defines initial positions of particles *************/
  /***** Locates parts_fracture number of particles on edge of fracture on flow_in zone
     with equal distance between each other *********************/
  int i,ii, j, k_current, k_new=0, numbf=1,  frc,firstn,lastn, flag_in=0,parts_fracture;
  struct inpfile initfile;
  double parts_dist=0;
  int  zonenumb_in=0, first_ind=0, last_ind=0;
  double ixmin=0, ixmax=0, iymin=0, iymax=0, izmin=0, izmax=0; 
  double ixmin2=0, ixmax2=0, iymin2=0, iymax2=0, izmin2=0, izmax2=0; 
  double px[4]={0.0, 0.0, 0.0, 0.0}, py[4]={0.0, 0.0, 0.0, 0.0};
  double hp_sum_flux1=0;
  FILE *fid_hplog = OpenFile ("hpham.log","w");
  fprintf(fid_hplog,"i, flux, Node_ID, npar, Frac_ID, Cell_ID, x, y, z \n");

  /* [hpham] Calculate number of fractures in the zone where initial particles are released ****/
  int numbf2=1;     // numbf2 is the number of fractures in the injection well zone 
  double hp_flux[nnode_iwel], hp_flux_w[nnode_iwel];     //flux at node
  double hp_total_flux=0;
  
  printf("\n [IPP.c] Number of nodes in the injection well zone: nnode_iwel = %d \n", nnode_iwel);      //hpham

  for (i=0; i<nnode_iwel; i++)      //hpham's comments: nnode_iwel - number of nodes in flow-in zone from boundary.zone
    {

      // sum of flux at current node. 
      for (j=0; j<node[nodeinj_wel[i]-1].numneighb; j++)
      {
//        hp_sum_flux1=hp_sum_flux1+fabs(node[nodeinj_wel[i]-1].flux[j]);
        hp_sum_flux1=hp_sum_flux1+(node[nodeinj_wel[i]-1].flux[j]);
      }
      hp_flux[i]=hp_sum_flux1;     //total flux at current node.
//      fprintf(fid_hplog," [IPP.c] i = %3d, flux = %7.4E, Node_ID = %6d, Frac_ID = %6d, Cell_ID = %6d, x = %8.4f, y = %8.4f, z = %8.4f \n",i+1, hp_flux[i], nodeinj_wel[i], node[nodeinj_wel[i]-1].fracture[0], node[nodeinj_wel[i]-1].cells, node[nodeinj_wel[i]-1].coord[0], node[nodeinj_wel[i]-1].coord[1], node[nodeinj_wel[i]-1].coord[2]);
//      printf(" [IPP.c] i = %3d, flux = %7.4E \n",i+1,hp_flux[i]);
//      printf(" [IPP.c] i = %3d, \n",i+1);

/*       if ((i>0)&&(node[nodeinj_wel[i]-1].fracture[0]!=node[nodeinj_wel[i-1]-1].fracture[0])) 
          numbf2++;  */
	  
    }     // end of for with i
      
  // Calculate total flux of all nodes
  for (i=0; i<nnode_iwel; i++)     
    {
      hp_total_flux=hp_total_flux+hp_flux[i];    
    }
//    printf("\n [IPP.c] total flux = %7.4E (m3/s) \n", hp_total_flux);     //hpham_comment

  //Calculate flux weights for each node.
  for (i=0; i<nnode_iwel; i++) 
    {
      hp_flux_w[i]=hp_flux[i]/hp_total_flux;     
//      printf("[IPP.c] node = %4d, flux = %8.4E, flux_weight = %8.4E \n",i+1, hp_flux[i], hp_flux_w[i]);     //hpham_comment
    }


/*   if ((numbf2==0)&&(nnode_iwel!=0))
      numbf2=1;
  printf("\n [IPP.c] %d fracture(s) in in-flow boundary zone \n",numbf2); */

      
  int res;
  initfile=Control_File("init_oneregion2:",16);
  res=strncmp(initfile.filename,"yes",3);
  if (res==0)
    {
      flag_in=6;      // Read inputs
      flag_w=1;
      /* 6th option: user specifies a region and all particles start */
      /* from a region inside a model domain                         */
  
      initfile = Control_Data("in_partn2:",10 );
      npart=initfile.flag;
//      /*  memory allocatation for particle structures */
//      particle=(struct contam*) malloc ((npart*2)*sizeof(struct contam));
	  printf("\n");
	  printf("\n [hpham] %d particles are injected in ER 20-1 #1 \n", npart);
      printf("\n");

      initfile = Control_Param("in_xmin2:",9 );
      ixmin2=initfile.param;    
      initfile = Control_Param("in_xmax2:",9 );
      ixmax2=initfile.param; 
      initfile = Control_Param("in_ymin2:",9 );
      iymin2=initfile.param; 
      initfile = Control_Param("in_ymax2:",9 );
      iymax2=initfile.param; 
      initfile = Control_Param("in_zmin2:",9 );
      izmin2=initfile.param; 
      initfile = Control_Param("in_zmax2:",9 );
      izmax2=initfile.param; 
  
//          initfile = Control_Data("in_region_loc:",14 );
//          zonenumb_in=initfile.flag;

//      printf(" xmin = %7.3f \n xmax = %7.3f \n ymin = %7.3f \n ymax = %7.3f \n zmin = %7.3f \n zmax = %7.3f \n", ixmin2, ixmax2, iymin2, iymax2, izmin2, izmax2 );

      /* define coordinations of region according to in-flow zone */
        }
      
//         printf(" x_well = %7.3f \n y_well = %7.3f \n well_radius = %7.3f \n top_well = %7.3f \n bot_well = %7.3f \n ", x0w, y0w, rwell,top_well, bot_well, opt_mwell_on );
      

      
      
  // Distribute the number of particles for each node (based on flux weight)      
  int hp_npar;     // the total number of particles placed at the injection well location.
  int hp_npar4node[nnode_iwel], hp_count=0;

  hp_npar=npart;

//  printf("\n The number of particles given in the input file is %d \n",hp_npar);      //hpham_comment

  particle=(struct contam*) malloc ((hp_npar*2)*sizeof(struct contam));
  for (i=0; i<nnode_iwel; i++)      
    {    
      hp_npar4node[i]=round(hp_npar*hp_flux_w[i]);
//      fprintf(fid_hplog," [IPP.c] i=%d, Node_ID %8d has %4d particles \n",i+1,nodeinj_wel[i], hp_npar4node[i]);
//      fprintf(fid_hplog," [IPP.c] i = %3d, flux = %7.4E, Node_ID = %6d, npar = %4d, Frac_ID = %6d, Cell_ID = %6d, x = %8.4f, y = %8.4f, z = %8.4f \n",i+1, hp_flux[i], nodeinj_wel[i], hp_npar4node[i],node[nodeinj_wel[i]-1].fracture[0], node[nodeinj_wel[i]-1].cells, node[nodeinj_wel[i]-1].coord[0], node[nodeinj_wel[i]-1].coord[1], node[nodeinj_wel[i]-1].coord[2]);  
    }
  
  // Recheck the number of particle (after rounding off)
  int hp_chk_npar=0;
	for (i=0; i<nnode_iwel; i++)      
    {    
      if (hp_npar4node[i]<0)
        hp_npar4node[i]=0;
      hp_chk_npar=hp_chk_npar+hp_npar4node[i];
//      printf("%10d %10d \n",i,hp_npar4node[i]);
      fprintf(fid_hplog,"%6d, %10.4E, %8d, %8d, %10d, %10d, %8.4f, %8.4f, %8.4f \n",i, hp_flux[i], nodeinj_wel[i], hp_npar4node[i],node[nodeinj_wel[i]-1].fracture[0], node[nodeinj_wel[i]-1].cells, node[nodeinj_wel[i]-1].coord[0], node[nodeinj_wel[i]-1].coord[1], node[nodeinj_wel[i]-1].coord[2]);  
//      printf("%6d, %10.4E, %8d, %8d, %10d, %10d, %8.4f, %8.4f, %8.4f \n",i, hp_flux[i], nodeinj_wel[i], hp_npar4node[i],node[nodeinj_wel[i]-1].fracture[0], node[nodeinj_wel[i]-1].cells, node[nodeinj_wel[i]-1].coord[0], node[nodeinj_wel[i]-1].coord[1], node[nodeinj_wel[i]-1].coord[2]);  
    }
  
  for (i=0; i<nnode_iwel; i++)      //hpham's comments: numbf2: the number of fractures
    {
      for (j=0; j<hp_npar4node[i]; j++)
	  {
        //particle[i].position[0]=node[nodeinj_wel[i]-1].coord[0];   
        //particle[i].position[1]=node[nodeinj_wel[i]-1].coord[1];
//        particle[hp_count].position[0]=node[nodeinj_wel[i]-1].coord_xy[0]+(drand48()-drand48())/10;     //hpham: place at the well location +- err_norm (0,0.5) m
        particle[hp_count].position[0]=node[nodeinj_wel[i]-1].coord_xy[0];
//		printf(" x = %8.4f, std_x = %12.4E \n", node[nodeinj_wel[i]-1].coord_xy[0], drand48()/20);
//        particle[hp_count].position[1]=node[nodeinj_wel[i]-1].coord_xy[1]+(drand48()-drand48())/10;     
        particle[hp_count].position[1]=node[nodeinj_wel[i]-1].coord_xy[1];
        particle[hp_count].velocity[0]=0.;
        particle[hp_count].velocity[1]=0.;

        particle[hp_count].fracture=node[nodeinj_wel[i]-1].fracture[0];
//        particle[hp_count].fracture=1;


        //printf("[InitialPartiPositions.c: ] node = %4d, fractures= %6d, %6d \n",nodeinj_wel[i]-1, node[nodeinj_wel[i]-1].fracture[0], node[nodeinj_wel[i]-1].fracture[1]);     //hpham: 
        //printf("[InitialPartiPositions.c: ] particle[i].fracture = %6d \n", particle[i].fracture);     //hpham:
        //printf("[IPP.c: ] node.coord_xy = %8.4f, %8.4f \n", node[nodeinj_wel[i]-1].coord_xy[0], node[nodeinj_wel[i]-1].coord_xy[1]);     //hpham:
        particle[hp_count].intcell=0;
        particle[hp_count].time=0.0; 
        particle[hp_count].fl_weight=0.0;
        particle[hp_count].cell=node[nodeinj_wel[i]-1].cells;    // 
		hp_count=hp_count+1;
	  }
    }
  fclose(fid_hplog);    //hpham
  k_new=hp_chk_npar;
  return k_new;
}

     

// Hello, testing Sublime sync here


   





  

int InitPos()     //hpham marked this line for tracking purpose.
{
  /***** Function defines initial positions of particles *************/
  /***** Locates parts_fracture number of particles on edge of fracture on flow_in zone
     with equal distance between each other *********************/
  int i, k_current, k_new=0, numbf=1,  frc,firstn,lastn, flag_in=0,parts_fracture;
  struct inpfile initfile;
  double parts_dist=0;
  int  zonenumb_in=0, first_ind=0, last_ind=0;
  double ixmin=0, ixmax=0, iymin=0, iymax=0, izmin=0, izmax=0; 
//  double ixmin2=0, ixmax2=0, iymin2=0, iymax2=0, izmin2=0, izmax2=0; 
  double xcyl_top2=0, ycyl_top2=0, zcyl_top2=5, xcyl_bot2=0, ycyl_bot2=0, zcyl_bot2=-5, rad_cyli2=0.3, x_init_par_cyl=0, y_init_par_cyl=0, z_init_par_cyl=0; 
  double px[4]={0.0, 0.0, 0.0, 0.0}, py[4]={0.0, 0.0, 0.0, 0.0};
//  FILE *fid_initpos_cyl = OpenFile ("ParticleInitCoordR_hpham.dat","w");
   
  /* calculate number of fractures in in-flow boundary face of domain ****/
  // printf("fract %d ", node[nodezonein[0]-1].fracture[0]);
  for (i=0; i<nzone_in; i++)     //hpham's comments: nzone_in - number of nodes in flow-in zone from boundary.zone
    {
      if ((i>0)&&(node[nodezonein[i]-1].fracture[0]!=node[nodezonein[i-1]-1].fracture[0]))
        {
      //    printf("fract %d ", node[nodezonein[i]-1].fracture[0]);
            numbf++;
        }
    }
 
  if ((numbf==0)&&(nzone_in!=0))
    numbf=1;
  printf(" %d fracture(s) in in-flow boundary zone \n",numbf); 
  
 
  
  
  double inter_p[numbf][4];
  int inter_fr[numbf];   
  int res;
  initfile=Control_File("init_nf:",8);
  res=strncmp(initfile.filename,"yes",3);
  printf(" res at line 237 = %d \n",res); 

  /* first option: the same number of particles on every boundary edge */
  if (res==0)
    {
      
      flag_in=1;
      flag_w=1; 
      initfile = Control_Data("init_partn:",11 );
      parts_fracture=initfile.flag;
      printf("\n  %d  particles per boundary fracture edge \n", parts_fracture);
      npart=(numbf+1)*parts_fracture;
       
 
      /*  memory allocatation for particle structures */
      particle=(struct contam*) malloc (npart*sizeof(struct contam));
    }
  else
    {     
      initfile=Control_File("init_eqd:",9);
      res=strncmp(initfile.filename,"yes",3);
      if (res==0)
    {
      flag_in=2;
      flag_w=1;
      /* second option: calculate total length of boundary edges */
      /* define the distance between particles and place particles */
      /* equidistant from each other on all edges */
     
      initfile = Control_Data("init_npart:",11 );
      parts_fracture=initfile.flag;
      npart=(numbf)*parts_fracture*2;
      /*  memory allocatation for particle structures */
      particle=(struct contam*) malloc ((npart)*sizeof(struct contam));
      // define a distance between particles
     
      frc=node[nodezonein[0]-1].fracture[0];
      firstn=nodezonein[0];
      lastn=nodezonein[0];
      double length2=0.0, t_length=0.0;
      for (i=0; i<nzone_in; i++)
        {
          if (node[nodezonein[i]-1].fracture[0]!=frc)
        {
          lastn=nodezonein[i-1];
          length2=pow((node[firstn-1].coord[0]-node[lastn-1].coord[0]),2)+pow((node[firstn-1].coord[1]-node[lastn-1].coord[1]),2)+pow((node[firstn-1].coord[2]-node[lastn-1].coord[2]),2);
          t_length=t_length+sqrt(length2);
 
          frc=node[nodezonein[i]-1].fracture[0];
          firstn=nodezonein[i];
        }
        }
      lastn=nodezonein[nzone_in-1];
      length2=pow((node[firstn-1].coord[0]-node[lastn-1].coord[0]),2)+pow((node[firstn-1].coord[1]-node[lastn-1].coord[1]),2)+pow((node[firstn-1].coord[2]-node[lastn-1].coord[2]),2);
      t_length=t_length+sqrt(length2);
 
          parts_dist=t_length/(parts_fracture*numbf);
      printf("\n  Particles placed on %f [m]  from each other  \n", parts_dist);
    }  
      
	  
	  else
      {
      initfile=Control_File("init_oneregion:",15);
      res=strncmp(initfile.filename,"yes",3);
      if (res==0)
        {
          flag_in=3;      // Read inputs
          flag_w=1;
          /* third option: user specifies a region and all particles start */
          /* from the edges that fit inside the region */
      
          initfile = Control_Data("in_partn:",9 );
          npart=initfile.flag;
          /*  memory allocatation for particle structures */
          particle=(struct contam*) malloc ((npart*2)*sizeof(struct contam));
          printf("\n Initially particles have the same starting region \n");
      
          initfile = Control_Param("in_xmin:",8 );
          ixmin=initfile.param;    
          initfile = Control_Param("in_xmax:",8 );
          ixmax=initfile.param; 
          initfile = Control_Param("in_ymin:",8 );
          iymin=initfile.param; 
          initfile = Control_Param("in_ymax:",8 );
          iymax=initfile.param; 
          initfile = Control_Param("in_zmin:",8 );
          izmin=initfile.param; 
          initfile = Control_Param("in_zmax:",8 );
          izmax=initfile.param; 
      
          initfile = Control_Data("in-flow-boundary:",17 );
          zonenumb_in=initfile.flag;
          /* define coordinations of region according to in-flow zone */
          if ((zonenumb_in==1)||(zonenumb_in==2)) 
        { 
          px[0]=ixmin;
          py[0]=iymin;
          px[1]=ixmin;
          py[1]=iymax;
          px[2]=ixmax;
          py[2]=iymax;
          px[3]=ixmax;
          py[3]=iymin;
        }
       
          if ((zonenumb_in==3)||(zonenumb_in==5)) 
        { 
          px[0]=izmin;
          py[0]=iymin;
          px[1]=izmin;
          py[1]=iymax;
          px[2]=izmax;
          py[2]=iymax;
          px[3]=izmax;
          py[3]=iymin;
        }
       
          if ((zonenumb_in==4)||(zonenumb_in==6)) 
        { 
          px[0]=ixmin;
          py[0]=izmin;
          px[1]=ixmin;
          py[1]=izmax;
          px[2]=ixmax;
          py[2]=izmax;
          px[3]=ixmax;
          py[3]=izmin;
        }
        }
        
        
        else



      {     // hpham added these lines to read PTDFN_control.dat opt6 (region inside domain)
      int count_inj=1;
	  initfile=Control_File("init_oneregion2:",16);
      res=strncmp(initfile.filename,"yes",3);
      if (res==0)
        {
          flag_in=6;      // Read inputs
          flag_w=1;
          FILE *fid_initpos_cyl = OpenFile ("ParticleInitCoordR.dat","w");
          /* 6th option: user specifies a region and all particles start */
          /* from a region inside a model domain                         */
      
          initfile = Control_Data("in_partn2:",10 );
          npart=initfile.flag;
          /*  memory allocatation for particle structures */
          particle=(struct contam*) malloc ((npart*2)*sizeof(struct contam));

		      printf("\n hpham: Initially particles are uniformly distributed in a vertical cylinder with: \n");
      
          initfile = Control_Param("xcyl_top:",9 );
          xcyl_top2=initfile.param;    
          initfile = Control_Param("ycyl_top:",9 );
          ycyl_top2=initfile.param; 
          initfile = Control_Param("zcyl_top:",9 );
          zcyl_top2=initfile.param; 
          
          initfile = Control_Param("xcyl_bot:",9 );
          xcyl_bot2=initfile.param; 
          initfile = Control_Param("ycyl_bot:",9 );
          ycyl_bot2=initfile.param; 
          initfile = Control_Param("zcyl_bot:",9 );
          zcyl_bot2=initfile.param; 

          initfile = Control_Param("rad_cyli:",9 );
          rad_cyli2=initfile.param; 
		  
//          initfile = Control_Data("in_region_loc:",14 );
//          zonenumb_in=initfile.flag;

          printf(" xtop = %7.3f \n ytop = %7.3f \n ztop = %7.3f \n xbot = %7.3f \n ybot = %7.3f \n zbot = %7.3f \n", xcyl_top2, ycyl_top2, zcyl_top2, xcyl_bot2, ycyl_bot2, zcyl_bot2);

//       Generate file ... ParticleInitCoordR.dat ----------------------------		  
		  
          fprintf(fid_initpos_cyl,"npar  %6d \n", npart);
		 if (npart < 1000)
              printf("\n Error: Need at least %d particles \n", npart);
		  for (i=1; i < npart+1; i++) 
			{
				x_init_par_cyl=xcyl_top2+(drand48()-drand48())*rad_cyli2/2;     //hpham: place at (xcyl_top2, xyyl_top2) +- err_norm (0,0.5) m
				y_init_par_cyl=ycyl_top2+(drand48()-drand48())*rad_cyli2/2;     //hpham: place at (xcyl_top2, xyyl_top2) +- err_norm (0,0.5) m
				z_init_par_cyl=zcyl_top2 - (count_inj-1)*(zcyl_top2-zcyl_bot2)/1000;     //hpham: place at (xcyl_top2, xyyl_top2) +- err_norm (0,0.5) m
				
				if ( i % (npart/1000) == 0 )
				{
					count_inj = count_inj + 1;
				}
				fprintf(fid_initpos_cyl,"%8.4f %8.4f %8.4f \n", x_init_par_cyl, y_init_par_cyl, z_init_par_cyl);  
			}
			fclose(fid_initpos_cyl);    //hpham
		    printf("\n Generated file  ParticleInitCoordR.dat \n"); 


            FILE *fid_lgi = OpenFile ("definedist.lgi","w");
			fprintf(fid_lgi,"read / avs / full_mesh.inp / mo2 \n");
			fprintf(fid_lgi, "\n");
			fprintf(fid_lgi, "cmo / create/ mo1  \n");
			fprintf(fid_lgi, "cmo/readatt/ mo1/ npar / 1,0,0 / ParticleInitCoordR.dat  \n");
			fprintf(fid_lgi, "cmo/readatt/ mo1/ xic, yic, zic/ 1,0,0 / ParticleInitCoordR.dat  \n");
			fprintf(fid_lgi, "cmo / set_id/mo2 / node/ n_num \n");
			fprintf(fid_lgi, "cmo / addatt/mo1 /idnum/ VINT/scalar / nnodes  \n");
			fprintf(fid_lgi, "interpolate/voronoi/mo1,idnum/1 0 0/mo2,n_num  \n");
			fprintf(fid_lgi, "dump / avs2/ClosestNodeR.inp / mo1/0 0 1 0  \n");

			fprintf(fid_lgi,"finish \n");
			fprintf(fid_lgi,"\n");
            fclose(fid_lgi);

//			os.system(lagritpath + " <definedist.lgi >distance.out")
            printf("\n hpham: Make sure you define the path to LAGRIT here at line 446. \n"); 
            system("/projects/DFN/apps/LAGRIT/lagrit_ulin3.2 + <definedist.lgi >distance.out");




			InitInMatrix(); 
			k_new=npart;  
        }
        


















		
      else
        {
         // if particles will be set randomly over all fractures surface
          initfile=Control_File("init_random:",12);
          res=strncmp(initfile.filename,"yes",3);
          if (res==0)      // yes, run this option.
         {
          flag_in=4;      // hpham noted: Particles released randomly over all fracture surfaces.
          initfile = Control_Data("in_randpart:",12 );
          npart=initfile.flag;      //total number of particles
          /*  memory allocatation for particle structures */
          particle=(struct contam*) malloc ((npart+1)*sizeof(struct contam));
          printf("\n hpham: flag_in=4. Initial particles will be distributed randomly over all fracture surfaces \n");
     
          double random_number=0, sum_aperture=0.0; 
          unsigned int currentcell, k_curr=0;
      
          do
            {
             random_number=drand48();      // Equal the number of triangular elements. How? 
             currentcell=random_number*ncells;

             printf("k_curr = %d, random_number = %f, ncells = %d, Current cell = %d \n ", k_curr, random_number, ncells, currentcell );     //hpham
//             printf("Current cell = %d \n", currentcell);     //hpham

             if ((currentcell!=0) && (((node[cell[currentcell-1].node_ind[0]-1].typeN<200)||(node[cell[currentcell-1].node_ind[0]-1].typeN>250))&& ((node[cell[currentcell-1].node_ind[1]-1].typeN<200)||(node[cell[currentcell-1].node_ind[1]-1].typeN>250)) && ((node[cell[currentcell-1].node_ind[2]-1].typeN<200)||(node[cell[currentcell-1].node_ind[2]-1].typeN>250))))
             {
              particle[k_curr].velocity[0]=0.;
              particle[k_curr].velocity[1]=0.;
              particle[k_curr].fracture=cell[currentcell-1].fracture;
              particle[k_curr].cell=currentcell;
              particle[k_curr].time=0.0;    
              Moving2Center (k_curr, currentcell);
              int insc;
              insc=InsideCell (currentcell);     //function returns 1 if particle is inside cell
              sum_aperture=sum_aperture+node[cell[currentcell-1].node_ind[0]-1].aperture;
              k_curr++;
             }
            }     // end of do loop.
          while(k_curr!=npart);
            k_new=k_curr;

          for (i=0; i<npart; i++)
            {
              particle[i].fl_weight=node[cell[particle[i].cell-1].node_ind[0]-1].aperture/sum_aperture;  
            }
        } //end if flag_in=4








          else
        {
          // if particles are set randomly in rock matrix
          initfile=Control_File_Optional("init_matrix:",12);
          res=strncmp(initfile.filename,"yes",3);
          if (res==0)
            {
              printf(" Initially particles are placed in rock matrix randomly. ");
              printf(" The closest cells to initial particles positions "); 
              printf(" will be set as starting point in DFN. ");
              flag_in=5;    
              InitInMatrix(); 
              k_new=npart;
            } //end if/else flag_in=5
		

		else
         {
          // if particles are uniformly distributed over 
		  // all fracture surfaces in a sphere.
          initfile=Control_File_Optional("init_sphere:",12);
          if (initfile.flag>0)
            {
              res=strncmp(initfile.filename,"yes",3);
              if (res==0)
            {
              printf(" Initially particles are uniformly distributed over ");
              printf(" all fracture surfaces in a sphere."); 
              flag_in=7;    
              InitInSphere(); 
              k_new=npart;
            }
            } 
           
         } // end flag_in = 7 - particles placed in sphere.
		
		
		
		
		} //end if flag_in=5
		 
        }    //end if/else flag_in=4
      
    } //end if /else flag_in=6
    } // end if flag_in=3
    } //end if/else flag_in=2    



    
  /* define beginning and end of edge */

  frc=node[nodezonein[0]-1].fracture[0];
  firstn=nodezonein[0];
  lastn=nodezonein[0];

  
  double thirdcoor=0.0;
  int firstcoor=0, secondcoor=0, frc_count=0;
  if (node[nodezonein[1]-1].fracture[0]==frc)
    {
      if (fabs(node[nodezonein[0]-1].coord[0])-fabs(node[nodezonein[1]-1].coord[0])<1e-10)
    {   
      firstcoor=1;
      secondcoor=2;
    } 
      else
    {
      firstcoor=0;
      if (fabs(node[nodezonein[0]-1].coord[1])-fabs(node[nodezonein[1]-1].coord[1])<1e-10)
        secondcoor=2;
      else
        secondcoor=1;
    } 
    }
  
  /*** hpham marked this: Loop on all nodes in zone: define the boundary nodes for each fracture ****/
  k_current=0;
  //  fprintf(inp,"%d\n", nzone_in);
  for (i=0; i<nzone_in-1; i++)     //nzone_in - number of nodes in flow-in zone from boundary.zone
    {
      //    fprintf(inp, " %d %d %d  %d %d %d\n", i,nodezonein[i], node[nodezonein[i]-1].fracture[0], node[nodezonein[i]-1].fracture[1], firstn, lastn);
      if (node[nodezonein[i]-1].fracture[0]==frc)
    {
   
      if (node[nodezonein[i]-1].coord[firstcoor]!=node[firstn-1].coord[firstcoor])
        {
          if (node[nodezonein[i]-1].coord[firstcoor]<node[firstn-1].coord[firstcoor])
        {
          firstn=nodezonein[i];
          first_ind=i;
        }
          if (node[nodezonein[i]-1].coord[firstcoor]>node[lastn-1].coord[firstcoor])
        {
          lastn=nodezonein[i];
          last_ind=i;
        }
          //      fprintf(inp,"first coord, %d %d\n",firstn, lastn);    
        }
      else
        {
          if (node[nodezonein[i]-1].coord[secondcoor]<node[firstn-1].coord[secondcoor])
        {
          firstn=nodezonein[i];
          first_ind=i;
                }
          if (node[nodezonein[i]-1].coord[secondcoor]>node[lastn-1].coord[secondcoor])
        {
          lastn=nodezonein[i];
          last_ind=i;
        }  
          //        fprintf(inp,"second coord, %d %d\n",firstn, lastn);      
        }
    

    }
      //printf("fracture in flow-in zone %d first node %d last node %d \n", frc, firstn, lastn);
      if ((node[nodezonein[i]-1].fracture[0]!=frc)||((i==nzone_in-2)))
    {
      //          printf("fracture in flow-in zone %d first node %d last node %d \n", frc, firstn, lastn);
          if ((i==nzone_in-2)  && (node[firstn-1].fracture[0]==node[nodezonein[nzone_in-1]-1].fracture[0]))
        {
          lastn=nodezonein[nzone_in-1];
          last_ind=nzone_in-1;
          
        }
      if (firstn!=lastn)
        {  
         
          if (flag_in==2)
        {
          
          k_new=InitParticles_eq (k_current, firstn, lastn, parts_dist, first_ind, last_ind);
          //        printf("fract first %d fract last %d number parts %d\n", node[firstn-1].fracture[0], node[lastn-1].fracture[0], k_new);
        }
          if (flag_in==1)
        k_new=InitParticles_np (k_current, firstn, lastn, parts_fracture, first_ind, last_ind);
          
          k_current=k_new;
       
          if (flag_in==3)
        {
            
          double cx1=0, cx2=0, cy1=0, cy2=0;
           
          /* define ends points of fracture edge, then calculate an intersection with starting region */ 
          inter_p[frc_count][0]=1e-10;
          inter_p[frc_count][1]=1e-10;
          inter_p[frc_count][2]=1e-10;
          inter_p[frc_count][3]=1e-10;
          if ((zonenumb_in==1)||(zonenumb_in==2)) 
            {
              cx1=node[firstn-1].coord[0];
              cy1=node[firstn-1].coord[1];
              if ((cx1>ixmin) && (cx1<ixmax) &&(cy1>iymin)&&(cy1<iymax))
            {
              inter_p[frc_count][0]=cx1;
              inter_p[frc_count][1]=cy1; 
            } 
              cx2=node[lastn-1].coord[0];
              cy2=node[lastn-1].coord[1];
              if ((cx2>ixmin) && (cx2<ixmax) &&(cy2>iymin)&&(cy2<iymax))
            {
              inter_p[frc_count][0]=cx2;
              inter_p[frc_count][1]=cy2; 
            } 
                
              thirdcoor= node[firstn-1].coord[2]; 
            }
            
          if ((zonenumb_in==3)||(zonenumb_in==5)) 
            {
              cx1=node[firstn-1].coord[2];
              cy1=node[firstn-1].coord[1];
              if ((cx1>izmin) && (cx1<izmax) &&(cy1>iymin)&&(cy1<iymax))
            {
              inter_p[frc_count][0]=cx1;
              inter_p[frc_count][1]=cy1; 
            }  
               
              cx2=node[lastn-1].coord[2];
              cy2=node[lastn-1].coord[1];
              if ((cx2>izmin) && (cx2<izmax) &&(cy2>iymin)&&(cy2<iymax))
            {
              inter_p[frc_count][0]=cx2;
              inter_p[frc_count][1]=cy2; 
            }
              thirdcoor= node[firstn-1].coord[0];   
            }
            
          if ((zonenumb_in==4)||(zonenumb_in==6)) 
            {
              cx1=node[firstn-1].coord[0];
              cy1=node[firstn-1].coord[2];
              if ((cx1>ixmin) && (cx1<ixmax) &&(cy1>izmin)&&(cy1<izmax))
            {
              inter_p[frc_count][0]=cx1;
              inter_p[frc_count][1]=cy1; 
            }   
              cx2=node[lastn-1].coord[0];
              cy2=node[lastn-1].coord[2];
              if ((cx2>ixmin) && (cx2<ixmax) &&(cy2>izmin)&&(cy2<izmax))
            {
              inter_p[frc_count][0]=cx2;
              inter_p[frc_count][1]=cy2; 
            } 
              thirdcoor= node[firstn-1].coord[1];  
            }
            
          double pr1, pr2, pr3, pr4, p_x, p_y;
          int ii;
     
          /* define intersection points of boundary fracture edges and starting region sides*/ 
          for (ii=0; ii<4; ii++)
            {
              int kk;
              kk=ii+1;
              if (ii==3)
            kk=0;
              pr1=(px[ii]-cx1)*(py[kk]-cy1)-(py[ii]-cy1)*(px[kk]-cx1);
              pr2=(cx1-px[ii])*(cy2-py[ii])-(cy1-py[ii])*(cx2-px[ii]);
              pr3=(px[ii]-cx2)*(py[kk]-cy2)-(py[ii]-cy2)*(px[kk]-cx2);
              pr4=(cx1-px[kk])*(cy2-py[kk])-(cy1-py[kk])*(cx2-px[kk]);

              if ((pr1*pr3<0)&&(pr2*pr4<0))
            {
             
              pr1=cx1*cy2-cy1*cx2;
              pr2=px[ii]*py[kk]-py[ii]*px[kk];
              pr3=(cx1-cx2)*(py[ii]-py[kk])-(cy1-cy2)*(px[ii]-px[kk]);
              p_x=((px[ii]-px[kk])*pr1-(cx1-cx2)*pr2)/pr3;
              p_y=((py[ii]-py[kk])*pr1-(cy1-cy2)*pr2)/pr3;
              if (inter_p[frc_count][0]==1e-10)
                {
                  inter_p[frc_count][0]=p_x;
                  inter_p[frc_count][1]=p_y;
                }
              else
                {  
                  inter_fr[frc_count]=frc; 
                  inter_p[frc_count][2]=p_x;
                  inter_p[frc_count][3]=p_y;   
                  frc_count++; 
                  break;
                }
       
            }
             
            }

   
        } //end of flag
            }
            
      if (i<nzone_in-2)
        {
          frc=node[nodezonein[i]-1].fracture[0];
         
          firstn=nodezonein[i];
          first_ind=i;
          lastn=nodezonein[i];
          last_ind=i;
          //          fprintf(inp,"p %d %d %d %d\n", i, frc, firstn, lastn); 
          if (node[nodezonein[i+1]-1].fracture[0]==frc)
        {
          if (abs(node[nodezonein[i]-1].coord[0])-abs(node[nodezonein[i+1]-1].coord[0])<1e-10)
            {   
              firstcoor=1;
              secondcoor=2;
            } 
          else
            {
              firstcoor=0;
   
              if (abs(node[nodezonein[i]-1].coord[1])-abs(node[nodezonein[i+1]-1].coord[1])<1e-10)
      
            secondcoor=2;
              else
            secondcoor=1;
            } 
        }
        }
        
       
    }
    }     // end of for (i=0; i<nzone_in-1; i++)     //nzone_in - number of nodes in flow-in zone from boundary.zone



    //   printf("fracture in flow-in zone %d first node %d last node %d \n", frc, firstn, lastn);
  if ((flag_in==3)&&(frc_count==0))
    {
      printf("\n There is no fracture crosses the given range. Try to increase the range. \n");
      printf("\n Program is terminated. \n");
      exit(1);  
    } 

    
  /* hpham's highlight: Place particles in the starting region. Every fracture (or part of fracture edge) 
                        will have the same amount of particles */ 
  if (flag_in==3)
    {
      for (i=0; i<frc_count; i++)
        {
            parts_fracture=(int) npart/frc_count;
            k_new=InitParticles_ones (k_current, inter_p, inter_fr[i], parts_fracture, i, thirdcoor, zonenumb_in, first_ind, last_ind);     //hpham
            k_current=k_new;
            //   printf(" %d intersects at %f %f %f %f\n",i, inter_p[i][0], inter_p[i][1], inter_p[i][2], inter_p[i][3]);
        }
    }










  // fclose(inp);
  if (flag_in==0)
    {
      printf("\n There is no option specified for particles initial positions! \n");
      printf("\n Program is terminated. \n");
      exit(1);
    }
    printf("\n k_new = %d \n", k_new);
//  k_new = nnode_iwel;     //hpham added (may not correct. ). 
  return k_new;
}




/////////////////////////////////////////////////////////////////////////////

int InitCell ()
{
  /** Function defines particle's initial cell ******/ 

//  printf(" np = %4d, particle[np].fracture = %6d \n", np, particle[np].fracture);

//  printf(" fracture[particle[np].fracture-1].numbcells = %d \n", fracture[particle[np].fracture-1].numbcells);

  int i, curcel, insc=0;  
  for (i=0; i<fracture[particle[np].fracture-1].numbcells; i++)
    {
      curcel=fracture[particle[np].fracture-1].firstcell+i;
      insc=InsideCell (curcel);

      if (insc==1)
    {
//      printf("Particle %d is in initial cell number is %d in fracture %d\n",np+1, particle[np].cell, particle[np].fracture);
      break;
    }
    }
  return insc;  
}
////////////////////////////////////////////////////////////////////////////

int InitParticles_np (int k_current, int firstn, int lastn, int parts_fracture, int first_ind, int last_ind)
{
  /***********function defines particle's initial positions inside one fracture ***/
  /*** here the same number of particles equdists inside one edge *****************/
  double deltax, deltay;
  int j, pf;
  pf=parts_fracture;
  deltax=(node[lastn-1].coord_xy[0]-node[firstn-1].coord_xy[0]);
  deltay=(node[lastn-1].coord_xy[1]-node[firstn-1].coord_xy[1]);
  
  // printf("%d fr. first last coordinates, [%f %f],[%f %f] \n", node[firstn-1].fracture[0], node[lastn-1].coord_xy[0],node[firstn-1].coord_xy[0],node[lastn-1].coord_xy[1],node[firstn-1].coord_xy[1]);
  for (j=0; j<pf; j++)
    { 
      particle[k_current].position[0]=node[firstn-1].coord_xy[0]+(deltax/pf)*(j)+deltax/(2.0*pf);
      particle[k_current].position[1]=node[firstn-1].coord_xy[1]+(deltay/pf)*(j)+deltay/(2.0*pf);
      //    printf("%f %f\n", particle[k_current].position[0],particle[k_current].position[1]);
 
      // first particle will be on boundary node
      //      particle[k_current].position[0]=node[firstn-1].coord_xy[0]+(deltax/(pf-1))*(j);
      //     particle[k_current].position[1]=node[firstn-1].coord_xy[1]+(deltay/(pf-1))*(j);

      particle[k_current].velocity[0]=0.;
      particle[k_current].velocity[1]=0.;
      particle[k_current].fracture=node[firstn-1].fracture[0];
      //      printf("%d %d %d %d \n", k_current, firstn, lastn, node[firstn-1].fracture[0]);
      particle[k_current].intcell=0;
      particle[k_current].time=0.0; 
      if (flag_w==1)
    particle[k_current].fl_weight=0.0;
      else
    particle[k_current].fl_weight=0.;
   

      k_current++;
   
      if (k_current>npart)
    { 
      printf(" \n Number of particles with allocated memory is less than number of particles set up initially. \n");
      printf(" Increase the number of particles. Program is terminated. \n");
      exit(1);
      
    }
    }
  return k_current;
}
////////////////////////////////////////////////////////////////////////////

int InitParticles_eq (int k_current, int firstn, int lastn, double parts_dist, int first_ind, int last_ind)
{
  /***********function defines particle's initial positions inside one fracture ***/
  /***** here the given distance between particles dictates how many particles **/
  /******************will be placed in one fracture*****************************/
  double deltax, deltay, edgelength, eqdist_x, eqdist_y;
  unsigned int j, pf;
  
  deltax=(node[lastn-1].coord_xy[0]-node[firstn-1].coord_xy[0]);
  deltay=(node[lastn-1].coord_xy[1]-node[firstn-1].coord_xy[1]);
  edgelength=sqrt(deltax*deltax+deltay*deltay);
  pf=(int)(edgelength/parts_dist);
  if (pf<2)
    {
      pf=1;
      eqdist_x=deltax/2.0;
      eqdist_y=deltay/2.0;
    }
  else
    {
      eqdist_x=deltax/pf;
      eqdist_y=deltay/pf;
    }
  // printf("edge %f parts_dist %f eqdist %f %f pf %d \n", edgelength, parts_dist, eqdist_x, eqdist_y, pf);
  // printf("%d fr. first last coordinates, [%f %f],[%f %f] \n", node[firstn-1].fracture[0], node[lastn-1].coord_xy[0],node[firstn-1].coord_xy[0],node[lastn-1].coord_xy[1],node[firstn-1].coord_xy[1]);
  for (j=0; j<pf; j++)
    { 
  
      particle[k_current].position[0]=node[firstn-1].coord_xy[0]+eqdist_x*(j)+eqdist_x/2.0;
      particle[k_current].position[1]=node[firstn-1].coord_xy[1]+eqdist_y*(j)+eqdist_y/2.0;
 
 
      // first particle will be on boundary node
      //      particle[k_current].position[0]=node[firstn-1].coord_xy[0]+(deltax/(pf-1))*(j);
      //      particle[k_current].position[1]=node[firstn-1].coord_xy[1]+(deltay/(pf-1))*(j);

      particle[k_current].velocity[0]=0.;
      particle[k_current].velocity[1]=0.;
      particle[k_current].fracture=node[firstn-1].fracture[0];
      particle[k_current].intcell=0;
      particle[k_current].time=0.0; 
      if (flag_w==1)
    particle[k_current].fl_weight=0.0;
      else
    particle[k_current].fl_weight=0.;
    
 
      k_current++;
   
      if (k_current>npart)
    { 
      printf("\n Number of particles with allocated memory is less than number of particles set up initially. \n");
      printf(" Increase the number of particles. Program is terminated. \n");
      exit(1);
    }
    }
  return k_current;
}
////////////////////////////////////////////////////////////////////////////
//k_new=InitParticles_ones (k_current,        inter_p,          inter_fr[i],    parts_fracture,     i,         thirdcoor,     zonenumb_in,     first_ind,     last_ind); 
int InitParticles_ones (int k_current, double inter_p[][4], int fracture_n, int parts_fracture, int ii, double thirdcoor, int zonenumb_in, int first_ind, int last_ind)
{
  /***********function defines particle's initial positions  ***/
  /****************** in one starting point*****************/
  int j=fracture_n-1;
  double x_1=0, y_1=0, z_1=0, x_2=0, z_2=0, y_2=0;
  double x1cor=0.0, y1cor=0.0, z1cor=0, x2cor=0.0, y2cor=0.0, z2cor=0;
  printf("k_current = %d | fracture_n = %d | parts_fracture = %d | ii = %d | thirdcoor = %7.4f | zonenumb_in = %d | first_ind = %d | last_ind = %d \n", k_current, fracture_n, parts_fracture, ii, thirdcoor, zonenumb_in, first_ind, last_ind);
  if ((zonenumb_in==1)||(zonenumb_in==2)) 
    {
      x1cor=inter_p[ii][0];
      y1cor=inter_p[ii][1];
      z1cor=thirdcoor;
      x2cor=inter_p[ii][2];
      y2cor=inter_p[ii][3];
      z2cor=thirdcoor;
  
    }
            
  if ((zonenumb_in==3)||(zonenumb_in==5)) 
    {
      x1cor=thirdcoor;
      x2cor=thirdcoor;
      z1cor=inter_p[ii][0];
      z2cor=inter_p[ii][2];
      y1cor=inter_p[ii][1];
      y2cor=inter_p[ii][3];
              
    }
            
  if ((zonenumb_in==4)||(zonenumb_in==6)) 
    {
      y1cor=thirdcoor;
      y2cor=thirdcoor;
      x1cor=inter_p[ii][0];
      x2cor=inter_p[ii][2];
      z1cor=inter_p[ii][1];
      z2cor=inter_p[ii][3];

    }
         
  if (fracture[j].theta!=0.0)
    {
    
      x_1=fracture[j].rot2mat[0][0]*x1cor+fracture[j].rot2mat[0][1]*y1cor+fracture[j].rot2mat[0][2]*z1cor;
      y_1=fracture[j].rot2mat[1][0]*x1cor+fracture[j].rot2mat[1][1]*y1cor+fracture[j].rot2mat[1][2]*z1cor; 
      z_1=fracture[j].rot2mat[2][0]*x1cor+fracture[j].rot2mat[2][1]*y1cor+fracture[j].rot2mat[2][2]*z1cor; 
      x_2=fracture[j].rot2mat[0][0]*x2cor+fracture[j].rot2mat[0][1]*y2cor+fracture[j].rot2mat[0][2]*z2cor;
      y_2=fracture[j].rot2mat[1][0]*x2cor+fracture[j].rot2mat[1][1]*y2cor+fracture[j].rot2mat[1][2]*z2cor; 
      z_2=fracture[j].rot2mat[2][0]*x2cor+fracture[j].rot2mat[2][1]*y2cor+fracture[j].rot2mat[2][2]*z2cor; 
      
    }
  else 
    {
      /* if angle =0 and fracture is parallel to xy plane, we use the same x and y coordinates */
      x_1=x1cor;
      y_1=y1cor;
      z_1=z1cor;
      x_2=x2cor;
      y_2=y2cor;
      z_2=z2cor;
      
    }
  double deltax, deltay;
  unsigned int  pf;
  pf=parts_fracture;
  deltax=(x_2-x_1);
  deltay=(y_2-y_1); 
  
  for (j=0; j<pf; j++)
    { 
  
      particle[k_current].position[0]=x_1+(deltax/pf)*(j)+deltax/(2.0*pf);
   
      particle[k_current].position[1]=y_1+(deltay/pf)*(j)+deltay/(2.0*pf);
      
      particle[k_current].velocity[0]=0.;
      particle[k_current].velocity[1]=0.;
      particle[k_current].fracture=fracture_n;
      particle[k_current].intcell=0;
      particle[k_current].time=0.0; 
      particle[k_current].fl_weight=0.0;
      k_current++;
   
      if (k_current>npart)
    { 
      printf(" \n Number of particles with allocated memory is less than number of particles set up initially. \n");
      printf(" Increase the number of particles. Program is terminated. \n");
      exit(1);
      
    }
    }
  return k_current;
}








////////////////////////////////////////////////////////////////////////////
void FlowInWeight(int numberpart)
{
  /*** function defines weights of particles based on *********************/
  /*** in-flow flux boundary cells ****************************************/
  int ind1=0,ind2=0,ver1=0, ver2=0, ver3=0,incell=0, n1in=0, n2in=0, jj;
  int ins;
  double sumflux1=0, sumflux2=0, particleflux[numberpart], totalflux=0;
  for (np=0; np<numberpart; np++) 
    {
      ins=0;
      ins=InitCell();
      incell=particle[np].cell;
      if (incell!=0)
        {
      incell=particle[np].cell;
      ver1=cell[incell-1].node_ind[0];
      ver2=cell[incell-1].node_ind[1];
      ver3=cell[incell-1].node_ind[2];
      n1in=0;
      n2in=0;
      if (node[ver1-1].typeN>=300)
            {
          n1in=ver1;
          ind1=0;
            }
      if (node[ver2-1].typeN>=300)
            {
          if (n1in==0)
        {
          n1in=ver2;
          ind1=1;
                }  
          else
                {
          n2in=ver2;
          ind2=1;
                } 
        }
      if (node[ver3-1].typeN>=300)
            {
          if (n1in==0)
        {
          n1in=ver3;
          ind1=2;
                }  
          else
                {
          n2in=ver3;
          ind2=2;
                }  
                
        }  
      //         printf("%d %d %d %d \n", np+1, n1in, n2in, incell);
      if ((n1in!=0) && (n2in!=0))
        {
          sumflux1=0;
          sumflux2=0;
          for (jj=0; jj<node[n1in-1].numneighb; jj++)
                  
            sumflux1=sumflux1+fabs(node[n1in-1].flux[jj]);
               
          for (jj=0; jj<node[n2in-1].numneighb; jj++)
                  
            sumflux2=sumflux2+fabs(node[n2in-1].flux[jj]);
                  
               
          particleflux[np]=particle[np].weight[ind1]*sumflux1+particle[np].weight[ind2]*sumflux2;
          totalflux=totalflux+particleflux[np];   
                  
                  
          //              printf("%d   %5.12e %5.12e  %5.12e  %5.12e  %5.12e\n", np+1,particle[np].weight[ind1],sumflux1, particle[np].weight[ind2],sumflux2,particleflux[np] );
                  
        }   
      else
        {
          int ncent=0;
          ncent=n1in+n2in;
          if (ncent!=0)
        {
          sumflux1=0;
               
          for (jj=0; jj<node[ncent-1].numneighb; jj++)
                  
            sumflux1=sumflux1+fabs(node[ncent-1].flux[jj]);
               
               
          particleflux[np]=sumflux1;
          totalflux=totalflux+particleflux[np];   
                  
                  
          //             printf("%d   %5.12e \n", np+1,particleflux[np] );
        }
        }
               
    }
    }
  
  for (np=0; np<numberpart; np++) 
    {
      particle[np].fl_weight=particleflux[np]/totalflux;
      // printf("%d   %5.12e %5.12e  %5.12e  \n", np+1,particleflux[np], totalflux, particle[np].fl_weight);
    }
  

 
  return;
}

//////////////////////////////////////////////////////////////////////////
void InitInMatrix()
{
  /**** function read files with data for particles initially placed in rock matrix***/
  struct inpfile inputfile;
  
  inputfile = Control_File("inm_coord:",10 );

  FILE *mc= OpenFile (inputfile.filename,"r");
  
  printf("\n OPEN AND READ FILE: %s \n \n", inputfile.filename);

  inputfile = Control_File("inm_nodeID:",10 );

  FILE *mn= OpenFile (inputfile.filename,"r");
  
  printf("\n OPEN AND READ FILE: %s \n \n", inputfile.filename);
 
  //FILE *mdf=OpenFile("distance_time.dat","w");
  int i, ii;
  unsigned int number, no, npart_name; //hpham add npart_name
  char cs;
  if (fscanf(mn,"%d  %d %d %d %d \n", &no, &no, &number, &no, &no)!=5)
             printf("error");

  printf("\n hpham: number= %d \n", number); 

  for (i=0; i<(number+1); i++)
    {
      do 
    cs=fgetc(mn); 
      while (cs!='\n');
    }
      
   
  //  if (number !=npart)
  //  {
  //    printf(" The numbers of particles in input files doesn't match. \n");  
  //    printf(" Program is terminated. \n");
  //    exit(1);
  //   }
          if (fscanf(mc," %s %d \n", &npart_name, &npart)!=2)   //hpham modified
             printf("error reading ParticleInitCoordR.dat"); 

  /*  memory allocatation for particle structures */
        
  particle=(struct contam*) malloc ((npart+1)*sizeof(struct contam));
     

  double xp, yp, zp,  sum_distance=0.0, xp2,yp2,zp2;
  
  printf("\n hpham: npar= %d \n", npart);
  
  for (ii=0; ii<npart; ii++)
    {
//      printf("hpham: ipar = %d \n",ii);
      if (fscanf(mn,"%d  %d   %d  %d  %d  %f %d \n",&no, &no, &no, &no, &no, &no, &number)!=7) // hpham 7 cols instead of 6 columns in the org version
           printf("error at line 1202 \n");     // when read ClosestNodeR.inp

      if (fscanf(mc," %lf %lf %lf \n ", &xp, &yp, &zp)!=3)
             printf ("error reading file ParticleInitCoordR.dat \n");
//      printf("hpham: ipar, xp, yp, zp = %d %f %f %f \n",ii, xp, yp, zp);                     


//      printf("hpham: ipar, xp2, yp2, zp2 = %d %f %f %f %d\n",ii, xp2, yp2, zp2,number);


      particle[ii].velocity[0]=0.;
      particle[ii].velocity[1]=0.;
      particle[ii].fracture=node[number-1].fracture[0];
      particle[ii].cell=0;
      particle[ii].position[0]=node[number-1].coord_xy[0];
      particle[ii].position[1]=node[number-1].coord_xy[1];
      particle[ii].time=0.0;    

    }       
  fclose(mc);
  fclose(mn);

  return;
}









//////////////////////////////////////////////////////////////////////////
void InitInSphere()
{
  /**** function read files with data for particles initially placed in sphere ***/
  struct inpfile inputfile;
  
  inputfile = Control_File("insphr_coord:",13 );

  FILE *mc= OpenFile (inputfile.filename,"r");
  
  printf("\n HPham: OPEN AND READ FILE: %s \n \n", inputfile.filename);

  int i, ii;
  unsigned int number, no, frac_ID; //hpham add npart_name

 
  if (fscanf(mc," %d \n", &npart)!=1)   //hpham modified
	 printf("error reading npart in ParticleInitCoordR.dat opt 7: sphere"); 

  /*  memory allocatation for particle structures */      
  particle=(struct contam*) malloc ((npart+1)*sizeof(struct contam));
 
  printf("hpham: npar for opt 7 sphere = %d \n", npart);
  
//  double* distance=malloc ((npart+1)*sizeof(double));
  double xp, yp, zp,  sum_distance=0.0, xp2,yp2,zp2;
    
  for (ii=0; ii<npart; ii++)
    {
      if (fscanf(mc," %lf %lf %lf %d %d \n ", &xp, &yp, &zp, &number, &frac_ID )!=5)
             printf ("error reading file ParticleInitCoordR.dat \n");
//      printf("hpham: ipar, xp, yp, zp, ele_ID, frac_ID = %d %f %f %f %d %d \n",ii, xp, yp, zp, number, frac_ID);                     
//      printf("hpham: node[number-1].coord[0] = %f \n", node[number-1].coord[0]);
//	  xp2=(node[number-1].coord[0]-xp)*(node[number-1].coord[0]-xp);
//      yp2=(node[number-1].coord[1]-yp)*(node[number-1].coord[1]-yp);
//      zp2=(node[number-1].coord[2]-zp)*(node[number-1].coord[2]-zp);
//      distance[ii]=sqrt(xp2+yp2+zp2);
//      printf("hpham: ii, distance = %d %f \n",ii, distance[ii]);                     
	  
      particle[ii].velocity[0]=0.;
      particle[ii].velocity[1]=0.;

//	  particle[ii].fracture=node[number+1].fracture[0];
//      printf("hpham: number, frac_IDs = %d %d %d \n",number, node[number-1].fracture[0], frac_ID);

      particle[ii].fracture=frac_ID;
      particle[ii].cell=number;
      particle[ii].position[0]=xp;
      particle[ii].position[1]=yp;
      particle[ii].time=0.0;    
      // define a time that took for particle to reach fracture from a rock matrix 
//      particle[ii].time=TimeFromMatrix(distance[ii]);
//      sum_distance=sum_distance+distance[ii];
	  printf("hpham: ipar, xp, yp, ele_ID, frac_ID = %d %f %f %d %d \n",ii, particle[ii].position[0], particle[ii].position[1], particle[ii].cell, particle[ii].fracture);

    }       
  fclose(mc);
  //define weights according to distance from fracture
//  for (i=0; i<npart; i++)
//    {
//      particle[i].fl_weight=distance[i]/sum_distance;
//    }
  
//  free(distance);
  return;
}
///////////////////////////////////////////////////////////////////////////





























///////////////////////////////////////////////////////////////////////////
double TimeFromMatrix(double pdist)
{
  /* function defines the time that is required for particle to travel from rock 
     matrix to fracture. This will be particle's initial time */
  struct inpfile inputfile;
  double ptime=0.0, ptime1=0.0, ptime2=0.0;
  double randomnumber=0.0;
  randomnumber=drand48();
  // if (randomnumber==1.0)
  //      randomnumber=0.99999; 
  double mporosity=0.0;
  double mdiffcoeff=0.0;
  
  inputfile=Control_Param("inm_porosity:",13);
  mporosity=inputfile.param;
  inputfile=Control_Param("inm_diffcoeff:",14);
  mdiffcoeff=inputfile.param;
  //  printf("diff coeff %5.9e porosity %lf \n",mdiffcoeff, mporosity);
  
  double retardation_factor=1.0;
  double inverse_erfc=0.0;
  double z;
  z=1-randomnumber;
  inverse_erfc=0.5*sqrt(pi)*(z+(pi/12)*pow(z,3)+((7*pow(pi,2))/480)*pow(z,5)+((127*pow(pi,3))/40320)*pow(z,7)+((4369*pow(pi,4))/5806080)*pow(z,9)+((34807*pow(pi,5))/182476800)*pow(z,11));
  ptime2=(1.0/inverse_erfc)*(1.0/inverse_erfc);
  ptime1=(pdist*pdist)/(4.0*(mporosity*mdiffcoeff/retardation_factor));
  ptime=(ptime1*ptime2)/timeunit;
  
        
   
  return ptime;
} 
//////////////////////////////////////////////////////////////////////////

// Copyright 2019, Nevada System of Higher Education on Behalf of the Desert Research Institute, and 
// Copyright 2018, Triad National Secuirty, LLC. All rights reserved.

// This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser 
// General Public License as published by the Free Software Foundation; either version 3.0 of the License, 
// or (at your option) any later version.  This software requires you separately obtain dfnWorks under an 
// appropriate license from the Los Alamos National Laboratory (LANL), which is operated by 
// Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.