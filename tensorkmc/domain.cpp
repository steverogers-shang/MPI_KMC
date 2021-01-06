/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "domain.h"
#include "app.h"
#include "lattice.h"
#include "memory.h"
#include "error.h"

#include "style_region.h"

using namespace SPPARKS_NS;

#define DELTA 1

/* ---------------------------------------------------------------------- */

Domain::Domain(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  dimension = 3;
  xperiodic = yperiodic = zperiodic = 1;
  nonperiodic = 0;
  xperiodic = yperiodic = zperiodic = 1;
  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  nx = ny = nz = 0;
  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;

  box_exist = 0;

  lattice = NULL;
  nregion = maxregion = 0;
  regions = NULL;
}

/* ---------------------------------------------------------------------- */

Domain::~Domain()
{
  delete lattice;
  for (int i = 0; i < nregion; i++) delete regions[i];
  memory->sfree(regions);
  // MPI_Comm_free(&comm_3d);
}

/* ----------------------------------------------------------------------
   setup global box
   assumes boxlo/hi are already set
------------------------------------------------------------------------- */

void Domain::set_box()
{
  if (boxxlo >= boxxhi || boxylo >= boxyhi || boxzlo >= boxzhi)
    error->one(FLERR,"Box bounds are invalid");

  xprd = boxxhi - boxxlo;
  yprd = boxyhi - boxylo;
  zprd = boxzhi - boxzlo;
}

/* ----------------------------------------------------------------------
   create a lattice
   delete it if style = none
------------------------------------------------------------------------- */

void Domain::set_lattice(int narg, char **arg)
{
  if (lattice) delete lattice;
  lattice = new Lattice(spk,narg,arg);
  if (lattice->style == 0) {
    delete lattice;
    lattice = NULL;
  }
}

/* ----------------------------------------------------------------------
   create a new region 
------------------------------------------------------------------------- */

void Domain::add_region(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal region command");

  if (find_region(arg[0]) >= 0) error->all(FLERR,"Reuse of region ID");

  // extend Region list if necessary

  if (nregion == maxregion) {
    maxregion += DELTA;
    regions = (Region **) 
      memory->srealloc(regions,maxregion*sizeof(Region *),"domain:regions");
  }

  // create the Region

  if (strcmp(arg[1],"none") == 0) error->all(FLERR,"Invalid region style");

#define REGION_CLASS
#define RegionStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) \
    regions[nregion] = new Class(spk,narg,arg);
#include "style_region.h"
#undef REGION_CLASS

  else error->all(FLERR,"Invalid region style");

  nregion++;
}

/* ----------------------------------------------------------------------
   return region index if name matches existing region ID
   return -1 if no such region
------------------------------------------------------------------------- */

int Domain::find_region(char *name)
{
  for (int iregion = 0; iregion < nregion; iregion++)
    if (strcmp(name,regions[iregion]->id) == 0) return iregion;
  return -1;
}

/* ----------------------------------------------------------------------
   boundary settings from the input script 
------------------------------------------------------------------------- */

void Domain::set_boundary(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal boundary command");

  if (strcmp(arg[0],"n") == 0) xperiodic = 0;
  else if (strcmp(arg[0],"p") == 0) xperiodic = 1;
  else error->all(FLERR,"Illegal boundary command");
  if (strcmp(arg[1],"n") == 0) yperiodic = 0;
  else if (strcmp(arg[1],"p") == 0) yperiodic = 1;
  else error->all(FLERR,"Illegal boundary command");
  if (strcmp(arg[2],"n") == 0) zperiodic = 0;
  else if (strcmp(arg[2],"p") == 0) zperiodic = 1;
  else error->all(FLERR,"Illegal boundary command");

  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  nonperiodic = 0;
  if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) nonperiodic = 1;

  if (nonperiodic && app->appclass != App::LATTICE)
    error->all(FLERR,"Boundary command currently only supported by on-lattice apps");
}

/* ----------------------------------------------------------------------
   assign nprocs to 1d box as equal partitions
------------------------------------------------------------------------- */

void Domain::procs2domain_1d()
{
  if (user_procgrid[0] || user_procgrid[1] || user_procgrid[2]) {
    if (user_procgrid[1] != 1 || user_procgrid[2] != 1)
      error->all(FLERR,"App style proc count is not valid for 1d simulation");
    procgrid[0] = user_procgrid[0];
  } else {
    procgrid[0] = nprocs;
  }

  procgrid[1] = procgrid[2] = 1;

  myloc[0] = me;
  myloc[1] = myloc[2] = 0;

  subxlo = boxxlo + myloc[0] * xprd/procgrid[0];
  if (myloc[0] < procgrid[0]-1) 
    subxhi = boxxlo + (myloc[0]+1) * xprd/procgrid[0];
  else subxhi = boxxhi;

  subylo = boxylo;
  subyhi = boxyhi;
  subzlo = boxzlo;
  subzhi = boxzhi;
}

/* ----------------------------------------------------------------------
   assign nprocs to 2d box so as to minimize perimeter per proc
------------------------------------------------------------------------- */

void Domain::procs2domain_2d()
{
  int ipx,ipy;
  double boxx,boxy,surf;

  if (user_procgrid[0] || user_procgrid[1] || user_procgrid[2]) {
    if (user_procgrid[2] != 1)
      error->all(FLERR,"App style proc count is not valid for 2d simulation");
    procgrid[0] = user_procgrid[0];
    procgrid[1] = user_procgrid[1];

  } else {

    // loop thru all possible factorizations of nprocs
    // surf = perimeter of a proc sub-domain

    double bestsurf = 2.0 * (xprd+yprd);
 
    ipx = 1;
    while (ipx <= nprocs) {
      if (nprocs % ipx == 0) {
	ipy = nprocs/ipx;
	boxx = xprd/ipx;
	boxy = yprd/ipy;
	surf = boxx + boxy;
	if (surf < bestsurf) {
	  bestsurf = surf;
	  procgrid[0] = ipx;
	  procgrid[1] = ipy;
	}
      }
      ipx++;
    }
  }

  procgrid[2] = 1;
  // delete by smz 170104, use MPI_Cart_create instead.
  myloc[0] = me % procgrid[0];
  myloc[1] = me/procgrid[0];
  myloc[2] = 0;

  // int my2drank;  
  // int periods[3]={1,1,1};  

  // MPI_Cart_create(world,3,procgrid,periods,1,&comm_3d);  
  // MPI_Comm_rank(comm_3d, &my2drank);  
  // MPI_Cart_coords(comm_3d, my2drank, 3, myloc); 
  // if (my2drank != me) printf("me:%d, my2drank:%d\n", me, my2drank);
  
  subxlo = boxxlo + myloc[0] * xprd/procgrid[0];
  if (myloc[0] < procgrid[0]-1) 
    subxhi = boxxlo + (myloc[0]+1) * xprd/procgrid[0];
  else subxhi = boxxhi;

  subylo = boxylo + myloc[1] * yprd/procgrid[1];
  if (myloc[1] < procgrid[1]-1) 
    subyhi = boxylo + (myloc[1]+1) * yprd/procgrid[1];
  else subyhi = boxyhi;

  subzlo = boxzlo;
  subzhi = boxzhi;
}

/* ----------------------------------------------------------------------
   assign nprocs to 3d box so as to minimize surface area per proc
------------------------------------------------------------------------- */

void Domain::procs2domain_3d()
{
  int ipx,ipy,ipz,nremain;
  double boxx,boxy,boxz,surf;

  if (user_procgrid[0] || user_procgrid[1] || user_procgrid[2]) {
    procgrid[0] = user_procgrid[0];
    procgrid[1] = user_procgrid[1];
    procgrid[2] = user_procgrid[2];

  } else {

    double bestsurf = 2.0 * (xprd*yprd + yprd*zprd + zprd*xprd);
  
    // loop thru all possible factorizations of nprocs
    // surf = surface area of a proc sub-domain

    ipx = 1;
    while (ipx <= nprocs) {
      if (nprocs % ipx == 0) {
	nremain = nprocs/ipx;
	ipy = 1;
	while (ipy <= nremain) {
	  if (nremain % ipy == 0) {
	    ipz = nremain/ipy;
	    boxx = xprd/ipx;
	    boxy = yprd/ipy;
	    boxz = zprd/ipz;
	    surf = boxx*boxy + boxy*boxz + boxz*boxx;
	    if (surf < bestsurf) {
	      bestsurf = surf;
	      procgrid[0] = ipx;
	      procgrid[1] = ipy;
	      procgrid[2] = ipz;
	    }
	  }
	  ipy++;
	}
      }
      ipx++;
    }
  }

  int XN, XP, YN, YP, ZN, ZP;
  // delete by smz 170104, use MPI_Cart_create instead, myloc is different.
  myloc[0] = me % procgrid[0];
  myloc[1] = (me/procgrid[0]) % procgrid[1];
  myloc[2] = me / (procgrid[0]*procgrid[1]);
  // printf("1.me:%d, (%d, %d, %d)\n", me , myloc[0], myloc[1], myloc[2]);
  int px, py, pz;
  px = myloc[0];
  py = myloc[1];
  pz = myloc[2];
  ipx = procgrid[0];
  ipy = procgrid[1];
  ipz = procgrid[2];

  XN = me - 1;
  XP = me + 1;
  if (px == 0) XN = me - 1 + ipx;
  if (px == ipx - 1) XP = me + 1 - ipx;

  YN = me - ipx;
  YP = me + ipx;
  if (py == 0) YN = me + (-1 + ipy) * ipx;
  if (py == ipy - 1) YP = me + (1 - ipy) * ipx;

  ZN = me - ipx * ipy;
  ZP = me + ipx * ipy;
  if (pz == 0) ZN = me + (-1 + ipz) * ipx * ipy;
  if (pz == ipz -1) ZP = me + (1 - ipz) * ipx *ipy;

  // printf("Man [%d %d %d] %d: %d %d %d %d %d %d\n", myloc[0], myloc[1], myloc[2], me, XN, XP, YN, YP, ZN, ZP);

  // int my2drank;  
  // int periods[3]={1,1,1};    
  // MPI_Cart_create(world,3,procgrid,periods,1,&comm_3d);  
  // MPI_Comm_rank(comm_3d, &my2drank);  
  // MPI_Cart_coords(comm_3d, my2drank, 3, myloc); 

  // int direct=0, displ=1;
  // MPI_Cart_shift(comm_3d, direct, displ, &XN, &XP);
  // direct=1;
  // MPI_Cart_shift(comm_3d, direct, displ, &YN, &YP);
  // direct=2;
  // MPI_Cart_shift(comm_3d, direct, displ, &ZN, &ZP);  

  procneigh[0][0] = XN;
  procneigh[0][1] = XP;
  procneigh[1][0] = YN;
  procneigh[1][1] = YP;
  procneigh[2][0] = ZN;
  procneigh[2][1] = ZP;

  subxlo = boxxlo + myloc[0] * xprd/procgrid[0];
  if (myloc[0] < procgrid[0]-1) 
    subxhi = boxxlo + (myloc[0]+1) * xprd/procgrid[0];
  else subxhi = boxxhi;

  subylo = boxylo + myloc[1] * yprd/procgrid[1];
  if (myloc[1] < procgrid[1]-1) 
    subyhi = boxylo + (myloc[1]+1) * yprd/procgrid[1];
  else subyhi = boxyhi;

  subzlo = boxzlo + myloc[2] * zprd/procgrid[2];
  if (myloc[2] < procgrid[2]-1) 
    subzhi = boxzlo + (myloc[2]+1) * zprd/procgrid[2];
  else subzhi = boxzhi;
}

/* ----------------------------------------------------------------------
   wrap xyz coords back into central domain
 ------------------------------------------------------------------------- */

void Domain::pbcwrap(double* xyz)
{

  // x coord

  if (xperiodic) {
    while (xyz[0] < boxxlo)
      xyz[0] += xprd;
    while (xyz[0] >= boxxhi)
      xyz[0] -= xprd;
  }
  
  // y coord
  
  if (yperiodic) {
    while (xyz[1] < boxylo)
      xyz[1] += yprd;
    while (xyz[1] >= boxyhi)
      xyz[1] -= yprd;
  }
  
  // z coord
  
  if (zperiodic) {
    while (xyz[2] < boxzlo)
      xyz[2] += zprd;
    while (xyz[2] >= boxzhi)
      xyz[2] -=   zprd;
  }
}

/* ----------------------------------------------------------------------
   shift xyz2 to periodic image closest to xyz1
 ------------------------------------------------------------------------- */

void Domain::pbcshift(double* xyz1, double* xyz2)
{

  // x coord

  if (xperiodic) {
    while ((xyz2[0] - xyz1[0])*2.0 < xprd)
	 xyz2[0] += xprd;
    while ((xyz2[0] - xyz1[0])*2.0 > xprd)
	 xyz2[0] -= xprd;
  }

  // y coord

  if (yperiodic) {
    while ((xyz2[1] - xyz1[1])*2.0 < yprd)
	 xyz2[1] += yprd;
    while ((xyz2[1] - xyz1[1])*2.0 > yprd)
	 xyz2[1] -= yprd;
  }

  // z coord

  if (zperiodic) {
    while ((xyz2[2] - xyz1[2])*2.0 < zprd)
	 xyz2[2] += zprd;
    while ((xyz2[2] - xyz1[2])*2.0 > zprd)
	 xyz2[2] -= zprd;
  }

}

/* ----------------------------------------------------------------------
   integer shifts for periodic image of xyz2 closest to xyz1
   only works for (neighbor distance)*2 < box length
 ------------------------------------------------------------------------- */

void Domain::set_pbcflags(double* xyz1, double* xyz2, int* pbcflags)
{

  // x coord

  if (xperiodic) {
    if ((xyz2[0] - xyz1[0])*2.0 < -xprd)
	 pbcflags[0] = 1;
    else if ((xyz2[0] - xyz1[0])*2.0 > xprd)
	 pbcflags[0] = -1;
    else
	 pbcflags[0] = 0;
  }

  // y coord

  if (yperiodic) {
    if ((xyz2[1] - xyz1[1])*2.0 < -yprd)
	 pbcflags[1] = 1;
    else if ((xyz2[1] - xyz1[1])*2.0 > yprd)
	 pbcflags[1] = -1;
    else
	 pbcflags[1] = 0;
  }

  // z coord

  if (zperiodic) {
    if ((xyz2[2] - xyz1[2])*2.0 < -zprd)
	 pbcflags[2] = 1;
    else if ((xyz2[2] - xyz1[2])*2.0 > zprd)
	 pbcflags[2] = -1;
    else
	 pbcflags[2] = 0;
  }
}

