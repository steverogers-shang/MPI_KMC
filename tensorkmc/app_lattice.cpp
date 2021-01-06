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

#include "spktype.h"
#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_lattice.h"
#include "comm_lattice.h"
#include "solve.h"
#include "domain.h"
#include "lattice.h"
#include "random_mars.h"
#include "random_park.h"
#include "cluster.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"


using namespace SPPARKS_NS;


#define DELTA 100000  //modified by bdwu 170227 yuanlaiwei 10000

enum{NOSWEEP,RANDOM,RASTER,COLOR,COLOR_STRICT};

/* ---------------------------------------------------------------------- */

AppLattice::AppLattice(SPPARKS *spk, int narg, char **arg) : App(spk,narg,arg)
{
  appclass = LATTICE;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // default settings

  sectorflag = 0;
  nset = 0;
  set = NULL;
  nstop = 1.0;
  tstop = 0.0;

  sweepflag = NOSWEEP;
  ranapp = NULL;
  ranstrict = NULL;
  siteseeds = NULL;
  sitelist = NULL;

  allow_update = 0;

  temperature = 0.0;

  propensity = NULL;
  // i2site = NULL;

  comm = NULL;

  nlocal = nghost = nmax = 0;
  // owner = NULL; delete by smz 170210
  // index = NULL;

  maxneigh = 0;
  numneigh = NULL;
  neighbor = NULL;

  dt_sweep = 0.0;
  naccept = nattempt = 0;
  nsweeps = 0;
  
  update_only = 0;
  
}

/* ---------------------------------------------------------------------- */

AppLattice::~AppLattice()
{
  Solve *s;
  for (int i = 0; i < nset; i++) {
    s = free_set(i);
    delete s;
  }
  delete [] set;

  delete ranapp;
  delete ranstrict;
  memory->destroy(siteseeds);
  memory->destroy(sitelist);

  delete comm;

  // memory->destroy(owner); delete by smz 170210
  // memory->destroy(index);

  memory->destroy(numneigh);
  memory->destroy(neighbor);
}

/* ---------------------------------------------------------------------- */

void AppLattice::input(char *command, int narg, char **arg)
{
  if (strcmp(command,"sector") == 0) set_sector(narg,arg);
  else if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else input_app(command,narg,arg);
}

/* ---------------------------------------------------------------------- */

void AppLattice::init()
{
  // error checks deleted by bdwu
  
  // if sectors, set number of sectors
  if (nprocs > 1 && sectorflag == 0 && solve)
    error->all(FLERR,"Cannot use KMC solver in parallel with no sectors");

  int dimension = domain->dimension;
  nsector = 1;
  if (sectorflag) {

  nsector = 8; } // dimension=3
  //deleted by bdwu 170228

  // if coloring, determine number of colors
  // setup test for create_set
  // check periodicity against lattice extent

  // ncolors = 1;
  // bothflag=0;
  //deleted by bdwu 170228
  

  // create sets based on sectors and coloring
  // set are either all sectors or all colors or both
  // for both, first sets are entire sections, remaining are colors in sectors
  // if new nset is same as old nset, pass each set's solver to create_set,
  //   so it can reuse solver and its RNG,
  //   so consecutive runs are identical to one continuous run

  
  int nsetold = nset;
  Solve **sold = new Solve*[nsetold];

  for (int i = 0; i < nset; i++) sold[i] = free_set(i);
  delete [] set;
  if(nsector > 1){
    nset = nsector;
    set = new Set[nset];
    for (int i = 0; i < nset; i++) {
      if (nset == nsetold) create_set(i,i+1,sold[i]);
      else create_set(i,i+1,NULL);
    }
  }
  else{
    nset = 1;
    set = new Set[nset];
    if (nset == nsetold) create_set(0,0,sold[0]);
    else create_set(0,0,NULL);
  } 

  if (nset != nsetold)
    for (int i = 0; i < nsetold; i++) delete sold[i];
  delete [] sold;
  
 
  // initialize mask array
  //deleted by bdwu 170228

  // setup RN generators, only on first init
  // ranapp is used for all options except sweep color/strict
  // setup ranapp so different on every proc
  // if color/strict, initialize per-lattice site seeds

  if (ranapp == NULL) {
    ranapp = new RandomPark(ranmaster->uniform());
    double seed = ranmaster->uniform();
    ranapp->reset(seed,me,100);
  }

  delete ranstrict;
  memory->destroy(siteseeds);
  ranstrict = NULL;
  siteseeds = NULL;

  // initialize comm, both for this proc's full domain and sectors
  // recall comm->init in case sectoring has changed
  if (me == 0) {
    fprintf(screen, "AppLattice: begin comm->init...");
    fprintf(logfile, "AppLattice: begin comm->init...");
    timer->current();
  }
  if (comm == NULL) comm = new CommLattice(spk);
  comm->init(nsector,delpropensity,delevent,NULL);
  if (me == 0) {
    fprintf(screen, "AppLattice: after comm->init...");
    fprintf(logfile, "AppLattice: after comm->init...");
    timer->current();
  }

  // set sweep function ptr

  // app-specific initialization, after general initialization
  init_app();

  // initialize output

  output->init(time);
}

/* ---------------------------------------------------------------------- */

void AppLattice::setup()
{
  // app-specific setup, before propensities are computed
  if (domain->me == 0) {
    fprintf(screen,"AppLattice::setup app...");
    fprintf(logfile,"AppLattice::setup app...");
    timer->current();
  }
  setup_app();
  // initialize propensities for KMC solver within each set
  // comm insures ghost sites are up to date
  if (domain->me == 0) {
    fprintf(screen,"AppLattice::after setup app...");
    fprintf(logfile,"AppLattice::after setup app...");
    timer->current();
  }
  if (solve) {
    // communicate ev er es
    // comm->all();
    comm->perform_all();


    // int global_vac = 0, proc_vac = 0;
    // int *proc_vacs = new int[nprocs];

    for (int i = 0; i < nset; i++) {
      for(int p = 0; p < set[i].nlocal; p++) {
        set[i].propensity[p] = 0.0;
      }
      count_vac(set[i].nlocal, set[i].site2i, &set[i].vacant);

      std::map<int ,int>::iterator l_it;
      int k = 0;
      for ( l_it = set[i].vacant.begin( ); l_it != set[i].vacant.end( ); l_it++ ) {
        set[i].propensity[k++] = site_propensity(l_it->second);
       //shanghui output:
       //	printf("shanghui test in app_lattice (setup)%d \n ", l_it->second);
      }
      set[i].solve->init(set[i].nlocal, set[i].propensity);
    }

  }

  // convert per-sector time increment info to KMC params

  if (sectorflag && solve) {
    if (tstop > 0.0) {
      Ladapt = false;
      dt_kmc = tstop;
    }

    if (nstop > 0.0) {
      Ladapt = true;
      double pmax = 0.0;
      for (int i = 0; i < nset; i++) {
      	int ntmp = set[i].solve->get_num_active();
      	if (ntmp > 0) {
      	  double ptmp = set[i].solve->get_total_propensity();
      	  ptmp /= ntmp;
      	  pmax = MAX(ptmp,pmax);
      	}
      }
      double pmaxall;
      MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_MAX,world);
      if (pmaxall > 0.0) dt_kmc = nstop/pmaxall;
      else dt_kmc = stoptime-time;
    }

    dt_kmc = MIN(dt_kmc,stoptime-time);
  }
  //deleted by bdwu 170228
  
  setup_end_app();  // second stage of app-specific setup

  nextoutput = output->setup(time); // setup future output
  
}

/* ---------------------------------------------------------------------- */
int kk_icount=0;
void AppLattice::iterate()
{
  timer->barrier_start(TIME_LOOP);
  if (solve) {
    if (sectorflag == 0) 
      iterate_kmc_global(stoptime);
    else
      iterate_kmc_sector(stoptime);    //add by bdwu 160703
  }
  kk_icount++;
  //printf("kk iterate %d\n",kk_icount);

  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   KMC solver on entire domain
   can only be invoked in serial
 ------------------------------------------------------------------------- */

void AppLattice::iterate_kmc_global(double stoptime)
{
  int isite;
  //int gloablcc = 0;
  
  // global KMC runs with one set
  // save ptr to system solver

  Solve *hold_solve = solve;
  solve = set[0].solve;
  propensity = set[0].propensity;

  vacant = &(set[0].vacant);
  //border = set[iset].border;
  //nborder = set[iset].nborder;
  current_set = 0;

  int done = 0;
  int notinset = 0;
  while (!done) {
    //gloablcc++;
    //printf("KKK run %d !\n",gloablcc);
    timer->stamp();
    isite = solve->event(&dt_step);
    timer->stamp(TIME_SOLVE);

    if (isite >= 0) {
      time += dt_step;
      if (time <= stoptime) {
	    site_event(isite,&notinset,ranapp);
	    naccept++;
	    timer->stamp(TIME_APP);
      } else {
      done = 1;
      time = stoptime;
      }
    } else {
      done = 1;
      time = stoptime;
    }

    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }

  // restore system solver

  solve = hold_solve;
}

/* ----------------------------------------------------------------------
   KMC solver on sectors
   can be invoked in serial or parallel
 ------------------------------------------------------------------------- */

void AppLattice::iterate_kmc_sector(double stoptime)
{
  if (me == 0) timer->current();
  int i,isite,done,nsites, *border, nborder;
  double dt,timesector;
  double pmax,pmaxall;
  std::map <int, int>::iterator Iter;
  // save ptr to system solver

  Solve *hold_solve = solve;
  
  xmid = 0.5 * (domain->subxlo + domain->subxhi);
  ymid = 0.5 * (domain->subylo + domain->subyhi);
  zmid = 0.5 * (domain->subzlo + domain->subzhi);
  
  int alldone = 0;
  int kk_allcount=0;
  int myid;//kk
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);//kk

  while (!alldone) {
    int notinset = 0;
    if (Ladapt) pmax = 0.0;
    for (int iset = 0; iset < nset; iset++) {
      timer->stamp();
      comm->set_iset(iset);

      solve = set[iset].solve;
      propensity = set[iset].propensity;
      // i2site = set[iset].i2site;
      vacant = &(set[iset].vacant);
      border = set[iset].border;
      nborder = set[iset].nborder;
      current_set = iset;

      // update vacancy number and propensity
      // necessary since outside sites may have changed
      nsites = 0;
      for (int m = 0; m < nborder; m++) {
		    i = border[m];
		    if(type(i) == 0) {
          int flag = 0, num = vacant->size();
		      for ( Iter = vacant->begin( ); Iter != vacant->end( ); Iter++ ) {
				    if(i == Iter->second) {flag = 1; break;}
          }
          if(flag == 0)  vacant->insert(std::pair<int, int>(num, i)); 
		    }
      }

      int *bsites = new int[vacant->size()];  

      for ( Iter = vacant->begin( ); Iter != vacant->end( ); Iter++ ) {
	      i = Iter->first;
        int site = Iter->second;
        propensity[i] = site_propensity(site);
        //shanghui output:
	//printf("shanghui test in app_lattice(iterate_kmc_sector) %d \n", site);
        bsites[nsites++] = i;
      }

 
      solve->update(nsites,bsites,propensity);
      timer->stamp(TIME_UPDATE);
      // pmax = maximum sector propensity per site
      if (Ladapt) {
        int ntmp = solve->get_num_active();
        if (ntmp > 0) {
          double ptmp = solve->get_total_propensity();
          ptmp /= ntmp;
          pmax = MAX(ptmp,pmax);
        }
      }
      // execute events until sector time threshhold reached
      done = 0;
      timesector = 0.0;
      int kk_seccount =0;
      while (!done) {
      kk_seccount++;
      //if(myid==0)
      //printf("kk!!! kk_allcount:%d kk_seccount:%d\n", kk_allcount,kk_seccount); 
		    timer->stamp();

		    isite = solve->event(&dt); 
        //if(myid==0)
        //printf("isite %d timesector %.3f, dt %.3f, dt_kmc %.3f\n", isite,timesector,dt,dt_kmc);
              
		    timer->stamp(TIME_SOLVE);
    		if (isite < 0) done = 1;
    		else {
    		  timesector += dt;
          //printf("timesector %.3f, dt %.3f, dt_kmc %.3f\n", timesector,dt,dt_kmc);
    		  if (timesector >= dt_kmc) done = 1; 
    		  else {
    		    site_event(isite,&notinset,ranapp);
    		    naccept++;
    		  }
    		  timer->stamp(TIME_APP);
    		}
            //if(me == 0){
               // fprintf(screen, " kk_allcount %d iset %d nset %d kk_seccount %d AppLattice: in the loop %lf  %lf %d %d\n", kk_allcount,iset,nset, kk_seccount,timesector, dt_kmc, isite, done);
            //    fprintf(logfile, "AppLattice: in the loop %lf  %lf %d %d\n", timesector, dt_kmc, isite, done);
            //}
      }
      delete[] bsites;
      if (nprocs > 1) {
    		communicate();
    		timer->stamp(TIME_COMM);
      }
    }

    if (allow_update) user_update(dt_kmc);

    // keep looping until overall time threshhold reached
 
    nsweeps++;
    time += dt_kmc;
    //alldone = 1;
    if (time >= stoptime) alldone = 1;
    kk_allcount++;
    //printf("kk alldone:%d\n", alldone);
    if (alldone || time >= nextoutput)
      nextoutput = output->compute(time,alldone);
    timer->stamp(TIME_OUTPUT);

    //if(me == 0){
    //    fprintf(screen, "AppLattice: in the loop %lf %lf %lf\n", time, dt_kmc, stoptime);
    //    fprintf(logfile, "AppLattice: in the loop %lf %lf %lf\n", time, dt_kmc, stoptime);
    //}

    // recompute dt_kmc if adaptive, based on pmax across all sectors

    if (Ladapt) {
      MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_MAX,world);
      if (pmaxall > 0.0) dt_kmc = nstop/pmaxall;
      else dt_kmc = stoptime-time;
      dt_kmc = MIN(dt_kmc,stoptime-time);
    }
  }

  // restore system solver
  solve = hold_solve;
}

/* ----------------------------------------------------------------------
   rejection KMC solver
 ------------------------------------------------------------------------- */
//delete by bdwu 160703


/* ----------------------------------------------------------------------
   iterate the user_update routine only
   app is responsible for doing communciation in user_update()
------------------------------------------------------------------------- */
//delete by bdwu 160703


/* ---------------------------------------------------------------------- */
//deleted by bdwu 170228
// void AppLattice::sweep_nomask_nostrict(int n, int *list)
// void AppLattice::sweep_mask_nostrict(int n, int *list)
// void AppLattice::sweep_nomask_strict(int n, int *list)
// void AppLattice::sweep_mask_strict(int n, int *list)

/* ---------------------------------------------------------------------- */

void AppLattice::input_app(char *command, int narg, char **arg)
{
  error->all(FLERR,"Unrecognized command");
}


/* ---------------------------------------------------------------------- */

void AppLattice::set_sector(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal sector command");

  nsector_user = 0;
  if (strcmp(arg[0],"yes") == 0) sectorflag = 1;
  else if (strcmp(arg[0],"no") == 0) sectorflag = 0;
  else {
    sectorflag = 1;
    nsector_user = atoi(arg[0]);
    if (nsector_user != 2 && nsector_user != 4 && nsector_user != 8)
      error->all(FLERR,"Illegal sector command");
  }

  nstop = 1.0;
  tstop = 0.0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"nstop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal sector command");
      nstop = atof(arg[iarg+1]);
      if (nstop <= 0.0) error->all(FLERR,"Illegal sector command");
      tstop = 0.0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"tstop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal sector command");
      tstop = atof(arg[iarg+1]);
      if (tstop <= 0.0) error->all(FLERR,"Illegal sector command");
      nstop = 0.0;
      iarg += 2;
    } else error->all(FLERR,"Illegal sector command");
  }
}
/* ---------------------------------------------------------------------- */
//delete by bdwu 170228
//void AppLattice::set_sweep(int narg, char **arg)

/* ---------------------------------------------------------------------- */

void AppLattice::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal temperature command");
  temperature = atof(arg[0]);
  if (temperature != 0.0) t_inverse = 1.0/temperature;
}

/* ---------------------------------------------------------------------- */
//delete by bdwu 170228
//void AppLattice::set_update_only(int narg, char **arg)


/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppLattice::stats(char *strtmp)
{
  char big[8],format[64];
  strcpy(big,BIGINT_FORMAT);

  bigint naccept_all;
  MPI_Allreduce(&naccept,&naccept_all,1,MPI_SPK_BIGINT,MPI_SUM,world);

  if (solve) {
    sprintf(format,"%%10g %%10%s %%10d %%10d",&big[1]);
    sprintf(strtmp,format,time,naccept_all,0,nsweeps);
  } else {
    bigint nattempt_all;
    MPI_Allreduce(&nattempt,&nattempt_all,1,MPI_SPK_BIGINT,MPI_SUM,world);
    sprintf(format,"%%10g %%10%s %%10%s %%10d",&big[1],&big[1]);
    sprintf(strtmp,format,time,naccept_all,nattempt_all-naccept_all,nsweeps);
  }
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppLattice::stats_header(char *strtmp)
{
  sprintf(strtmp,"%10s %10s %10s %10s","Time","Naccept","Nreject","Nsweeps");
}

/* ----------------------------------------------------------------------
   create a subset of owned sites
   insure all ptrs in Set data struct are allocated or NULL
   isector = 0 = all sites (no sector)
   isector > 1 = sites within a sector
   icolor = 0 = all sites (no color)
   icolor > 1 = sites of a certain color
 ------------------------------------------------------------------------- */

void AppLattice::create_set(int iset, int isector, Solve *oldsolve)
{
  double unitdis = 2 / domain->lattice->getLatconst() ;
  int subx = round((domain->subxhi - domain->subxlo) * unitdis);
  int suby = round((domain->subyhi - domain->subylo) * unitdis);
  int subz = round((domain->subzhi - domain->subzlo)  * unitdis);
  if(domain->me == 0){
    printf("kk ......... %d %d %d\n",subx,suby,subz );
  }
  
  if (subx< 12 || suby < 12 || subz < 12)
      error->all(FLERR,"sector side can't be less than 6!");
  // sector boundaries

  double xmid = 0.5 * (domain->subxlo + domain->subxhi);
  double ymid = 0.5 * (domain->subylo + domain->subyhi);
  double zmid = 0.5 * (domain->subzlo + domain->subzhi);

  // count sites in subset

  int flag,iwhich,jwhich,kwhich,msector,mcolor;

  int delcolor = delevent + delpropensity;

  int n = 0;
  for (int i = 0; i < nlocal; i++) {
    flag = 1;

    if (isector > 0) {
      if (xyz[i][0] < xmid) iwhich = 0;
      else iwhich = 1;
      if (xyz[i][1] < ymid) jwhich = 0;
      else jwhich = 1;
      if (xyz[i][2] < zmid) kwhich = 0;
      else kwhich = 1;

      if (nsector == 8)
        msector = 4*kwhich + 2*jwhich + iwhich + 1;
      if (isector != msector) flag = 0;
    }

    if (flag) n++;
  }

  set[iset].nlocal = n;

  // setup site2i for sites in set

  memory->create(set[iset].site2i,n,"app:site2i");

  n = 0;
  for (int i = 0; i < nlocal; i++) {
    flag = 1;

    if (isector > 0) {
      if (xyz[i][0] < xmid) iwhich = 0;
      else iwhich = 1;
      if (xyz[i][1] < ymid) jwhich = 0;
      else jwhich = 1;
      if (xyz[i][2] < zmid) kwhich = 0;
      else kwhich = 1;

      if (nsector == 8)
        msector = 4*kwhich + 2*jwhich + iwhich + 1;

      if (isector != msector) flag = 0;
    }

    if (flag) set[iset].site2i[n++] = i;
  }

  // setup i2site for sites in set, only for KMC solver
  // i2site = 0 to nsite-1 for owned points in set, else -1

  // delete by smz 170302. don't need i2site
  // if (solve) {
  //   (int *) memory->create(set[iset].i2site,nlocal+nghost,"app:i2site");
  //   for (int i = 0; i < nlocal+nghost; i++) set[iset].i2site[i] = -1;
  //   for (int i = 0; i < set[iset].nlocal; i++) 
  //     set[iset].i2site[set[iset].site2i[i]] = i;
  // } else set[iset].i2site = NULL;

  // allocate propensity array for set

  memory->create(set[iset].propensity,set[iset].nlocal,"app:propensity");

  // allocate KMC solver for set
  // reuse old solver if provided, else delete it

  if (solve) {
    if (oldsolve) set[iset].solve = oldsolve;
    else set[iset].solve = solve->clone();
  } else {
    set[iset].solve = NULL;
    delete oldsolve;
  }

}

/* ----------------------------------------------------------------------
   free memory inside a set
   return Solver so caller can reuse it if desired
 ------------------------------------------------------------------------- */

Solve *AppLattice::free_set(int iset)
{
  memory->destroy(set[iset].border);
  memory->destroy(set[iset].propensity);
  memory->destroy(set[iset].site2i);
  // memory->destroy(set[iset].i2site);
  set[iset].vacant.clear();
  return set[iset].solve;
}

/* ----------------------------------------------------------------------
   create list of border sites for a set
   border site = site in set with a 1 to Nlayer neighbor outside the set
   neighbor can be another owned site (outside set) or a ghost
   border = lattice index of the sites
 ------------------------------------------------------------------------- */
// used in create_set(), not used deleted by bdwu 170228
//int AppLattice::find_border_sites(int isector)

  
/* ----------------------------------------------------------------------
   unset all mask values of owned sites in iset whose propensity
     could change due to events on sites one neighbor outside the set
   border list stores indices of these sites
   their mask value may be out of date, due to state change in other sets
 ------------------------------------------------------------------------- */
// not used  deleted by bdwu 170228
//void AppLattice::boundary_clear_mask(int iset)


/* ----------------------------------------------------------------------
   push new site onto stack and assign new id
 ------------------------------------------------------------------------- */
// used in diag_cluster.cpp
void AppLattice::push_new_site(int i, int* cluster_ids, int id,
					  std::stack<int>* cluststack)
{
  // This can be used to screen out unwanted spin values
  // int isite = iarray[0][i];

  cluststack->push(i);
  cluster_ids[i] = id;
}

/* ----------------------------------------------------------------------
   push connected neighbors of this site onto stack
     and assign current id
   ghost neighbors are masked by id = -1
   previously burned sites are masked by id > 0
 ------------------------------------------------------------------------- */
// used in diag_cluster.cpp
void AppLattice::push_connected_neighbors(int i, int* cluster_ids, int id,
					  std::stack<int>* cluststack)
{
  int ii;
  int isite = iarray[0][i];

  for (int j = 0; j < numneigh[i]; j++) {
    ii = neighbor[i][j];
    if (iarray[0][ii] == isite && cluster_ids[ii] == 0) {
      cluststack->push(ii);
      cluster_ids[ii] = id;
      domain->pbcshift(xyz[i],xyz[ii]);
    }
  }
}

/* ----------------------------------------------------------------------
   add cluster id of connected ghost sites to neighbor list of cluster
 ------------------------------------------------------------------------- */
//used in diag_cluster
void AppLattice::connected_ghosts(int i, int* cluster_ids, 
				  Cluster* clustlist, int idoffset)
{
  int iclust;
  int ii;
  int isite = iarray[0][i];
  int pbcflags[3];

  // check if this was a site that was ignored

  if (cluster_ids[i] == 0) return;

  iclust = cluster_ids[i]-idoffset;
  
  // add ghost cluster to neighbors of local cluster

  for (int j = 0; j < numneigh[i]; j++) {
    ii = neighbor[i][j];
    if (iarray[0][ii] == isite && ii >= nlocal) {
      clustlist[iclust].add_neigh(cluster_ids[ii]);
      domain->set_pbcflags(xyz[i],xyz[ii],pbcflags);
      clustlist[iclust].add_pbcflags(cluster_ids[ii],pbcflags);
    }
  }
}

/* ----------------------------------------------------------------------
   grow per-site arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n 
------------------------------------------------------------------------- */

void AppLattice::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  memory->grow(id,nmax,"app:id");
  memory->grow(xyz,nmax,3,"app:xyz");
  // memory->grow(owner,nmax,"app:owner"); delete by smz 0210
  // memory->grow(index,nmax,"app:index");

  memory->grow(numneigh,nmax,"app:numneigh");
  // if (maxneigh) memory->grow(neighbor,nmax,maxneigh,"app:neighbor"); //delete by smz 161107

  for (int i = 0; i < ninteger; i++)
    memory->grow(iarray[i],nmax,"app:iarray");
  for (int i = 0; i < ndouble; i++)
    memory->grow(darray[i],nmax,"app:darray");

  grow_app();
}

/* ----------------------------------------------------------------------
   add an owned site
   called from create_sites command
   grow arrays if necessary
 ------------------------------------------------------------------------- */
// used in create_sites.cpp read_sites.cpp
void AppLattice::add_site(tagint n, double x, double y, double z)
{
  if (nlocal == nmax) grow(0);

  id[nlocal] = n;
  xyz[nlocal][0] = x;
  xyz[nlocal][1] = y;
  xyz[nlocal][2] = z;

  

  // owner[nlocal] = me; delete by smz 170210
  // index[nlocal] = nlocal;

  for (int i = 0; i < ninteger; i++) iarray[i][nlocal] = 0;
  for (int i = 0; i < ndouble; i++) darray[i][nlocal] = 0.0;

  nlocal++;
}

/* ----------------------------------------------------------------------
   add a ghost site
   called from create_sites command
   grow arrays if necessary
 ------------------------------------------------------------------------- */
// used in ghosts_from_connectivity
void AppLattice::add_ghost(tagint n, double x, double y, double z)
{
  if (nlocal+nghost == nmax) grow(0);

  int m = nlocal + nghost;

  id[m] = n;
  xyz[m][0] = x;
  xyz[m][1] = y;
  xyz[m][2] = z;

  // owner[m] = proc; delete by smz 170210
  // index[m] = index_owner;

  for (int i = 0; i < ninteger; i++) iarray[i][m] = 0;
  for (int i = 0; i < ndouble; i++) darray[i][m] = 0.0;

  nghost++;
}

/* ----------------------------------------------------------------------
   set neighbor connectivity for owned site I
   nvalues = # of neighbors
   called from read_sites command
 ------------------------------------------------------------------------- */
//used in read_sites.cpp
void AppLattice::add_neighbors(int i, int nvalues, char **values)
{
  numneigh[i] = nvalues;
  for (int m = 0; m < nvalues; m++)
    neighbor[i][m] = atoi(values[m]);
}

/* ----------------------------------------------------------------------
   set values for owned site I
   called from read_sites command
 ------------------------------------------------------------------------- */
//used in read_sites.cpp
void AppLattice::add_values(int i, char **values)
{
  for (int m = 0; m < ninteger; m++) iarray[m][i] = atoi(values[m]);
  for (int m = 0; m < ndouble; m++) darray[m][i] = atof(values[m+ninteger]);
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   lo,hi = inclusive bounds
   5 possibilities:
     (1) i = i to i, (2) * = lo to hi,
     (3) i* = i to hi, (4) *j = lo to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */
//not used
void AppLattice::add_value(int i, int type, int index, char *value)
{
  if (type == 0) iarray[index-1][i] = atoi(value);
  else if (type == 1) darray[index-1][i] = atof(value);
}
void AppLattice::bounds(char *str, int lo, int hi, int &nlo, int &nhi)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = MAX(atoi(str),lo);
    nhi = MIN(atoi(str),hi);
  } else if (strlen(str) == 1) {
    nlo = lo;
    nhi = hi;
  } else if (ptr == str) {
    nlo = lo;
    nhi = MIN(atoi(ptr+1),hi);
  } else if (strlen(ptr+1) == 0) {
    nlo = MAX(atoi(str),lo);
    nhi = hi;
  } else {
    nlo = MAX(atoi(str),lo);
    nhi = MIN(atoi(ptr+1),hi);
  }
}

/* ----------------------------------------------------------------------
   print connectivity stats
 ------------------------------------------------------------------------- */
//used in read_sites.cp
void AppLattice::print_connectivity()
{
  int i;

  tagint min = maxneigh;
  tagint max = 0;

  for (i = 0; i < nlocal; i++) {
    min = MIN(min,numneigh[i]);
    max = MAX(max,numneigh[i]);
  }

  tagint minall,maxall;
  MPI_Allreduce(&min,&minall,1,MPI_SPK_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&max,&maxall,1,MPI_SPK_TAGINT,MPI_MAX,world);

  tagint *count = new tagint[maxall+1];
  tagint *countall = new tagint[maxall+1];

  for (i = 0; i <= maxall; i++) count[i] = 0;

  for (i = 0; i < nlocal; i++) count[numneigh[i]]++;
  MPI_Allreduce(count,countall,maxall+1,MPI_SPK_TAGINT,MPI_SUM,world);

  if (me == 0)
    for (i = minall; i <= maxall; i++) {
      if (screen)
	fprintf(screen,"  " TAGINT_FORMAT " sites have %d neighbors\n",\
		countall[i],i);
      if (logfile)
	fprintf(logfile,"  " TAGINT_FORMAT " sites have %d neighbors\n",
		countall[i],i);
    }

  delete [] count;
  delete [] countall;
}
