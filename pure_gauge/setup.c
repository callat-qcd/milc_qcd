/******** setup.c *********/
/* MIMD version 7 */
#define IF_OK if(status==0)


#include <string.h>
#include <ctype.h>

/* Forward declarations */
#ifdef ANISOTROPY
static int dirstring2index (char savebuf[], int *status);
#endif
#include "pure_gauge_includes.h"
int initial_set();

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int  setup()   {
int prompt;

        /* print banner, get volume, nflavors, seed */
    prompt=initial_set();
        /* initialize the node random number generator */
    initialize_prn(&node_prn,iseed,volume+mynode());
        /* Initialize the layout functions, which decide where sites live */
    setup_layout();
        /* allocate space for lattice, set up coordinate fields */
    make_lattice();
        /* set up neighbor pointers and comlink structures */
    make_nn_gathers();

    return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
int prompt,status;
    /* On node zero, read lattice size, seed, nflavors and send to others */
    if(mynode()==0){
        /* print banner */
#ifndef ANISOTROPY
        printf("Pure gauge SU3\n");
#else
        printf("Anisotropic pure gauge SU3\n");
#endif
#ifdef RMD_ALGORITHM
        printf("Microcanonical simulation with refreshing\n");
#endif
#ifdef HMC_ALGORITHM
        printf("Microcanonical simulation with refreshing\n");
#endif
        printf("MIMD version 7\n");
        printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
#ifdef HMC_ALGORITHM
        printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef ORA_ALGORITHM
        printf("Overrelaxed/quasi-heat bath algorithm\n");
#endif
	time_stamp("start");
        status=get_prompt(stdin, &prompt);
	IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
	IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
	IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
	IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );
	IF_OK status += get_i(stdin, prompt,"iseed", &par_buf.iseed );

	if(status>0) par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(mynode()==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

if( par_buf.stopflag != 0 )
  normal_exit(0);


    nx=par_buf.nx;
    ny=par_buf.ny;
    nz=par_buf.nz;
    nt=par_buf.nt;
    iseed=par_buf.iseed;
    
    this_node = mynode();
    number_of_nodes = numnodes();
    volume=(size_t)nx*ny*nz*nt;
    return(prompt);
}

/* read in parameters and coupling constants    */
int readin(int prompt) {
/* read in parameters for su3 monte carlo       */
/* argument "prompt" is 1 if prompts are to be given for input */

int status;
#ifdef ORA_ALGORITHM
 int status2;
char savebuf[128];
#endif

    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

        printf("\n\n");
	status=0;
    
	/* warms, trajecs */
	IF_OK status += get_i(stdin, prompt,"warms", &par_buf.warms );
	IF_OK status += get_i(stdin, prompt,"trajecs", &par_buf.trajecs );
    
	/* trajectories between propagator measurements */
	IF_OK status += 
	    get_i(stdin, prompt,"traj_between_meas", &par_buf.propinterval );
    
	/* get couplings and broadcast to nodes	*/
	/* beta */
#ifndef ANISOTROPY
	IF_OK status += get_f(stdin, prompt,"beta", &par_buf.beta );
#else
        /* beta[0] - isotropic, beta[1] - anisotropic */
        IF_OK status += get_vf(stdin, prompt,"beta", par_buf.beta, 2 );
        /* Direction of anisotropy */
        IF_OK status += get_s(stdin, prompt,"ani_dir",savebuf);
        par_buf.ani_dir = dirstring2index( savebuf, &status);
#endif

#if ( defined HMC_ALGORITHM || defined RMD_ALGORITHM )
        /* microcanonical time step */
	IF_OK status +=
            get_f(stdin, prompt,"microcanonical_time_step", &par_buf.epsilon );
#endif  
        /*microcanonical steps per trajectory */
	IF_OK status += get_i(stdin, prompt,"steps_per_trajectory", &par_buf.steps );
    
#ifdef ORA_ALGORITHM
        /*qhb steps per trajectory */
	IF_OK status += get_i(stdin, prompt,"qhb_steps", &par_buf.stepsQ );
#endif   

        /* find out what kind of starting lattice to use */
	IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
	    par_buf.startfile );

 
#ifdef ORA_ALGORITHM
       IF_OK if (prompt==1) printf(
           "enter 'no_gauge_fix', 'landau_gauge_fix', or 'coulomb_gauge_fix'\n");
       IF_OK status2=scanf("%s",savebuf);
       IF_OK {
	 if(strcmp("coulomb_gauge_fix",savebuf) == 0 ){
	   par_buf.fixflag = COULOMB_GAUGE_FIX;
	   if(this_node==0)printf("fixing to coulomb gauge\n");
	 }
	 else if(strcmp("landau_gauge_fix",savebuf) == 0 ) {
	   par_buf.fixflag = LANDAU_GAUGE_FIX;
	   if(this_node==0)printf("fixing to landau gauge\n");
	 }
	 else if(strcmp("no_gauge_fix",savebuf) == 0 ) {
	   par_buf.fixflag = NO_GAUGE_FIX;
	   if(this_node==0)printf("NOT fixing the gauge\n");
	 }
	 else{
           printf("error in input: fixing_command %s is invalid\n",savebuf); 
	   status++;
	 }
       }
#endif
 
        /* find out what to do with lattice at end */
	IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
	    par_buf.savefile );
	IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				      par_buf.stringLFN );
 
        /* send parameter structure */
	if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
    } /* end if(this_node==0) */

    /* Node 0 broadcasts parameter buffer to all other nodes */
    broadcast_bytes((char *)&par_buf,sizeof(par_buf));

    if( par_buf.stopflag != 0 )
      normal_exit(0);

    warms = par_buf.warms;
    trajecs = par_buf.trajecs;
    steps = par_buf.steps;
    stepsQ = par_buf.stepsQ;
    propinterval = par_buf.propinterval;
    startflag = par_buf.startflag;
    fixflag = par_buf.fixflag;
    saveflag = par_buf.saveflag;
    epsilon = par_buf.epsilon;
#ifndef ANISOTROPY
    beta = par_buf.beta;
#else
    ani_dir=par_buf.ani_dir;

    beta[0] = par_buf.beta[0];
    beta[1] = par_buf.beta[1];
#endif
    strcpy(startfile,par_buf.startfile);
    strcpy(savefile,par_buf.savefile);
    strcpy(stringLFN, par_buf.stringLFN);

    /* Do whatever is needed to get lattice */
    if( startflag != CONTINUE )
      startlat_p = reload_lattice( startflag, startfile );

    return(0);
} /*readin()*/

#ifdef ANISOTROPY
static int dirstring2index (char savebuf[], int *status) {
  short mydir;

  if ( savebuf[0] >= 'A' && savebuf[0] <= 'Z' )
    mydir = tolower(savebuf[0]);
  else
    mydir = savebuf[0];
  switch(mydir) {
    case XUP:
    case '0':
    case 'x': mydir = XUP; break;
    case YUP:
    case '1':
    case 'y': mydir = YUP; break;
    case ZUP:
    case '2':
    case 'z': mydir = ZUP; break;
    case TUP:
    case '3':
    case 't': mydir = TUP; break;
    default:
      node0_printf("Expecting direction \
as x,y,z,t, X,Y,Z,T, or 0,1,2,3;  instead %c\n", mydir);
     (*status)++;
  }
  return mydir;
}
#endif
