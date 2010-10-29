/* A U T O C O R R E L A T O R Version 1.0.0, 26.07.2010 by PK */

/* Autocorrelator code for msd */

/**

The program computes autocorrelation function using the blocking 
algorithm described in Frenkel and Smit: Understanding molecular
simulation, Academic Press, Sand Diego, 2002

It can be found in Chapter 4.4.2: Order-n algorithm to measure 
correlations, which is page 90 in the issue from 2002

The algorithm as presented in the book can be directly used 
to analyse correlation functions which are linearly additive, 
such as scalar product, vectors, etc.

In a modified version, it can also be used for non-linearly
additive quantities such as mean-square displacement (msd).

In the current version, the correlator algorithm is used to compute
msd. In future, instead of functions operation() and dump_res(),
pointer to functions should be used to enable computation of
an arbitrary quantity. The correction for the systematic error
in the case of msd is implemented in the dmp_res function.
It is however by no means guaranteed, that a similar correction
can be used for other non-linearly additive quantities.

PK, 26.07.2010

*/

/* 
TODO

* the correction of the systematic error is not really well readable but it works !!!

*/

// in case the program should only be used as a library, undefine TEST
#include "autocorr.h"
#ifndef TEST
#define TEST 0
#endif
#ifndef DEBUG 
#define DEBUG 0
#endif

#if TEST
// This is just for testing the correlator as a standalone program and should be removed later on
// global variables
int nc=2;			// groupping number
double dt=1.0;			// dt
long int n_vals=0;			// number of values
int n_part=0;			// number of particles
int gs=8;			// group size
char *inif_name="corr.in";	// input file name
char *outf_name="corr.out";	// output file name
int linear_flag=0;		// if set to zero, do blocking correlation, otherwise do linear correlation

// temporary stuff for debugging
int *n_vals_hist; // histogram to note how many values were added to a particlar level of results
int level_shift; // histogram to note how many values were added to a particlar level of results

void error (char *s) /* error messages */
{
        fprintf (stderr, "\ncorrelator: %s\n\t the program was terminated\n\n", s);
        fflush (stdout);
        exit (1);
}

// read the data in ascii format -- for testing only
double **read_data(char *inif_name, int n_vals, int n_part){
	FILE *f=fopen(inif_name, "r");
	if (f==NULL) error("cannot open inif\n");
	int val=0, part=0;
	double **d;
	d=(double**)malloc(n_vals*sizeof(double*));
	for(val=0;val<n_vals;val++) d[val]=(double*)malloc(3*n_part*sizeof(double)); 
	for(val=0;val<n_vals;val++) {
	  for(part=0;part<n_part;part++) {
	    fscanf(f,"%lf %lf %lf", &d[val][3*part],&d[val][3*part+1],&d[val][3*part+2]);
	  }
	}
	fclose(f);
	return d;
}

void print_help()
{
	printf ("\nusage: corr [options]\n \
	options:\n \
	-f <filename>    input file (default: corr.in)\n \
	-o <filename>    output file (default: corr.out)\n \
	-n <integer>     number of stored configs\n \
	-p <integer>     number of particles stored in configs\n \
	-b <integer>     number of data to be blocked in one level (default: 2)\n \
	-g <integer>     group size (default: 8)\n \
	-t <double>      time step (channel width), (default: 1.0)\n \
	-A               read data (times of the photon arrivals) as a single ASCII column rather than t3r2ascii output\n \
	-l               use linear version of the algorithm instead of the blocking one (for testing of the results)\n \
	-h               help\n\n \
	Expected input data format: the output of t3r2ascii without the header (can be changed with -A option)\n\n");
	exit(1);
}

int parse_options (int argc, char *argv[])
{
	if (argc==1) print_help();
	while (--argc && (*++argv)[0] == '-') 
	{
		char c;
		c = *++argv[0];
		switch (c) 
		{
			case 'n' :
				n_vals=atoi(*++argv);
				break;
			case 'p' :
				n_part=atoi(*++argv);
				break;
			case 'b' :
				nc=atoi(*++argv);
				break;
			case 'g' :
				gs=atoi(*++argv);
				break;
			case 't' :
				dt=atof(*++argv);
				break;
			case 'f':
				inif_name=(*++argv);
				break;	
			case 'o':
				outf_name=(*++argv);
				break;	
			case 'l':
				linear_flag=1; 
				++argc;
				break;	
			case 'h' :
				print_help();
				break;
			default:
				printf ("\ncorr: illegal option %s\ntry corr -h for help\
					\n\n", argv[0]);
				exit (1);
		
		}
		--argc;
	}
	return 0;
}

int get_n_vals (char *inif_name) {
	int i=0;
	FILE *inif=fopen(inif_name, "r");
	if (inif==NULL) error("cannot open inif\n");
	char c=fgetc(inif);
	while (c!=EOF) {
		c=getc(inif);
		if(c=='\n') i++;
	}
	fclose(inif);
	return i;
}

// print the results to a file
int dump_res(FILE *f, double *result, int n_res, int nc, int gs, double dt, int n_vals, int n_part) {
	int act_n_vals=n_vals*n_part; // number of data values contributing to the average at a particular point
	int i;
	int d=0; // difference in current blocking level
	int shift=nc*gs-1;
	int t_i=0;
	int dt_i=1;
	
	// correcting the systematic error in the result
	int c_index=0; // index in the results field to be used for correction of the systematic error 
	double *fac; // correction factor
	int *fac_t_i; // time corresponding to the correction factor
	int n_fac=n_res/gs-1;
	int fac_index=1; // first factor is 0.0 by default
	int t_i_fac=1;
	fac=(double*)malloc(n_fac*sizeof(double));
	fac_t_i=(int*)malloc(n_fac*sizeof(double));
	for(i=0;i<n_fac;i++) { fac[i]=0.0; fac_t_i[i]=0; }

	fprintf(f,"#time\t product(acf)\n");
	for(i=0; i<n_res; i++){
		// correct for the systematic error due to coarse-graining
		act_n_vals=n_part*(n_vals-d);
		if(t_i_fac==t_i && fac_index<n_fac) {
			fac[fac_index]=0.5*(result[i]/act_n_vals+fac[c_index]);
			fac_t_i[fac_index]=t_i;
			fac_index++;
			t_i_fac*=nc;
		} else if(t_i > t_i_fac && fac_index<n_fac) {
			fprintf(stderr,"autocorr unexpected error: impossible to aply correction for the systematic error, t_i=%lf > t_i_fac=%lf\n",t_i*dt,t_i_fac*dt);
			exit(1);
		}
		#if DEBUG
		printf("adding %lf from t=%lf to the result at t=%lf\n",fac[c_index],fac_t_i[c_index]*dt,t_i*dt);
		#endif
		fprintf(f, "%e\t %e\t %d \t %d\n", t_i*dt, result[i]/act_n_vals+fac[c_index], act_n_vals, n_vals_hist[i]);
		t_i+=dt_i;
		if(i==shift){ // if we have crossed the block border, change dt
		      #if DEBUG
			printf("change c_index from %d",c_index);
		      #endif
			c_index++;
		      #if DEBUG
			printf("to %d\n",c_index);
		      #endif
			dt_i*=nc;
			//shift+=shift;
			shift+=gs;
			n_vals/=nc;
			d=gs;
		      #if DEBUG
			printf("changing shift to %d at t=%f and i=%d\n", shift,t_i*dt,i); fflush(stdout);
		      #endif
		} else { d++; }
	}
	fflush(f);
	return 0;
}

// initialise the array where results will be stored
double* init_res(int n_vals){
	int i;
	double *result;
	result=(double*)malloc(n_vals*sizeof(double));
	for(i=0; i<n_vals; i++) result[i]=0.0;
	//dump_res(result, n_vals);
	return result;
}

#endif

int get_l_max (long int n_vals, int nc, int gs){
	int i=0; 
	n_vals/=gs; // divide it by the basic groupping gs
	while(n_vals/nc){ 
		n_vals/=nc;
		i++;
	}
	if(i==0) return 1; // at least zero-th blocking level is always possible
	return i; 
	// in principle one more level of blocking should be possible 
	// but it will be only incompletely filled with values which may be confusing
}

// operation to be performed on the data
double my_operation (double **data, int pid, int t, int tau) {
  double dx, dy, dz;
  dx=data[t][3*pid  ]-data[t+tau][3*pid  ];
  dy=data[t][3*pid+1]-data[t+tau][3*pid+1];
  dz=data[t][3*pid+2]-data[t+tau][3*pid+2];
  return (dx*dx+dy*dy+dz*dz);
}


int simple_correlate(double *result, int n_res, double **data, int n_part, int n_vals, int t, int tau_min, int tau_max_plus) {
	/* result points to such a position in the results, that corresponds to act_dt=dt
	   dt is the current dt corresponding to the blocking level; 
	   simple_correlate does not bother about the blocking 
	   */

	int tau=tau_min; // actual dt (between the two processed values)
	int pid; // particle id
	double product; // "product" of the operation performed on data (whatever this operation may be in future, e.g. scalar product of displacement squared or a simple product or ...)
	// avoid the rounding error
	while (tau<tau_max_plus && (t+tau)<n_vals) {
		for(pid=0; pid<n_part; pid++) {
			product=my_operation(data, pid, t, tau); // product of my operation
			result[tau]+=product; // add the output to the results
			#if TEST
			n_vals_hist[level_shift+tau]+=1; // write down that we have added the value
			if(tau>(n_res-1)) {
			  printf("Warning: tau: %d > %d n_res\n",tau,n_res); 
			  exit(1);
			}
			#endif
		}
		tau++; // go to the next lag time
	}
	return 0; // currently the return value has no meaning
}

// perform the blocking of data
int block_data(double **data, int nc, int n_vals, int n_part) {
	int pos=0; // position of the beginning of current block in the original data
	int b_pos=0; // position within the block
	int new_pos=0; // postion after the blocking, always starts at 0 position
	double dx, dy, dz;
	int pid=0;
	// run through all particles
	for(pid=0;pid<n_part;pid++) {
	  pos=0; new_pos=0;
	// we lose at most nc-1 data poins which should not hurt but may be corrected for
	  while(pos+nc-1<n_vals) { 
	      dx=0.0; dy=0.0; dz=0.0;
	      for(b_pos=0; b_pos<nc; b_pos++) {
	        dx+=data[pos+b_pos][3*pid];
	        dy+=data[pos+b_pos][3*pid+1];
	        dz+=data[pos+b_pos][3*pid+2];
	      }
	     data[new_pos][3*pid]=dx/nc;
	     data[new_pos][3*pid+1]=dy/nc;
	     data[new_pos][3*pid+2]=dz/nc;
	     pos+=b_pos;
	     new_pos++;
	  }
	}
	return new_pos; // number of new data values
}

// master function for block correlation algorithm
int block_correlate (double **data, int n_vals, int n_part, int gs, int nc, int l_max, double *result, int n_res) {
	//int pos=0; // our actual position in the data arrays
	int l=0;// current blocking level of the actual data
	//int new_n_vals=n_vals; // number of values in the current (shrunk) data array
	int t; // configuration number aka time 
	int tau_max_plus=nc*gs; // max lag time
	int tau_min=0; // min lag time
	#if TEST
	level_shift=0; // write down current blocking level to a global variable 
	#endif

	// the actual correlation
	for(l=0; l<l_max; l++){
		// perform the simple correlation on this level
		for(t=0; t<n_vals; t++) simple_correlate(result, n_res, data, n_part, n_vals, t, tau_min, tau_max_plus);
		// and the blocking 
		n_vals=block_data(data, nc, n_vals,n_part); // returns number of data values after blocking
		// shift the position in results
		#if TEST
		level_shift+=gs; // write down current blocking level to a global variable 
		#endif
		n_res-=gs;
		result+=gs;
		if(l==0) { 
			tau_min=gs;
			tau_max_plus=nc*gs;
		}
	}
	return 0; // return the pointer to the beginning of the results array
}

#if TEST

int main (int argc, char *argv[]) {
	
	int i; // index variable
	int l_max=0; // maximum blocking level
	int n_res=0; // number of points in the results array
	parse_options (argc, argv); // parse input options

	// open the input and output files
	FILE *outf=fopen(outf_name, "w");
	if (outf==NULL) error("cannot open output file\n");

	// calculate the number of input values if it was not changed by parse_options()
	if(!n_vals) n_vals=get_n_vals(inif_name);
	fprintf(outf, "# number of input values: %ld\n", n_vals);

	double **data; // data array
	double *result; // result array
/*	#if !TEST
	long int n_vals; // number of configs stored
	int n_part; // number of particles per config
	#endif */
	// read the data
	data=read_data(inif_name, n_vals, n_part);  
	
	// print some of the paremeters into the output file
	fprintf(outf, "#n_vals=%ld\n# max_t=%e, dt=%e\n", n_vals, n_vals*dt, dt);
	fflush(outf);
	
	l_max=get_l_max(n_vals, nc, gs);
	if(l_max==1) { n_res=n_vals; } 
	else { n_res=gs*(l_max+1); }
	result=init_res(n_res);
	n_vals_hist=(int *)malloc(n_res*sizeof(int));
	for(i=0;i<n_res;i++) n_vals_hist[i]=0; // set the initial values
	// the actual calculation
	block_correlate(data, n_vals, n_part, gs, nc, l_max, result, n_res);
	fprintf (outf, "#maximum blocking level=%d, %d values in results\n", l_max, n_res);
	
	dump_res(outf, result, n_res, nc, gs, dt, n_vals, n_part); // dump the results
	fclose(outf); 
	return 0;
}

#endif
