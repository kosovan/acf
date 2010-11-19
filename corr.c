/* C O R R E L A T O R Version 0.9, 11.3.2006 */

/* This is really an ancient code and has nothing to do with autocorr anymore, except for sharing logical structure of the groupping */
#include<corr.h>

#define TEST 1

#ifndef TYPE_DATA
#define TYPE_DATA
typedef struct Data{
	double t; double I; // time and intensity
}data;
#endif

#if TEST
// global variables
int nc=2;			// groupping number
double dt=1.0;			// dt
int n_vals=0;			// number of values
int gs=8;			// group size
char *inif_name="corr.in";	// input file name
char *outf_name="corr.out";	// output file name
int ascii_flag=0;		// if set to zero, expect t3r2ascii input
int linear_flag=0;		// if set to zero, do blocking correlation, otherwise do linear correlation
long int max_n_vals=0; 
#endif

void error (char *s) /* error messages */
{
        fprintf (stderr, "\ncorrelator: %s\n\t the program was terminated\n\n", s);
        fflush (stdout);
        exit (1);
}

// the maximum number of values, if also all the empty channels were counted
long int get_max_n_vals(data *d, double dt, int n_vals){
	return (long int)floor((d[n_vals-1].t-d[0].t)/dt);
}

// calculate the maximum possible blocking level
int get_l_max (long int n_vals, int nc, int gs, long int max_n_vals, double dt, data *d){
	int i=0; 
	max_n_vals/=gs; // divide it by the basic groupping gs
	while(max_n_vals/nc){ 
		max_n_vals/=nc;
		i++;
	}
	return i;
}

// read the output from the t3r2ascii without the header
data *read_simple_ascii(char *inif_name, int n_vals) {
	FILE *f=fopen(inif_name, "r");
	if (f==NULL) error("cannot open inif\n");
	int i=0;
	data *d;
	d=(data*)malloc(n_vals*sizeof(data));
	while(i<n_vals){
		fscanf(f,"%lf", &d[i].t);
		d[i].I=1.0;
		i++;
	}
	fclose(f);
	return d;
}

// read the output from the t3r2ascii without the header
data *read_t3r_ascii(char *inif_name, int n_vals){
	FILE *f=fopen(inif_name, "r");
	if (f==NULL) error("cannot open inif\n");
	int i=0;
	data *d;
	d=(data*)malloc(n_vals*sizeof(data));
	while(i<n_vals){
		fscanf(f,"%*f %*f %lf %*f %*f", &d[i].t);
		d[i].I=1.0;
		i++;
	}
	fclose(f);
	return d;
}

// perform linear correlation of data1[0] with items in data2 in the interval (min_dt, max_dt)
int simple_correlate(data *result, data data1, data *data2, int pos2, int n_vals2, double dt, double max_dt, double min_dt){
	/* *result points to such a position in the results, that corresponds to act_dt=dt
	   dt is the current dt corresponding to the blocking level; 
	   simple_correlate does not bother about the blocking */

	double act_dt=0.0; // actual dt (between the two processed values)
	double t1=data1.t; // to make notation shorter
	int result_pos=0; // position in the results field
	int minpos2=pos2; // maximum position where still data2[minpos2].t<t1+dt
	double ival;
	// avoid the rounding error
	double dif=0.0001*dt; // stored dt may differ from the dt in the result where it is stored
	max_dt+=0.001*dt;
	min_dt-=0.001*dt;
	while(data2[pos2].t<t1 && pos2<n_vals2) pos2++; // move in the 2nd array in order to be beyond t1
	minpos2=pos2; // note where t2>=t1 for the first time
	for ( ;pos2<n_vals2; pos2++){
		act_dt=data2[pos2].t-t1; // calculate the time difference between the stored data
		if (act_dt<min_dt && minpos2<n_vals2){ // if it is less than min, go for the next
		       minpos2++; 
	       	       continue; 
		}
		if (act_dt>max_dt) break; // if it is more than max, we are finished
		ival=data1.I*(data2[pos2].I); // intensity value
		result_pos=(int)((act_dt-min_dt)/dt); // find proper position for the result
		// test if the result is stored in a proper place - for debugging only
		if(fabs(result[result_pos].t-act_dt)>dif) { 
			printf ("misplaced result: act_dt=%.20e, to be stored in t=%.20e\n\
min_dt=%e, act_dt=%e, dt=%e result_pos=%d\n", act_dt, result[result_pos].t, min_dt, act_dt, dt, result_pos);
			error("");
		} 
		result[result_pos].I+=ival; // add ival to the results
	}
	return minpos2; // remember where data2[minpos2].t<t1, to be used for next call of simple_correlate
}

// print the results to a file
int dump_res(FILE *f, data *res, int n_vals){
	int i;
	fprintf(f,"#time\t intensity\n");
	// the first value is crippled, print it out with a hash mark
	fprintf(f, "# %e\t %lf\n", res[0].t, res[0].I);
	for(i=1; i<n_vals; i++){
		fprintf(f, "%e\t %lf\n", res[i].t, res[i].I);
	}
	fflush(f);
	return 0;
}

// initialise the array where results will be stored
data* init_res(int n_vals, int nc, int gs, double dt){
	int i;
	data *result;
	result=(data*)malloc(n_vals*sizeof(data));
	double t=0.0;
	int shift=gs, act_shift=nc*shift+1;
	result[0].I=result[0].t=0.0;
	for(i=1; i<n_vals; i++){
		if(i==act_shift){ // if we have crossed the block border, change dt
			dt*=(double)nc;
			act_shift+=shift;
		}
		t+=dt;
		result[i].I=0.0;
		result[i].t=t;
	}
	//dump_res(result, n_vals);
	return result;
}

// initialise the array where results will be stored for the linear correlation
data* init_res_lin(double n_vals, double dt){
	int i;
	data *result;
	result=(data*)malloc(n_vals*sizeof(data));
	double t=0.0;
	for(i=0; i<n_vals; i++){
		result[i].I=0.0;
		result[i].t=t;
		t+=dt;
	}
	//dump_res(result, n_vals);
	return result;
}

// calculate the average intensity of the whole dataset
double avg_intensity(data *d, int n_vals, double dt, long int max_n_vals){
	int i=0;
	double result=0.0;
	while(i<n_vals){ // add all the intensities
		result+=d[i].I;
		i++;
	}
	result/=(double)max_n_vals; // divide them by the number of channels spanning the entire time interval
	return result;
}

// normalize results by dividing the by the square of the average intensity and by the number of values which were adde in each field
int normalize_res(data *res, int n_res, int n_vals, double avg, double dt, long int max_n_vals){
	int i;
	double n_act; // the actual number of the values which contribute to the results values
	avg*=avg; // we need the square of average intensity
	for (i=0; i<n_res; i++) {
		n_act=(double)max_n_vals-res[i].t/dt-1.0;
		if (res[i].I<0){
			fprintf (stderr, "res in %d %f < 0\n", i, res[i].I);
			error("error in results");
		}	       
		res[i].I/=avg*n_act;
	}
	return 0;
}

// perform the blocking of data
int block_data(data *d, double act_dt, int nc, int n_vals){
	int pos=0; // position in the original data
	int b_pos=0; // postion after the blocking, always starts at 0 position
	double min_t=0.0; // starting time for the block
	double max_dt=act_dt*(double)nc;
	double max_t=0.0; // and where it ends
	double temp=0.0;
	double sqrt_nc=sqrt(nc);
	while(pos<n_vals){ // run through the whole array
		min_t=max_t;
		max_t=min_t+max_dt;
		if(!(d[pos].t<max_t)) continue; // continue if there is no value for the current block
		// cumulate data in the block
		temp=0.0;
		while(d[pos].t<max_t && pos<n_vals){
			if(pos<b_pos) {
				fprintf(stderr, "pos=%d < b_pos=%d", pos, b_pos);
				error("wrong blocking");
			}
			temp+=d[pos++].I;
		}
		temp/=sqrt_nc;
		/* FIXME
		   					W A R N I N G 
		   I am not sure why the total intensity should be divided by the root of the number of channels
		   which are being goupped. I would expect to simply average it, i.e. divide by nc. Perhaps it is 
		   becasue in correlations it is multiplied and hence it doubles once again. 

		   But it works!!!

		 										PK 27.2.2006
		 */
		d[b_pos].I=temp;
		d[b_pos++].t=min_t;
	}
	return b_pos;
}

/* data_i = i-th data structure,
   n_vals_i = number of values in the i-th structure
   dt = the fundamental time interval (minimum distance between successive data points)
   nc = number of members in a blocking operation       
   gs = group size, the minimum correlation interval for the l-th level is gs*actutal_dt(l)
   l_max = maximum level of blocking
   cross_flag = flag if cross correlation is being done (otherwise data 1 & 2 point to the same array and double blocking must be avoided
*/

data* block_correlate (data *data1, int n_vals1, data *data2, int n_vals2, double dt, int gs, int nc, int l_max, int n_res, int cross_flag) {
	int pos1=0, pos2=0; // our actual position in the data arrays
	int l=0;// current blocking level of the actual data
	double act_dt=dt; // actual dt for the current use
	data *result; // array to store the results
	double max_dt=act_dt*((double)gs*nc-1);
	double min_dt=0.0;

	//prepare the results array
	result=init_res(n_res, nc, gs, dt);
	//dump_res(result, n_res);

	// the actual correlation
	for(l=0; l<l_max; l++){
		// perform the simple correlation on this level
		for(pos1=pos2=0; pos1<n_vals1; pos1++) {
			pos2=simple_correlate(result, data1[pos1], data2, pos2, n_vals2, act_dt, max_dt, min_dt);
		}
		// and the blocking 
		n_vals1=block_data(data1, act_dt, nc, n_vals1); // returns number of data values after blocking
	       	if(cross_flag) n_vals2=block_data(data2, act_dt, nc, n_vals2); // if(!cross_flag), double blocking would occur because data1 and daat2 point to the same array
		else n_vals2=n_vals1;
		//adjust the time step
		min_dt=max_dt+act_dt;
		act_dt*=nc;
		max_dt=min_dt+act_dt*((double)gs-1);
		
		// shift the position in results
		if(l==0) result+=nc*gs;
		else result+=gs;
	}
	return result-(l_max+1)*gs; // return the pointer to the beginning of the results array
}

// perform the linear corelation of the results, even though it takes so long - for consistency check of the blocking correlator
data *lin_correlate(data *data1, int n_vals1, data *data2, int n_vals2, double dt) {
	int pos1=0, pos2=0; // our actual position in the data arrays
	double max_dt=0.5*(data2[n_vals2-1].t-data1[0].t); // absolute maximum for dt
	double min_dt=0.0;
	int n_res=(int)ceil(max_dt/dt); // number of positions in the results
	data *result; // array to store the results
	//initialise the results
	result=init_res_lin(n_res, dt);
	// perform the simple correlation on this level
	for(pos1=pos2=0; pos1<n_vals1; pos1++) {
		if(!(pos1%100)){ // so that I know tat it is moving on
			putchar('#'); 
			fflush(stdout);
		}
		//printf("pos1=%d, div by 100=%d\n", pos1, (pos1%100));
		pos2=simple_correlate(result, data1[pos1], data2, pos2, n_vals2, dt, max_dt, min_dt);
	}
	return result;
}

void print_help()
{
	printf ("\nusage: corr [options]\n \
	options:\n \
	-f <filename>    input file (default: corr.in)\n \
	-o <filename>    output file (default: corr.out)\n \
	-n <integer>     number of input data values, use if you want to override the default automatic determination\n \
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
			case 'b' :
				nc=atoi(*++argv);
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
			case 'A':
				ascii_flag=1; 
				++argc;
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

#if TEST
int main (int argc, char *argv[]) {
	
	int l_max=0; // maximum blocking level
	int n_res=0; // number of points in the results array
	parse_options (argc, argv); // parse input options

	// open the input and output files
	FILE *outf=fopen(outf_name, "w");
	if (outf==NULL) error("cannot open output file\n");

	// calculate the number of input values if it was not changed by parse_options()
	if(!n_vals) n_vals=get_n_vals(inif_name);
	fprintf(outf, "# number of input values: %d\n", n_vals);

	data *data1, *data2, *result;	// arrays where input data and results are stored
	double avg=0.0; // average intensity
	int cross_flag=0; // if cross correlation is being performed
	int n_vals1=n_vals, n_vals2=n_vals; // FIXME - for cross-correlation, the 2 values will be independent
	// read the data
	if (ascii_flag) data1=read_simple_ascii(inif_name, n_vals);  
	else data1=read_t3r_ascii(inif_name, n_vals);  
	if(!cross_flag) data2=data1;
	else error ("cross flag set; the program is not ready for cross correlation yet");
	
	// calculate the maximum possible number of values (as if all the counts were non-zero)
	// FIXME for testing purposes, we define max_n_vals as a global variable
	//long int max_n_vals=get_max_n_vals(data1, dt, n_vals); 
	max_n_vals=get_max_n_vals(data1, dt, n_vals); 
	
	// print some of the paremeters into the output file
	fprintf(outf, "#max_n_vals=%ld\n# max_t=%e, min_t=%e\n", max_n_vals, data1[n_vals-1].t, data1[0].t);
	fflush(outf);
	
	if(linear_flag) { 
	  l_max=1;
	  n_res=max_n_vals;
	  result=lin_correlate(data1, n_vals1, data2, n_vals2, dt);
	} else { 
	  l_max=get_l_max(n_vals1, nc, gs, max_n_vals, dt, data1);
	  n_res=gs*(l_max+1);
	  // the actual calculation
	  result=block_correlate(data1, n_vals1, data2, n_vals2, dt, gs, nc, l_max, n_res, cross_flag);
	  normalize_res(result, n_res, n_vals, avg, dt, max_n_vals); // normalization
	}
	fprintf (outf, "#maximum blocking level=%d, %d values in results\n", l_max, n_res);
	avg=avg_intensity(data1, n_vals, dt, max_n_vals); // average intensity
	
	dump_res(outf, result, n_res); // dump the results
	fclose(outf);
	return 0;
}
#endif

// the working version of my correlator program
/* // once again for the linear correlation
   inif=fopen("./dummy2", "r");
   if (inif==NULL) error("cannot open inif\n");
   data1=read_acf_data(inif, n_vals, n_skips, dt);
   data2=data1;
   max_dt=0.5*(data2[n_vals2-1].t-data1[0].t); // absolute maximum for dt
   n_res=(int)floor(max_dt/dt/2.0); // number of points in the results array
   result=lin_correlate(data1, n_vals1, data2, n_vals2, dt);
   normalize_res(result, n_res, n_vals, avg, dt);
   dump_res(result, n_res);
*/
