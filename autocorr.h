#define TEST 1

#if TEST
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void error (char *s); /* error messages */

double **read_data(char *inif_name, int n_vals, int n_part);

void print_help();

int parse_options (int argc, char *argv[]);

int get_n_vals (char *inif_name);

int main (int argc, char *argv[]);

/* print the results to a file
@param n_vals	number of values in the time series -- to be used for normalization
@param n_part	number of particles in configs

*/
int dump_res(FILE *f, double *result, int n_res, int nc, int gs, double dt, int n_vals, int n_part);

// initialise the array where results will be stored
double* init_res(int n_vals);

#endif

/* Calculate the maximum possible blocking level for the entire data set

@param n_vals	number of values in the time series
@param nc	number of entries to be groupped together in one groupping step	
@param gs	group size
*/
int get_l_max (long int n_vals, int nc, int gs);

/* operation to be performed on the data
currently we will be computing only MSD
*/
double my_operation (double **data, int pid, int t, int tau);

/* Compute correlation without blocking of data[t] with items in data[t+tau_min] up to data[t+tau_max]

	@param result	array in which the result is stored
	@param n_res	number of entries in the results array to allow for overflow checking
	@param data	the actual input data (particle coordinates)
	@param n_part	number of particles in the data array
	@param t_max	number of configs in the data array
	@param t	index of the 1st entry where to start correlating
	@param tau_min	index of the 2nd entry where to start correlating
	@param tau_max	index of the 2nd entry BEFORE which to stop correlating
	   
	assume data stored in the form: data[config][3*pid=pid_x] 
	  				data[config][3*pid+1=pid_y] 
	  				data[config][3*pid+2=pid_z] 

*/ 
int simple_correlate(double *result, int n_res, double **data, int n_part, int n_vals, int t, int tau_min, int tau_max);

// perform the blocking of data
int block_data(double **data, int nc, int n_vals, int n_part);

/* master function for block correlation algorithm
   @param data		input data containing only coordinates of those particles stored in configs, which we want to analyze
   @param n_part	number of particles contained in data array
   @param n_vals	number of configs contained in data
   @param nc 		number of data entries joined in one blocking operation       
   @param gs 		group size, the minimum correlation interval for the l-th level is gs*tau_min(l)
   @param l_max 	maximum level of blocking
   @param res	 	array to store the results, double **res should expect data in the form res[0][t], res[1][acf(t)]
   @param n_res 	number of entries in results, space should be already allocated

*/
int block_correlate (double **data, int n_vals, int n_part, int gs, int nc, int l_max, double *result, int n_res);

