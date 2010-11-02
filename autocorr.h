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

// initialise the array where results will be stored
double* init_res(int n_vals);

#endif

/** Print the results to a file, including post-processing corrections
@param f
@param result
@param n_res
@param nc
@param gs
@param dt
@param n_vals	number of values in the time series -- to be used for normalization
@param n_part	number of particles in configs
@param result_jack	result array for jacknifing supposed to have as many entries as n_res
@param n_jack	number of jacknifing blocks
@param jack_max_pid	up to which particle ID Jacknifing should be performed (useful in case that n_part%n_jack!=0)

*/
int dump_res(FILE *f, double *result, int n_res, int nc, int gs, double dt, int n_vals, int n_part, double **result_jack, int n_jack, int jack_max_pid);

/** Calculate the standard deviation from the jacknife method
  @param result_jack	array where output of jacknife method is stored
  @param n_jack		number of blocks for jacknifing
  @param n_res		length of the results array
*/
int stdev_jack (double **result_jack, int n_jack, int n_res);

/** Calculate the maximum possible blocking level for the entire data set

@param n_vals	number of values in the time series
@param nc	number of entries to be groupped together in one groupping step	
@param gs	group size
*/
int get_l_max (long int n_vals, int nc, int gs);

/* operation to be performed on the data
currently we will be computing only MSD
*/
double my_operation (double **data, int pid, int t, int tau);

/** Compute correlation without blocking of data[t] with items in data[t+tau_min] up to data[t+tau_max]

@param result	array in which the result is stored
@param n_res	number of entries in the results array to allow for overflow checking
@param result_jack	result array for jacknifing supposed to have as many entries as n_res
@param n_jack	number of jacknifing blocks
@param data	the actual input data (particle coordinates)
@param n_part	number of particles in the data array
@param t_max	number of configs in the data array
@param t	index of the 1st entry where to start correlating
@param tau_min	index of the 2nd entry where to start correlating
@param tau_max	index of the 2nd entry BEFORE which to stop correlating
@param jack_max_pid	up to which particle ID Jacknifing should be performed (useful in case that n_part%n_jack!=0)
	   
Assume data stored in the form: \\
  data[config][3*pid=pid_x] \\
  data[config][3*pid+1=pid_y] \\
  data[config][3*pid+2=pid_z] \\
*/ 
int simple_correlate(double *result, int n_res, double **data, int n_part, int n_vals, int t, int tau_min, int tau_max, double **result_jack, int n_jack, int jack_max_pid);

// perform the blocking of data
int block_data(double **data, int nc, int n_vals, int n_part);

/** Master function for block correlation algorithm
  @param data		input data containing only coordinates of those particles stored in configs, which we want to analyze
  @param n_part	number of particles contained in data array
  @param n_vals	number of configs contained in data
  @param nc 		number of data entries joined in one blocking operation       
  @param gs 		group size, the minimum correlation interval for the l-th level is gs*tau_min(l)
  @param l_max 	maximum level of blocking
  @param result	 	array to store the results, double **res should expect data in the form res[0][t], res[1][acf(t)]
  @param n_res 	number of entries in results, space should be already allocated
  @param result_jack	result array for jacknifing
  @param n_jack	number of jacknifing blocks
  @param jack_max_pid	up to which particle ID Jacknifing should be performed (useful in case that n_part%n_jack!=0)
*/
int block_correlate (double **data, int n_vals, int n_part, int gs, int nc, int l_max, double *result, int n_res, double **result_jack, int n_jack, int jack_max_pid);
