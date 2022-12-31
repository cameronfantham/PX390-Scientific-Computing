#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void read_input(double *L, int *N, double *t_F, double *t_D, double *T_0, double *A, double *omega, double *P);
void write_output(double ctime, double x, double T);

int main(void) {
  // **********
  // Parameters
  // **********
  // Length of domain.
  double L;
    // Number of grid points.                  
  int N;
  //  Length of time to run simulation.
  double t_F;
  // Timestep for diagnostic output.
  double t_D;
  // Initial condition parameter.
  double T_0;
  //Boundary: Amplitude of oscillation.
  double A;
  //Boundary: Frequency of oscillation.
  double omega;
  //Nonlinearity exponent.
  double P;
  // Read in from file.
  read_input(&L, &N, &t_F, &t_D, &T_0, &A, &omega, &P);
  // Grid spacing.
  double dx = L/(N-1);
  // Maximum Local Derivative. (stable for  B = n*(T_0 + 1) where n >= 2)
  double B = 8*(T_0 + 1);
  // Time step.
  double dt = 0.5*dx*dx/(P*powl(fabs(B),P-1));
  // ************
  // Grid Storage
  // ************
  // Arrays for y at current and next timestep.
  double* T, *T_next;
  // Allocate memory according to size of N.
  T       = malloc(sizeof(double)*N);
  T_next  = malloc(sizeof(double)*N);

  // **************
  // Initialisation 
  // **************
  int k;
  double ctime = 0.0;
  double next_output_time = t_D;
  for(k=0;k<N;k++) {
    //Initial condition.
    T[k] = T_0*(1 - k*dx/L);
    T_next[k] = 0.0;
    write_output(ctime,k*dx,T[k]);
  }
  // *******************
  // Loop over timesteps 
  // *******************
  while (ctime<t_F){
    //Boundary condition in x.
    T[0] = A*(1 - cos(omega*ctime)) + T_0;
    T_next[0] = T[0];
    // If we go past next output step, reduce timestep so we exactly hit output time.
    double dt_0 = dt;
    if (ctime+dt_0>next_output_time) {
      dt_0 = next_output_time - ctime;
    }
    //Loop over interior points. 
    double dTdx;
    for (k=1; k<N-1; k++) {
      // Calculate T_next.
      dTdx = (T[k+1] - T[k-1])/(2*dx);
      T_next[k] = T[k] + P*powl(fabs(dTdx),P-1)*dt_0*(T[k+1]+T[k-1]-2*T[k])/(dx*dx);
    }
    // Boundary condition dT/dx = -T_0/L at x=L.
    dTdx = -T_0/L;
    T_next[N-1] = T_next[N-2] + dx*dTdx;
    
    // Swap values of y and y_next pointers to efficiently
    // 'copy' arrays without element-by-element assignment.
    double* T_temp;
    T_temp = T;
    T = T_next;
    T_next = T_temp;

    // Increment time.   
    ctime += dt_0;
    // Output values.
    if (ctime == next_output_time) {
      for (k=0; k<N; k++) {
	      write_output(ctime,k*dx,T[k]);
      }
      next_output_time += t_D;
    }
  }
  free(T);
  free(T_next);
}

void read_input(double *L, int *N, double *t_F, double *t_D, double *T_0, double *A, double *omega, double *P) {
   FILE *infile;
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file (reading)\n");
       exit(1);
   }
   if(8!=fscanf(infile,"%lf %d %lf %lf %lf %lf %lf %lf",L,N,t_F,t_D,T_0,A,omega,P)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}

void write_output(double ctime, double x, double T) {
   FILE *outfile;
   if(!(outfile=fopen("output.txt","a"))) {
       printf("Error opening file (writing)\n");
       exit(1);
   }
   fprintf(outfile,"%g %g %g \n",ctime,x,T);
   fclose(outfile);
}
