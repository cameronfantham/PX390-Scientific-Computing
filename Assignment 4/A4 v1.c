#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <lapacke.h>

void read_input(double *L, int *N, double *t_F, double *t_D, double *T_0, double *A, double *f, double *P);
void write_output(double ctime, double x, double T);

double G_func(double dTdx, double P) {
  if(dTdx>0) {
    return -powl(fabs(dTdx),P);
  } 
  if (dTdx==0) {
    return 0;
  }
  if (dTdx<0) {
    return powl(fabs(dTdx),P);
  }
  return 0;
}

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
  double f;
  //Nonlinearity exponent.
  double P;
  // Read in from file.
  read_input(&L, &N, &t_F, &t_D, &T_0, &A, &f, &P);
  // Grid spacing.
  double dx = L/(N-1);
  // Time step.
  double dt = 0.5*dx*dx/(P*powl(abs(B),P-1));
  // ************
  // Grid Storage
  // ************
  // Arrays for y at current and next timestep.
  double* T, *T_next;
  // Allocate memory according to size of N.
  T       = malloc(sizeof(double)*N);
  T_next  = malloc(sizeof(double)*N);

  int j;
  double x;

  // **************
  // Initialisation 
  // **************
  double ctime = 0.0;
  double next_output_time = t_D;
  for(j=0;j<N;j++) {
    x=j*dx;
    T[j] = T_0*(1 - x/L);
    T_next[j] = 0.0;
  }
  // *******************
  // Loop over timesteps 
  // *******************
  while (ctime<t_F){
    double dt_0 = dt;
    // If we would go past the next output step, reduce the timestep so we exactly hit the output time.
    if (ctime+dt_0>next_output_time) {
      dt_0 = next_output_time - ctime;
    }
    T[0] = A*(1 - cos(f*ctime)) + T_0;
    //loop over points 
    for (j=1; j<N; j++) {
      x = j*dx;
      int jm = j-1;
      // Need upwinding for stability.
      double dTdx;
      if (x==L) {
        dTdx = -T_0/L; 
      } else { 
        dTdx = (T[j] - T[jm])/dx;
      }
      T_next[j] = T[j] - dt_0*(G_func(dTdx,P))/dx;
    }
    
    // Swap values of y and y_next pointers to efficiently
    // 'copy' arrays without element-by-element assignment.
    double* T_temp;
    T_temp = T;
    T = T_next;
    T_next = T_temp;

    // Increment time.   
    ctime = ctime + dt_0;
    // Output values.
    if (ctime == next_output_time) {
      for (j=0; j<N; j++) {
	      x = j*dx;
	      write_output(ctime,x,T[j]);
      }
      next_output_time += t_D;
    }
  }
  free(T);
  free(T_next);
}

void read_input(double *L, int *N, double *t_F, double *t_D, double *T_0, double *A, double *f, double *P) {
   FILE *infile;
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(8!=fscanf(infile,"%lf %d %lf %lf %lf %lf %lf %lf",L,N,t_F,t_D,T_0,A,f,P)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}

void write_output(double ctime, double x, double T) {
   FILE *infile;
   if(!(infile=fopen("output.txt","w"))) {
       printf("Error opening file\n");
       exit(1);
   }
   fprintf(infile,"%lf %lf %lf",ctime,x,T);
   if(3!=fprintf(infile,"%lf %lf %lf",ctime,x,T)) {
       printf("Error writing parameters to file\n");
       exit(1);
   }
   fclose(infile);
}