// Bug Fix List
/*
1) #include <stdio.h>
2) x>=L/4 && x<L/2
3) added for loop to set initial y to 0:
   for(j=0;j<nx;j++) {
     y[j] = 0.0;
     y_next[j] = 0.0;
   }
4) y[jp] - y[j] changed to y[j] - y[jm]
5) removed unused variables jp and output
6) added if loop for periodicity:
   if (jm<0) {
     dydx = (y[j] - y[nx-1])/dx; 
     } else { 
       dydx = (y[j] - y[jm])/dx;
   }
5) y_next[j] = y[j] - dt*C*dydx + S_func(x,L)
   
   changed to

   y_next[j] = y[j] - dt0*(C*dydx - S_func(x,L))
6) added: "double* y_temp;
           y_temp = y;"
7) remove double from "ctime = ctime + dt0;"
8) "malloc(nx)" changed to "malloc(sizeof(double)*nx)" 
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void read_input(double *C, double *L, int *nx, double *t_F,double *t_out);

double S_func(double x, double L) {
  if(x>=L/4 && x<L/2) {
    return 1.0;
  } else {
    return 0.0;
  }
}

int main(void) {
  // **********
  // Parameters
  // **********
  // Number of grid points.                  
  int nx;
  // Length of domain.
  double L;
  // Equation coefficient.
  double C;
  //  Length of time to run simulation.
  double t_F;
  // How frequently in time to output.
  double output_timestep;
  // Read in from file.
  read_input(&C, &L, &nx, &t_F, &output_timestep);
  // Grid spacing.
  double dx = L/(nx-1);		      
  // Time step.
  double dt = 0.5*dx/C;
  // ************
  // Grid Storage 
  // ************
  // Arrays for y at current and next timestep.
  double* y, *y_next;  
  // Allocate memory according to size of nx.
  y       = malloc(sizeof(double)*nx);
  y_next  = malloc(sizeof(double)*nx);

  int j;
  double x;

  // **************
  // Initialisation 
  // **************
  double ctime = 0;
  double next_output_time = output_timestep;
  for(j=0;j<nx;j++) {
    y[j] = 0.0;
    y_next[j] = 0.0;
  }
  // *******************
  // Loop over timesteps 
  // *******************
  while (ctime<t_F){
    double dt0 = dt;
    // If we would go past the next output step, reduce the timestep so we exactly hit the output time.
    if (ctime+dt0>next_output_time) {
      dt0 = next_output_time - ctime;
    }
    
    //loop over points 
    for (j=0; j<nx; j++) {
      x = j*dx;
      int jm = j-1;
      // Need upwinding for stability.
      double dydx;
      if (jm<0) {
        dydx = (y[j] - y[nx-1])/dx; 
      } else { 
        dydx = (y[j] - y[jm])/dx;
      }
      y_next[j] = y[j] - dt0*(C*dydx - S_func(x,L));
    }
    
    // Swap values of y and y_next pointers to efficiently
    // 'copy' arrays without element-by-element assignment.
    double* y_temp;
    y_temp = y;
    y = y_next;
    y_next = y_temp;

    // Increment time.   
    ctime = ctime + dt0;
    // Output values.
    if (ctime == next_output_time) {
      for (j=0; j<nx; j++) {
	      x = j*dx;
	      printf("%g%g%g \n",ctime,x,y[j]);
      }
      next_output_time += output_timestep;
    }
  }
  free(y);
  free(y_next);
}

void read_input(double *C, double *L, int *nx, double *t_F, double *t_out) {
   FILE *infile;
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(5!=fscanf(infile,"%lf %lf %d %lf %lf",C,L,nx,t_F,t_out)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}