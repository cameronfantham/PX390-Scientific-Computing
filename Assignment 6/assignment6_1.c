#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>
#include <string.h>

void read_input(long *Nt, long *Nz, int *E, double *tf, long *Imin);
void read_function(double *Q11, double *Q22, double *H, double *S, double *R, long npoints);

struct band_mat{
  long ncol;        /* Number of columns in band matrix            */
  long nbrows;      /* Number of rows (bands in original matrix)   */
  long nbands_up;   /* Number of bands above diagonal           */
  long nbands_low;  /* Number of bands below diagonal           */
  double *array;    /* Storage for the matrix in banded format  */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix   */
  double *array_inv;/* Store the inverse if this is generated */
  int *ipiv;        /* Additional inverse information         */
};
typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows     = nbands_lower + nbands_upper + 1;
  bmat->ncol       = n_columns;
  bmat->nbands_up  = nbands_upper;
  bmat->nbands_low = nbands_lower;
  bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
  bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
  if (bmat->array==NULL||bmat->array_inv==NULL) {
    return 0;
  }
  /* Initialise array to zero */
  long i;
  for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    bmat->array[i] = 0.0;
  }
  return 1;
};


/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column);
}

void setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  /* Copy bmat array into the temporary store */
  int i,bandno;
  for(i=0;i<bmat->ncol;i++) {
    for (bandno=0;bandno<bmat->nbrows;bandno++) {
      bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
    }
    x[i] = b[i];
  }

  long nrhs = 1;
  long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  return info;
}

int printmat(band_mat *bmat) {
  long i,j;
  for(i=0; i<bmat->ncol;i++) {
    for(j=0; j<bmat->nbrows; j++) {
      printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
    }
  }
  return 0;
}

/*Check that a grid point has valid coordinates */
int is_valid(long j, long p, long J, long P) {
  return (j>=0)&&(j<J)&&(p>=0)&&(p<P);
}


/* Return the 1D element index corresponding to a particular grid point.
   We can rewrite this function without changing the rest of the code if
   we want to change the grid numbering scheme!
   Output: long integer with the index of the point
   Input:
   long j:  The X grid point index
   long k:  The Y grid point index
   long P:  The number of Y points.     //I HAVE REINDEXED THESE
*/
long indx( long j, long p, long P) { // column, row, num rows in matrix
  return j*P + p;
}

long indx_loop( long j, long p, long J, long P) { // column, row, num columns in matrix, num rows in matrix
	if (j < 0) {j += J;};
	if (j >= J) {j -= J;};
	if (p < 0) {p += P;};
	if (p >= P) {p -= P;};
	//printf("%ld %ld \n", j, p);
	return j*P + p;
}


/* Return the 2D point corresponding to a particular 1D grid index */
void gridp(long indx, long P, long *j, long *p) { //index, number of rows, column index, row index,
  *j = indx%P;
  *p = indx - (*j)*P;
}

/* An example of how to use the band matrix routines to solve a PDE:
   The equation solved is related to the steady state solution of the heat
   diffusion equation.
*/
int main() {
  band_mat bmat;
  /* We have a three-point stencil (domain of numerical dependence) of
     our finite-difference equations:
     1 point to the left  -> nbands_low = 1
     1       to the right -> nbands_up  = 1
  */

    //These values will all be read in from input files
	long Nt;
	long Nz;
	int E;
	double tf;
	long Imin;
	read_input(&Nt, &Nz, &E, &tf, &Imin); //Import values



	double npoints = Nt * Nz;


	//Importing in the diffusivity and source term functions
	double *Q11 = malloc(sizeof(double)*npoints);
	double *Q22 = malloc(sizeof(double)*npoints);
	double *H = malloc(sizeof(double)*npoints);
	double *S = malloc(sizeof(double)*npoints);
	double *R = malloc(sizeof(double)*npoints);
	read_function(Q11, Q22, H, S, R, npoints);


  	double dtheta = 0.2; // 0.5 is arbitrary placeholder
  	double dzeta = 0.2; // 0.5 is arbitrary placeholder

 	long j;
	long k;

	double *T, *T_next, *T_temp;  // arrays for y at current and next timestep
	T = malloc(sizeof(double)*npoints);
	T_next = malloc(sizeof(double)*npoints);
	memset(T, 0.0, sizeof(double)*npoints);


	if (E == 0) {
		//SOLVE FOR STEADY STATE
		long nbands_low = 15;
		long nbands_up = 15;
		init_band_mat(&bmat, nbands_low, nbands_up, npoints);
		double *b = malloc(sizeof(double)*(npoints));


		for(k=0; k<Nz; k++) { 			//column by column
			for (j=0; j<Nt; j++) {		//row by row

				double CT_jp_k; //coefficient of T at (j+1, k)
				double CT_jm_k; //coefficient of T at (j-1, k)
				double CT_j_kp; //coefficient of T at (j, k+1)
				double CT_j_km; //coefficient of T at (j, k-1)
				double CT_j_k; //coefficient of T at (j, k) (for RHS, not including LHS terms produced when Dt/dt != 0)


				double dQ11_dtheta = (Q11[indx_loop(j+1, k, Nz, Nt)] - Q11[indx_loop(j-1, k, Nz, Nt)]) / (2*dtheta);
				double dQ22_dzeta  = (Q22[indx_loop(j, k+1, Nz, Nt)] - Q22[indx_loop(j, k-1, Nz, Nt)]) / (2*dzeta);
				double Hval 	= H[indx_loop(j, k, Nz, Nt)];
				double Q11val 	= Q11[indx_loop(j, k, Nz, Nt)];
				double Q22val 	= Q22[indx_loop(j, k, Nz, Nt)];
				double Rval 	= R[indx_loop(j, k, Nz, Nt)];
				double Sval 	= S[indx_loop(j, k, Nz, Nt)];

				CT_jp_k = (Hval*dQ11_dtheta*(0.5/dtheta)) + Hval*Q11val*(1/(dtheta*dtheta));
				CT_jm_k = (Hval*dQ11_dtheta*(-0.5/dtheta)) + Hval*Q11val*(1/(dtheta*dtheta));
				CT_j_kp = (Hval*dQ22_dzeta*(0.5/dzeta)) + Hval*Q22val*(1/(dzeta*dzeta));
				CT_j_km = (Hval*dQ22_dzeta*(-0.5/dzeta)) + Hval*Q22val*(1/(dzeta*dzeta));
				CT_j_k = (Hval*Q11val*(-2.0/(dtheta*dtheta)) + Hval*Q22val*(-2.0/dzeta*dzeta) - Rval);



				setv(&bmat,indx_loop(j, k, Nz, Nt),	indx_loop(j, k, Nz, Nt),	CT_j_k);
				setv(&bmat,indx_loop(j, k, Nz, Nt),	indx_loop(j+1, k, Nz, Nt),	CT_jp_k);
				setv(&bmat,indx_loop(j, k, Nz, Nt),	indx_loop(j-1, k, Nz, Nt),	CT_jm_k);
				setv(&bmat,indx_loop(j, k, Nz, Nt),	indx_loop(j, k+1, Nz, Nt),	CT_j_kp);
				setv(&bmat,indx_loop(j, k, Nz, Nt),	indx_loop(j, k-1, Nz, Nt),	CT_j_km);


				b[indx_loop(j, k, Nz, Nt)] = -1 * Sval;


			}
		}


		solve_Ax_eq_b(&bmat, T, b);





	}

	else {
		//SOLVE FOR TIME-EVOLUTION CASE.

		double dt = 0.1; // arbitrary at the moment
		double stableTimestep = 0.1;
		long divisor = Imin; // Keep increasing this until endtime/divisor < minimum timestep for stability.
		long i;


		while (dt > stableTimestep) {
			divisor++;
			dt = tf/divisor;
		}


		for (i = 0; i < divisor; i++) {


			for(k=0; k<Nz; k++) { 			//column by column
				for (j=0; j<Nt; j++) {		//row by row

				double CT_jp_k; //coefficient of T at (j+1, k)
				double CT_jm_k; //coefficient of T at (j-1, k)
				double CT_j_kp; //coefficient of T at (j, k+1)
				double CT_j_km; //coefficient of T at (j, k-1)
				double CT_j_k; //coefficient of T at (j, k) (for RHS, not including LHS terms produced when Dt/dt != 0)

				double dQ11_dtheta = (Q11[indx_loop(j+1, k, Nz, Nt)] - Q11[indx_loop(j-1, k, Nz, Nt)]) / (2*dtheta);
				double dQ22_dzeta  = (Q22[indx_loop(j, k+1, Nz, Nt)] - Q22[indx_loop(j, k-1, Nz, Nt)]) / (2*dzeta);
				double Hval 	= H[indx_loop(j, k, Nz, Nt)];
				double Q11val 	= Q11[indx_loop(j, k, Nz, Nt)];
				double Q22val 	= Q22[indx_loop(j, k, Nz, Nt)];
				double Rval 	= R[indx_loop(j, k, Nz, Nt)];
				double Sval 	= S[indx_loop(j, k, Nz, Nt)];

				CT_jp_k = (Hval*dQ11_dtheta*(0.5/dtheta)) + Hval*Q11val*(1/(dtheta*dtheta));
				CT_jm_k = (Hval*dQ11_dtheta*(-0.5/dtheta)) + Hval*Q11val*(1/(dtheta*dtheta));
				CT_j_kp = (Hval*dQ22_dzeta*(0.5/dzeta)) + Hval*Q22val*(1/(dzeta*dzeta));
				CT_j_km = (Hval*dQ22_dzeta*(-0.5/dzeta)) + Hval*Q22val*(1/(dzeta*dzeta));
				CT_j_k = (Hval*Q11val*(-2.0/(dtheta*dtheta)) + Hval*Q22val*(-2.0/dzeta*dzeta) - Rval);


				double RHS = T[indx_loop(j, k, Nz, Nt)] + dt*( CT_jp_k*T[indx_loop(j+1, k, Nz, Nt)] +
        CT_jm_k*T[indx_loop(j-1, k, Nz, Nt)] + CT_j_kp*T[indx_loop(j, k+1, Nz, Nt)] +
        CT_j_km*T[indx_loop(j, k-1, Nz, Nt)] + Sval);

				T_next[indx_loop(j, k, Nz, Nt)] =  RHS / ( 1 - dt*CT_j_k );

				}
			}

			T_temp = T;
			T = T_next;
			T_next = T_temp;

		}
	}


	FILE *output_file;
	output_file	= fopen("output.txt", "w");
	// If issue with file, return 1
	if(output_file == NULL) {
		return 1;
	}
	//print numbers to file
	for(j=0; j<Nz; j++) {
		for(k=0; k<Nt; k++) {
			fprintf(output_file, "%lf \n", T[indx_loop(j, k, Nz, Nt)]);
		}
	}
	//close file and return 0
	fclose(output_file);

	return 0;
}





/*  Read in:
Num. of theta points,
Num. of zeta points,
Boolean steady state/time evolution (E)
final time (for time evolution)
Minimum iteration count
*/
void read_input(long *Nt, long *Nz, int *E, double *tf, long *Imin) {
   FILE *infile;
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(5!=fscanf(infile,"%ld %ld %d %lf %ld\n",Nt,Nz,E,tf,Imin)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}

//  Read in functions: Q11, Q22, H, S, R;
void read_function(double *Q11, double *Q22, double *H, double *S, double *R, long npoints) {
   FILE *infile;
   if(!(infile=fopen("coefficients.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   int i = 0;
   while (i < npoints) {
		fscanf(infile,"%lf %lf %lf %lf %lf \n", &Q11[i], &Q22[i], &H[i], &S[i], &R[i]);
		i++;
   }

   fclose(infile);
}
