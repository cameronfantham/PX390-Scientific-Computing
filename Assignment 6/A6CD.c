#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>
#include <math.h>

struct band_mat{
  long ncol;         /* Number of columns in band matrix         */
  long nbrows;       /* Number of rows (bands in original matrix)*/
  long nbands_up;    /* Number of bands above diagonal           */
  long nbands_low;   /* Number of bands below diagonal           */
  double *array;     /* Storage for the matrix in banded format  */
  /* Internal temporary storage for solving inverse problem     */
  long nbrows_inv;   /* Number of rows of inverse matrix         */
  double *array_inv; /* Store the inverse if this is generated   */
  int *ipiv;         /* Additional inverse information           */
};
typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory, and set the parameters */
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
  bmat->array      = (double *) calloc(bmat->nbrows*bmat->ncol,sizeof(double));
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  bmat->array_inv  = (double *) calloc((bmat->nbrows+bmat->nbands_low)*bmat->ncol,sizeof(double));
  bmat->ipiv       = (int *) calloc(bmat->ncol,sizeof(int));
  if (bmat->array==NULL||bmat->array_inv==NULL) {
    return 0;
  }
  return 1;
};
/* Get a pointer to a location in the band matrix, using
   the row and column indices of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indices out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}
/* Retrun the value of a location in the band matrix, using
   the row and column indices of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column);
}
/* Set the value of a location in the band matrix, using
   the row and column indices of the full matrix.           */
double setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
  return 0;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                         */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  int i,bandno;
  for(i=0;i<bmat->ncol;i++) {
    for (bandno=0;bandno<bmat->nbrows;bandno++) {
      bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];};
    x[i] = b[i];};
  long nrhs = 1;
  long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  return info;
}

int printmat(band_mat *bmat) {
  long j,p;
    for(j=0; j<bmat->ncol;j++) {
      for(p=0; p<bmat->nbrows; p++) {
	      printf("%ld %ld %g \n",j,p,bmat->array[bmat->nbrows*j + p]);};};
  return 0;
}
/* Return the 1D element index corresponding to a particular grid point, i */
long i_indx(long j, long p, long J, long P) {
  if (j <  0) j += J; /* Left Periodicity  */
  if (j >= J) j -= J; /* Right Periodicity */
  if (p <  0) p += P; /* Lower Periodicity */
  if (p >= P) p -= P; /* Upper Periodicity */
  return j*P + p;
}
/* Index for periodicity, m(i) */
long m_indx(long j, long p, long J, long P) {
  if (j <  0) j += J; /* Left Periodicity  */
  if (j >= J) j -= J; /* Right Periodicity */
  if (p <  0) p += P; /* Lower Periodicity */
  if (p >= P) p -= P; /* Upper Periodicity */
  long N = J*P, i = j*P + p;
  if(i<(double)N/2) return 2*i;
  else return 2*(N-i)-1;
}

int main(){
  //Parameters
  long i,m,j,p;    /*Indices for looping*/
  long ntheta;   /*Number of grid points in theta direction*/
  long nzeta;    /*Number of grid points in zeta direction*/
  long E;        /*If equal to 0, solve for steady state, else solve for time evolution*/
  double t_f;    /*Final time for time evolution mode*/
  long I_min;    /*Minimum number of timesteps (of equal length) in timestepping mode*/
  band_mat bmat; /*Banded matrix*/
  /*Read and scan in parameters*/
  FILE *input = fopen("input.txt","r");
  fscanf(input,"%ld %ld %ld %lf %ld",&ntheta,&nzeta,&E,&t_f, &I_min);
  fclose(input);
  long nbands_lower = 3*fmax(ntheta,nzeta);
  long nbands_upper = nbands_lower;
  long ncols = ntheta*nzeta; /*Total number of gridpoints*/
  double pi = 4*atan(1);
  double theta_len = 2*pi, zeta_len = 2*pi;
  double dtheta = theta_len/ntheta, dzeta = zeta_len/nzeta, dt = t_f/I_min;
  //Storage
  double *T =      calloc(ncols,sizeof(double));
  double *T_next = calloc(ncols,sizeof(double));
  double *T_temp = calloc(ncols,sizeof(double));
  double *b =      calloc(ncols,sizeof(double));
  double *Q11 =    calloc(ncols,sizeof(double));
  double *Q22 =    calloc(ncols,sizeof(double));
  double *H =      calloc(ncols,sizeof(double));
  double *S =      calloc(ncols,sizeof(double));
  double *R =      calloc(ncols,sizeof(double));
  //Initialisation
  init_band_mat(&bmat,nbands_lower,nbands_upper,ncols);
  FILE *coeff = fopen("coefficients.txt","r");
  for(j=0;j<ntheta;j++){
    for(p=0;p<nzeta;p++){
      i = i_indx(j,p,ntheta,nzeta);
      fscanf(coeff,"%lf %lf %lf %lf %lf",&Q11[i],&Q22[i],&H[i],&S[i],&R[i]);};};
  fclose(coeff);
  //STEADY STATE (Central Difference)
  if(E==0) {
    //Set band matrix values
    for(j=0;j<ntheta;j++) {
      for(p=0;p<nzeta;p++) {
        i = i_indx(j,p,ntheta,nzeta);
        m = m_indx(j,p,ntheta,nzeta);
        double A = Q11[i]/(dtheta*dtheta);
        double B = Q22[i]/(dzeta*dzeta);
        double C = (Q11[i_indx(j+1,p,ntheta,nzeta)]-Q11[i_indx(j-1,p,ntheta,nzeta)])/(4*dtheta*dtheta);
        double D = (Q22[i_indx(j,p+1,ntheta,nzeta)]-Q22[i_indx(j,p-1,ntheta,nzeta)])/(4*dzeta*dzeta);
        //j-1,p term
        setv(&bmat,m,m_indx(j-1,p,ntheta,nzeta),A-C);
        //j+1,p term
        setv(&bmat,m,m_indx(j+1,p,ntheta,nzeta),A+C);
        //j,p-1 term
        setv(&bmat,m,m_indx(j,p-1,ntheta,nzeta),B-D);
        //j,p+1 term
        setv(&bmat,m,m_indx(j,p+1,ntheta,nzeta),B+D);
        //j,p term
        setv(&bmat,m,m,-2*(A+B)-R[i]/H[i]);
        //Source term
        b[m] = -S[i]/H[i];
      }
    }
    //Solve matrix equation
    solve_Ax_eq_b(&bmat,T,b);
    //Print output to file
    FILE *outfile = fopen("output.txt","a");
    for(j=0; j<ntheta; j++) {
      for(p=0; p<nzeta; p++) {
        fprintf(outfile,"%lf %lf %lf \n",j*dtheta,p*dzeta,T[m_indx(j,p,ntheta,nzeta)]);};};
    fclose(outfile);
  }
  //TIME EVOLUTION (Central Differnce, Implicit Method)
  else {
    double ctime = 0; /* Current time */
    while (ctime<t_f) {
      //Set band matrix values
      for(j=0;j<ntheta;j++) {
        for(p=0;p<nzeta;p++) {
          i = i_indx(j,p,ntheta,nzeta);
          m = m_indx(j,p,ntheta,nzeta);
          double A = Q11[i]/(dtheta*dtheta);
          double B = Q22[i]/(dzeta*dzeta);
          double C = (Q11[i_indx(j+1,p,ntheta,nzeta)]-Q11[i_indx(j-1,p,ntheta,nzeta)])/(4*dtheta*dtheta);
          double D = (Q22[i_indx(j,p+1,ntheta,nzeta)]-Q22[i_indx(j,p-1,ntheta,nzeta)])/(4*dzeta*dzeta);
          //j-1,p term
          setv(&bmat,m,m_indx(j-1,p,ntheta,nzeta),A-C);
          //j+1,p term
          setv(&bmat,m,m_indx(j+1,p,ntheta,nzeta),A+C);
          //j,p-1 term
          setv(&bmat,m,m_indx(j,p-1,ntheta,nzeta),B-D);
          //j,p+1 term
	        setv(&bmat,m,m_indx(j,p+1,ntheta,nzeta),B+D);
          //j,p term
          setv(&bmat,m,m,-1/(dt*H[i])-(2*(A+B)+R[i]/H[i]));
          //RHS term
          b[m] = (-T[m]/dt-S[i])/H[i];
        }
      }
      //Solve matrix equation
      solve_Ax_eq_b(&bmat,T_next,b);
      //Increment time
      ctime += dt;
      //Swap values of pointers to copy arrays without element-by-element assignment
      T_temp = T_next;
      T_next = T;
      T = T_temp;
    }
    //Print output to file
    FILE *outfile = fopen("output.txt","a");
    for(j=0; j<ntheta; j++) {
      for(p=0; p<nzeta; p++) {
        fprintf(outfile,"%lf %lf %lf \n",j*dtheta,p*dzeta,T[m_indx(j,p,ntheta,nzeta)]);};};
    fclose(outfile);
  }
  //Free memory
  free(T);
  free(T_next);
  free(T_temp);
  free(b);
  free(Q11);
  free(Q22);
  free(H);
  free(S);
  free(R);
  free((&bmat)->array);
  free((&bmat)->array_inv);
  free((&bmat)->ipiv);
  return 0;
}
