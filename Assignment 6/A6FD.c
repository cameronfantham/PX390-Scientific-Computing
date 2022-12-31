#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>
#include <math.h>
#include <string.h>

struct band_mat{
  long ncol;        /* Number of columns in band matrix         */
  long nbrows;      /* Number of rows (bands in original matrix)*/
  long nbands_up;   /* Number of bands above diagonal           */
  long nbands_low;  /* Number of bands below diagonal           */
  double *array;    /* Storage for the matrix in banded format  */
  /* Internal temporary storage for solving inverse problem     */
  long nbrows_inv;  /* Number of rows of inverse matrix         */
  double *array_inv;/* Store the inverse if this is generated   */
  int *ipiv;        /* Additional inverse information           */
};
typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
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
  for (i=0;i<(bmat->nbrows+bmat->nbands_low)*bmat->ncol;i++) {
	bmat->array_inv[i] = 0.0;
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

/* Index for periodicity */
long indx(long j, long p, long J, long P) {
  if (j < 0) j += J;  /* Left Periodicity  */
  if (j >= J) j -= J; /* Right Periodicity */
  if (p < 0) p += P;  /* Lower Periodicity */
  if (p >= P) p -= P; /* Upper Periodicity */
  long N = J*P, i = j*P + p;
  if(i<(double)N/2) return 2*i;
  else return 2*(N-i)-1;
}
/* Reindex */
long reindx(long j, long p, long J, long P){
  if (j < 0) j += J;  /* Left Periodicity  */
  if (j >= J) j -= J; /* Right Periodicity */
  if (p < 0) p += P;  /* Lower Periodicity */
  if (p >= P) p -= P; /* Upper Periodicity */
  return j*P + p;
}
int main(){
  //Parameters
  long i,j,p;  /*Indices*/
  long ntheta; /*Number of grid points in theta direction*/
  long nzeta;  /*Number of grid points in zeta direction*/
  long E;      /*If equal to 0, solve for steady state. Else solve for time evolution*/
  double t_f;  /*Final time for time evolution mode*/
  long I_min;  /* Minimum number of timesteps (of equal length) in timestepping mode.*/
  FILE *input = fopen("input.txt","r");
  fscanf(input,"%ld %ld %ld %lf %ld",&ntheta,&nzeta,&E,&t_f, &I_min);
  fclose(input);
  band_mat bmat;
  long nbands_low = 2*nzeta-1;
  long nbands_up = nbands_low;
  long ncols = ntheta*nzeta;
  double pi = 4*atan(1);
  double theta_len = 2*pi, zeta_len = 2*pi;
  double dtheta = theta_len/ntheta, dzeta = zeta_len/nzeta, dt;
  //Storage
  double *T =      malloc(sizeof(double)*ncols);
  double *T_next = malloc(sizeof(double)*ncols);
  double *T_temp = malloc(sizeof(double)*ncols);
  double *b =      malloc(sizeof(double)*ncols);
  double *Q11 =    malloc(sizeof(double)*ncols);
  double *Q22 =    malloc(sizeof(double)*ncols);
  double *H =      malloc(sizeof(double)*ncols);
  double *S =      malloc(sizeof(double)*ncols);
  double *R =      malloc(sizeof(double)*ncols);
  //Initialisation
  init_band_mat(&bmat,nbands_low,nbands_up,ncols);
  memset(T, 0.0, sizeof(double)*ncols);
  memset(T_next, 0.0, sizeof(double)*ncols);
  memset(T_temp, 0.0, sizeof(double)*ncols);
  memset(b, 0.0, sizeof(double)*ncols);
  FILE *coeff = fopen("coefficients.txt","r");
  for(j=0;j<ntheta;j++){
    for(p=0;p<nzeta;p++){
      i = reindx(j,p,ntheta,nzeta);
      fscanf(coeff,"%lf %lf %lf %lf %lf",&Q11[i],&Q22[i],&H[i],&S[i],&R[i]);};};
  fclose(coeff);
  //STEADY STATE
  if(E==0) {
    for(j=0;j<ntheta;j++) {
      for(p=0;p<nzeta;p++) {
        double A = Q11[reindx(j,p,ntheta,nzeta)]/(dtheta*dtheta);
        double B = Q22[reindx(j,p,ntheta,nzeta)]/(dzeta*dzeta);
        double C = (Q11[reindx(j+1,p,ntheta,nzeta)]-Q11[reindx(j,p,ntheta,nzeta)])/(dtheta*dtheta);
        double D = (Q22[reindx(j,p+1,ntheta,nzeta)]-Q22[reindx(j,p,ntheta,nzeta)])/(dzeta*dzeta);
        //j-1,p
        setv(&bmat,indx(j,p,ntheta,nzeta),indx(j-1,p,ntheta,nzeta),A);
        //j+1,p
        setv(&bmat,indx(j,p,ntheta,nzeta),indx(j+1,p,ntheta,nzeta),A+C);
        //j,p-1
        setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p-1,ntheta,nzeta),B);
        //j,p+1
	      setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p+1,ntheta,nzeta),B+D);
        //j,p
        setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p,ntheta,nzeta),-2*(A+B)-(C+D)-R[reindx(j,p,ntheta,nzeta)]/H[reindx(j,p,ntheta,nzeta)]);
        //Source term
        b[reindx(j,p,ntheta,nzeta)] = -S[reindx(j,p,ntheta,nzeta)]/H[reindx(j,p,ntheta,nzeta)];
      }
    }
    solve_Ax_eq_b(&bmat,T,b);
    //Print output
    FILE *outfile = fopen("output.txt","a");
    for(j=0; j<ntheta; j++) {
      for(p=0; p<nzeta; p++) {
        fprintf(outfile,"%lf %lf %lf \n",j*dtheta,p*dzeta,T[reindx(j,p,ntheta,nzeta)]);};};
    fclose(outfile);
  }
  //TIME EVOLUTION
  else {
    double ctime = 0; /* Current time */
    while (ctime<t_f) {
      for(j=0;j<ntheta;j++) {
        for(p=0;p<nzeta;p++) {
          double A = Q11[indx(j,p,ntheta,nzeta)]/(dtheta*dtheta);
          double B = Q22[indx(j,p,ntheta,nzeta)]/(dzeta*dzeta);
          double C = (Q11[indx(j+1,p,ntheta,nzeta)]-Q11[indx(j,p,ntheta,nzeta)])/(dtheta*dtheta);
          double D = (Q22[indx(j,p+1,ntheta,nzeta)]-Q22[indx(j,p,ntheta,nzeta)])/(dzeta*dzeta);
          //j-1,p
          setv(&bmat,indx(j,p,ntheta,nzeta),indx(j-1,p,ntheta,nzeta),A);
          //j+1,p
          setv(&bmat,indx(j,p,ntheta,nzeta),indx(j+1,p,ntheta,nzeta),A+C);
          //j,p-1
          setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p-1,ntheta,nzeta),B);
          //j,p+1
	        setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p+1,ntheta,nzeta),B+D);
          //j,p
          dt = 0.001;
          setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p,ntheta,nzeta),-1/(dt*H[indx(j,p,ntheta,nzeta)])-(2*(A+B)-(C+D)+R[indx(j,p,ntheta,nzeta)]/H[indx(j,p,ntheta,nzeta)]));
          //Source term
          b[indx(j,p,ntheta,nzeta)] = (-T[indx(j,p,ntheta,nzeta)]/dt-S[indx(j,p,ntheta,nzeta)])/H[indx(j,p,ntheta,nzeta)];
        }
      }
      solve_Ax_eq_b(&bmat,T_next,b);
      //Swap values of pointers to copy arrays without element-by-element assignment
      T_temp = T_next;
      T_next = T;
      T = T_temp;
      //Increment time
      ctime += dt;
    }
    //Output values
    FILE *outfile;
    outfile=fopen("output.txt","a");
    for(j=0; j<ntheta; j++) {
      for(p=0; p<nzeta; p++) {
        fprintf(outfile,"%lf %lf %lf %lf \n",t_f,j*dtheta,p*dzeta,T[indx(j,p,ntheta,nzeta)]);};};
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
