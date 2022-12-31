#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>

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

int main(){
  //Parameters
  long ncols; /*Number of Grid Points*/
  double L;   /*Domain Length*/
  double v;   /*Advection Velocity*/
  double tau; /*Decay Rate*/
  FILE *input;
  input=fopen("input.txt","r");
  fscanf(input,"%lf %ld %lf %lf",&L,&ncols,&v,&tau);
  fclose(input);
  double dx = L/(ncols+1);
  band_mat bmat1,bmat2;
  long nbands_low=1,nbands_up=1;
  //Storage
  double *b1,*b2,*A,*B,*k,*S;
  b1=malloc(ncols*sizeof(double));
  b2=malloc(ncols*sizeof(double));
  A=malloc(ncols*sizeof(double));
  B=malloc(ncols*sizeof(double));
  k=malloc(ncols*sizeof(double));
  S=malloc(ncols*sizeof(double));
  //Initialisation
  init_band_mat(&bmat1,nbands_low,nbands_up,ncols);
  init_band_mat(&bmat2,nbands_low,nbands_up,ncols);
  FILE *coeff;
  coeff=fopen("coefficients.txt","r");
  long i;
  for(i=0;i<ncols;i++) {
    b1[i]=0;
    b2[i]=0;
    A[i]=0;
    B[i]=0;
    fscanf(coeff,"%lf %lf",&k[i],&S[i]);};
  fclose(coeff);
  //Set and solve for A (backwards difference)
  for(i=0;i<ncols;i++) {
    if(i>0)       {setv(&bmat1,i,i-1,k[i-1]+v*dx);
    setv(               &bmat1,i,i,-k[i]-k[i-1]-v*dx-tau*dx*dx);};
    if(i<ncols-1) {setv(&bmat1,i,i+1,k[i]);};
    b1[i]=-S[i]*dx*dx;};
  setv(&bmat1,0,0,k[0]-v*dx); /*boundary condition*/
  solve_Ax_eq_b(&bmat1,A,b1);
  //Set and solve for B (backwards difference)
  for(i=0;i<ncols;i++) {
    if(i>0)       {setv(&bmat2,i,i-1,k[i-1]+v*dx);
    setv(               &bmat2,i,i,-k[i]-k[i-1]-v*dx);};
    if(i<ncols-1) {setv(&bmat2,i,i+1,k[i]);};
    b2[i]=-A[i]*tau*dx*dx;};
  setv(&bmat2,0,0,k[0]-v*dx); /*boundary condition*/
  solve_Ax_eq_b(&bmat2,B,b2);
  //Print output
  FILE *outfile;
  outfile=fopen("output.txt","a");
  fprintf(outfile,"%lf %lf %lf \n",0.0,0.0,0.0); /*x=0 boundary condition*/
  for(i=0;i<ncols;i++) {
    fprintf(outfile,"%lf %lf %lf \n",(i+1)*dx,A[i],B[i]);};
  fprintf(outfile,"%lf %lf %lf \n",L,0.0,0.0);   /*x=L boundary condition*/
  fclose(outfile);
  //Free memory
  free(b1);
  free(b2);
  free(A);
  free(B);
  free(k);
  free(S);
  free((&bmat1)->array);
  free((&bmat1)->array_inv);
  free((&bmat1)->ipiv);
  free((&bmat2)->array);
  free((&bmat2)->array_inv);
  free((&bmat2)->ipiv);
  return 0;
}
