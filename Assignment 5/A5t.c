#include <stdlib.h>
#include <stdio.h>

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
  return 1;
};


/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indices out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column);
}
/* Set the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
  return 0;
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
  band_mat bmat1;
  band_mat bmat2;
  long ncols = 10;
  double dx = 1/(ncols+1);
  long nbands_low = 1;  
  long nbands_up  = 1;
  init_band_mat(&bmat1, nbands_low, nbands_up, ncols);
  init_band_mat(&bmat2, nbands_low, nbands_up, ncols);
  long i;
  /* Loop over the equation number and set the matrices'
     values equal to the coefficients of the grid values 
     note boundaries treated with special cases           */
  for(i=0; i<ncols; i++) {
    if(i>0)       {setv(&bmat1,i,i-1,-1.0/(dx*dx));};
    setv(               &bmat1,i,i,   2.0/(dx*dx));
    if(i<ncols-1) {setv(&bmat1,i,i+1,-1.0/(dx*dx));}; 
  }
  for(i=0; i<ncols; i++) {
    if(i>0)       {setv(&bmat2,i,i-1,-1.0/(dx));};
    setv(               &bmat2,i,i,   0.0/(dx));
    if(i<ncols-1) {setv(&bmat2,i,i+1,-1.0/(dx));}; 
  }
  printmat(&bmat1);
  printmat(&bmat2);
  return 0;
}