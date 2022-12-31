#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>
#include <math.h>

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
  long j,p;
    for(j=0; j<bmat->ncol;j++) {
      for(p=0; p<bmat->nbrows; p++) {
	printf("%ld %ld %g \n",j,p,bmat->array[bmat->nbrows*j + p]);
      }
    }
  return 0;
}
/*Check that a grid point has valid coordinates */
int is_valid(long j, long p, long J, long P) {
  return (j>=0)&&(j<J)&&(p>=0)&&(p<P);
}
/* Index */
long indx(long j, long p, long J, long P){
  return j*P + p;
}
int main(){
  //Parameters
  long i,j,p;
  long ntheta; /*Number of grid points in theta direction*/
  long nzeta;  /*Number of grid points in zeta direction*/
  long E;  /*If equal to 0, solve for steady state. Else solve for time evolution*/
  long I_min;  /* Minimum number of timesteps (of equal length) in timestepping mode.*/
  double t_f;  /*Final time for time evolution mode*/
  FILE *input;
  input=fopen("input.txt","r");
  fscanf(input,"%ld %ld %ld %lf %ld",&ntheta,&nzeta,&E,&t_f, &I_min);
  fclose(input);
  band_mat bmat;
  long nbands_low = 2*nzeta-1;
  long nbands_up = nbands_low;
  long ncols = ntheta*nzeta;
  double pi = 4*atan(1);
  double theta_len = 2*pi, zeta_len = 2*pi;
  double dtheta = theta_len/ntheta;
  double dzeta = zeta_len/nzeta;
  double dt = t_f/I_min;
  //Storage
  double *T,*T_next,*T_temp,*b,*Q11,*Q22,*H,*S,*R;
  T =      malloc(ncols*sizeof(double));
  T_next = malloc(ncols*sizeof(double));
  T_temp = malloc(ncols*sizeof(double));
  b =      malloc(ncols*sizeof(double));
  Q11 =    malloc(ncols*sizeof(double));
  Q22 =    malloc(ncols*sizeof(double));
  H =      malloc(ncols*sizeof(double));
  S =      malloc(ncols*sizeof(double));
  R =      malloc(ncols*sizeof(double));
  //Initialisation
  init_band_mat(&bmat,nbands_low,nbands_up,ncols);
  FILE *coeff;
  coeff=fopen("coefficients.txt","r");
  for(j=0;j<ntheta;j++){
    for(p=0;p<nzeta;p++){
      i = indx(j,p,ntheta,nzeta);
      T[i]=0;
      T_next[i]=0;
      b[i]=0;
      fscanf(coeff,"%lf %lf %lf %lf %lf",&Q11[i],&Q22[i],&H[i],&S[i],&R[i]);};};
  fclose(coeff);
  //Steady state
  if(E==0) {
    long j0 = -ntheta/2, p0 = -nzeta/2;
    for(j=0;j<ntheta;j++){
      for(p=0;p<nzeta;p++){
        double A,B,C,D;
        A = Q11[indx(j,p,ntheta,nzeta)]/(dtheta*dtheta);
        B = Q22[indx(j,p,ntheta,nzeta)]/(dzeta*dzeta);
        if(j!=j0||p!=p0){
          //j-1,p
          if(is_valid(j-1,p,ntheta,nzeta)) {
	    setv(&bmat,indx(j,p,ntheta,nzeta),indx(j-1,p,ntheta,nzeta),A);
          } else {  /* Periodicity */
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j-1+ntheta,p,ntheta,nzeta),A);};
          //j+1,p
          if(is_valid(j+1,p,ntheta,nzeta)) {
            C = (Q11[indx(j+1,p,ntheta,nzeta)]-Q11[indx(j,p,ntheta,nzeta)])/(dtheta*dtheta);
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j+1,p,ntheta,nzeta),C+A);
          } else {  /* Periodicity */
            C = (Q11[indx((j+1)-ntheta,p,ntheta,nzeta)]-Q11[indx(j,p,ntheta,nzeta)])/(dtheta*dtheta);
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j+1-ntheta,p,ntheta,nzeta),C+A);};
          //j,p-1
          if(is_valid(j,p-1,ntheta,nzeta)) {
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p-1,ntheta,nzeta),B);
          } else {  /* Periodicity */
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p-1+nzeta,ntheta,nzeta),B);};
          //j,p+1
          if(is_valid(j,p+1,ntheta,nzeta)) {
            D = (Q22[indx(j,p+1,ntheta,nzeta)]-Q22[indx(j,p,ntheta,nzeta)])/(dzeta*dzeta);
	    setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p+1,ntheta,nzeta),D+B);
          } else {  /* Periodicity */
            D = (Q22[indx(j,(p+1)-nzeta,ntheta,nzeta)]-Q22[indx(j,p,ntheta,nzeta)])/(dzeta*dzeta);
	    setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p+1-nzeta,ntheta,nzeta),D+B);};
          //j,p
          if(is_valid(j+1,p,ntheta,nzeta)&&is_valid(j,p+1,ntheta,nzeta)){
            C = (Q11[indx(j+1,p,ntheta,nzeta)]-Q11[indx(j,p,ntheta,nzeta)])/(dtheta*dtheta);
            D = (Q22[indx(j,p+1,ntheta,nzeta)]-Q22[indx(j,p,ntheta,nzeta)])/(dzeta*dzeta); 
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p,ntheta,nzeta),-2*A-2*B-C-D-R[indx(j,p,ntheta,nzeta)]/H[indx(j,p,ntheta,nzeta)]);
          } if(!is_valid(j+1,p,ntheta,nzeta)&&is_valid(j,p+1,ntheta,nzeta)) {  /* Periodicity */
            C = (Q11[indx(j+1-ntheta,p,ntheta,nzeta)]-Q11[indx(j,p,ntheta,nzeta)])/(dtheta*dtheta);
            D = (Q22[indx(j,p+1,ntheta,nzeta)]-Q22[indx(j,p,ntheta,nzeta)])/(dzeta*dzeta); 
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p,ntheta,nzeta),-2*A-2*B-C-D-R[indx(j,p,ntheta,nzeta)]/H[indx(j,p,ntheta,nzeta)]);
          } if(is_valid(j+1,p,ntheta,nzeta)&&!is_valid(j,p+1,ntheta,nzeta)) {  /* Periodicity */
            C = (Q11[indx(j+1,p,ntheta,nzeta)]-Q11[indx(j,p,ntheta,nzeta)])/(dtheta*dtheta);
            D = (Q22[indx(j,p+1-nzeta,ntheta,nzeta)]-Q22[indx(j,p,ntheta,nzeta)])/(dzeta*dzeta); 
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p,ntheta,nzeta),-2*A-2*B-C-D-R[indx(j,p,ntheta,nzeta)]/H[indx(j,p,ntheta,nzeta)]);};
          /* Uniform source term in heat equation */
          b[indx(j,p,ntheta,nzeta)] = -S[indx(j,p,ntheta,nzeta)]/H[indx(j,p,ntheta,nzeta)];
        }
      }
    }
    solve_Ax_eq_b(&bmat,T,b);
    //Print output
    FILE *outfile;
    outfile=fopen("output.txt","a");
    for(j=0; j<ntheta; j++) {
      for(p=0; p<nzeta; p++) {   
        fprintf(outfile,"%lf %lf %lf \n",j*dtheta,p*dzeta,T[indx(j,p,ntheta,nzeta)]);};};
    fclose(outfile);
  }
  //Time evolution
  else {
    double ctime = 0; /* Current time */
    while (ctime<t_f){
      for(j=0;j<ntheta;j++){
        for(p=0;p<nzeta;p++){
          double A,B,C,D;
          A = Q11[indx(j,p,ntheta,nzeta)]*H[indx(j,p,ntheta,nzeta)]/(dtheta*dtheta);
          B = Q22[indx(j,p,ntheta,nzeta)]*H[indx(j,p,ntheta,nzeta)]/(dzeta*dzeta);
          //j-1,p
          if(is_valid(j-1,p,ntheta,nzeta)) {
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j-1,p,ntheta,nzeta),A*dt);
          } else {  /* Periodicity */
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j-1+ntheta,p,ntheta,nzeta),A*dt);};
          //j+1,p
          if(is_valid(j+1,p,ntheta,nzeta)) {
            C = (Q11[indx(j+1,p,ntheta,nzeta)]-Q11[indx(j,p,ntheta,nzeta)])*H[indx(j,p,ntheta,nzeta)]/(dtheta*dtheta);
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j+1,p,ntheta,nzeta),(C+A)*dt);
          } else {  /* Periodicity */
            C = (Q11[indx(j+1-ntheta,p,ntheta,nzeta)]-Q11[indx(j,p,ntheta,nzeta)])*H[indx(j,p,ntheta,nzeta)]/(dtheta*dtheta);
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j+1-ntheta,p,ntheta,nzeta),(C+A)*dt);};
          //j,p-1
          if(is_valid(j,p-1,ntheta,nzeta)) {
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p-1,ntheta,nzeta),B*dt);
          } else {  /* Periodicity */
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p-1+nzeta,ntheta,nzeta),B*dt);};
          //j,p+1
          if(is_valid(j,p+1,ntheta,nzeta)) {
            D = (Q22[indx(j,p+1,ntheta,nzeta)]-Q22[indx(j,p,ntheta,nzeta)])*H[indx(j,p,ntheta,nzeta)]/(dzeta*dzeta);
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p+1,ntheta,nzeta),(D+B)*dt);
          } else {  /* Periodicity */
            D = (Q22[indx(j,p+1-nzeta,ntheta,nzeta)]-Q22[indx(j,p,ntheta,nzeta)])*H[indx(j,p,ntheta,nzeta)]/(dzeta*dzeta);
	    setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p+1-nzeta,ntheta,nzeta),(D+B)*dt);};
          //j,p
          if(is_valid(j+1,p,ntheta,nzeta)&&is_valid(j,p+1,ntheta,nzeta)){
            C = (Q11[indx(j+1,p,ntheta,nzeta)]-Q11[indx(j,p,ntheta,nzeta)])*H[indx(j,p,ntheta,nzeta)]/(dtheta*dtheta);
            D = (Q22[indx(j,p+1,ntheta,nzeta)]-Q22[indx(j,p,ntheta,nzeta)])*H[indx(j,p,ntheta,nzeta)]/(dzeta*dzeta); 
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p,ntheta,nzeta),1+(-2*A-2*B-C-D-R[indx(j,p,ntheta,nzeta)])*dt);
          } if(!is_valid(j+1,p,ntheta,nzeta)&&is_valid(j,p+1,ntheta,nzeta)) {  /* Periodicity */
            C = (Q11[indx(j+1-ntheta,p,ntheta,nzeta)]-Q11[indx(j,p,ntheta,nzeta)])*H[indx(j,p,ntheta,nzeta)]/(dtheta*dtheta);
            D = (Q22[indx(j,p+1,ntheta,nzeta)]-Q22[indx(j,p,ntheta,nzeta)])*H[indx(j,p,ntheta,nzeta)]/(dzeta*dzeta); 
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p,ntheta,nzeta),1+(-2*A-2*B-C-D-R[indx(j,p,ntheta,nzeta)])*dt);
          } if(is_valid(j+1,p,ntheta,nzeta)&&!is_valid(j,p+1,ntheta,nzeta)) {  /* Periodicity */
            C = (Q11[indx(j+1,p,ntheta,nzeta)]-Q11[indx(j,p,ntheta,nzeta)])*H[indx(j,p,ntheta,nzeta)]/(dtheta*dtheta);
            D = (Q22[indx(j,p+1-nzeta,ntheta,nzeta)]-Q22[indx(j,p,ntheta,nzeta)])*H[indx(j,p,ntheta,nzeta)]/(dzeta*dzeta); 
            setv(&bmat,indx(j,p,ntheta,nzeta),indx(j,p,ntheta,nzeta),1+(-2*A-2*B-C-D-R[indx(j,p,ntheta,nzeta)])*dt);};
            /* Uniform source term in heat equation */
            b[indx(j,p,ntheta,nzeta)] = T_next[indx(j,p,ntheta,nzeta)]-S[indx(j,p,ntheta,nzeta)]*dt;   
          }
        }
      }
      solve_Ax_eq_b(&bmat,T,b);
      // Swap values of pointers to copy arrays without element-by-element assignment 
      T_temp = T;
      T = T_next;
      T_next = T_temp;
      // Increment time
      ctime += dt;   
      // Output values
      FILE *outfile;
      outfile=fopen("output.txt","a");
      for(j=0; j<ntheta; j++) {
        for(p=0; p<nzeta; p++) {   
          fprintf(outfile,"%lf %lf %lf %lf \n",ctime,j*dtheta,p*dzeta,T[indx(j,p,ntheta,nzeta)]);};};
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