#include <stdio.h>
#include <math.h>
void en_pot(double *posx, double *posy, double *posz, 
long ncharges, double *res){
    int i,j;
    double a,b,c;
    double ep=0;
    for(i=0;i<ncharges;i++)
        for(j=0;j<ncharges;j++)
            if(i!=j){
               a=pow(posx[i]-posx[j],2);
               b=pow(posy[i]-posy[j],2);
               c=pow(posz[i]-posz[j],2);
               ep+=pow(a+b+c,-0.5);
            }       
    *res=0.5*ep;
}

int main(){
  double d = 0;
  double *arr = &d;
  double posx[50] = {0,0};
  double posy[50] = {0.001,0};  
  double posz[50] = {0,0};
  printf("%lf\n%p\n", *arr, arr);
  en_pot(posx,posy,posz,2, arr);
  printf("%lf\n%p\n", *arr, arr);
}
