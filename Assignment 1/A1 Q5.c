#include <stdio.h>
int order_of_magnitude(double a){
    int c=0;
    if(a>=1){
        do{
            a=a/10;
            c++;
        }while(a>=10);
    }
    if(a<1){
        do{
            a=a*10;
            c--;
        }while(a<1);
    }
    return c;
}
