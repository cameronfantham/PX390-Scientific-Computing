#include <stdio.h>
int factors(long k){
    int b=1;
    int factors = 0;
    do{
        if(k%b==0.0){
            factors++;
        }
        b++;
    }while(b<k+1);
    printf("%d",factors);
    return 0;
}