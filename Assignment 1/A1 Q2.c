#include<stdio.h>
#include<math.h>
int powers(int a,int b){
    long result = powl(a,b);
    printf("%li",result);
    return 0;
}