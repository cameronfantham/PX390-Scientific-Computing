#include <stdio.h>
#include <stdlib.h>
void  concat_arrays(long **newarrptr, long *arr1, int nelements1, long *arr2, int nelements2){
  long *arr = malloc(sizeof(int)*4*(nelements1+nelements2));
  int i,j, k;
  for(i=0;i<nelements1;i++)
    arr[i]=arr1[i];
  for(j=nelements1,i=0;i<nelements2;i++,j++)
    arr[j]=arr2[i];
  *newarrptr=arr;
}
 int main(){
   long *newarrptr;
   long arr1[] = {1,2,3};
   long arr2[] = {4,5,6};

   concat_arrays(&newarrptr, arr1, 3, arr2, 3);
   for(int i=0;i<6;i++)
      printf("%lo\n",*(newarrptr+i));
 }