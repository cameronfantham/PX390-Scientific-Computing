#include <stdio.h>
int  *reverse_order(int *arr, int arr_length){
  int i, j;
  int rev_arr[50];
  for (i=arr_length-1, j=0; i>=0; i--, j++)
      rev_arr[j] = arr[i];
  printf("%d\n",rev_arr[2]);
  *arr=*rev_arr;    
  return arr;
}

int main(){
  int arrl=3;
  int arr[]={3,4,-4};
  int *rarr = reverse_order(arr,arrl);
  int i;
  for(i=0;i<arrl;i++)
      printf("%d\n",*(rarr+i));
  return 0;
}