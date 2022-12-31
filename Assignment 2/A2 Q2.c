#include <stdio.h>
int arr_product(int *arr, int len){ 
    int result = 1; 
    for (int i = 0; i < len; i++) 
        result = result * arr[i]; 
    return result; 
}