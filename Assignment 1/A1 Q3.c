#include <stdio.h>
#include <string.h>
void ascii_box(int size){
    int i=1;
    int j=1;
    int k=1;
    scanf("%d", &size);
    char src1[50], src2[50], dest1[50], dest2[50];
    strcpy(src1,  "X");
    strcpy(src2, "O");
    strcpy(dest1, "X");
    strcpy(dest2, "X");
    do{
        strcat(dest1, src1);
        i++;
    }while(i<size);
    printf("%s\n", dest1);
    do{
        strcat(dest2, src2);
        j++;
    }while(j<size-1);
    strcat(dest2, src1);
    do{
        printf("%s\n", dest2);
        k++;
    }while(k<size-1);
    printf("%s\n", dest1);
    return 0;
}