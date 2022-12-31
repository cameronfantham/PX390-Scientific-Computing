#include <stdio.h>
#include <stdlib.h>
int writenumstofile(int num1, int num2){
    /*open file*/
    FILE *f = fopen("output.txt", "w");
    /*display error message*/
    if (f == NULL){
        return 1;
        exit(1);
    }
    /* print integers */
    scanf("%d%d",&num1,&num2);
    fprintf(f, "%d %d", num1, num2);
    /*close file*/
    fclose(f);
    return 0;
}