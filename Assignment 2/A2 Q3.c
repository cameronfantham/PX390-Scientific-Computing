#include <stdio.h>
void remvowels(char *string){
    int count = 0;
    for (int i = 0; string[i]; i++) {
        if (string[i] != 'a' && string[i] != 'e' && string[i] != 'i'
        &&  string[i] != 'o' && string[i] != 'u'){
            string[count++] = string[i];
        }
    }
    string[count] = '\0';
}

int main(){
  char str[]="ba ba ba ba ba boran";
  remvowels(str);
  printf("S:%s",str);
}