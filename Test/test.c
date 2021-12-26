#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define min(x, y) (((x) < (y)) ? (x) : (y))

int main(){
    printf("%d", min(5, 10));
    return 0;
}