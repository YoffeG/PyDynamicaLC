#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

void int2str(long n, char *str)
{
    int len = strlen(str);
    int i;
    for( i = 0; i < len; i++ )
    {
        long mask = pow(10, i+1); /* マスクする値。12345を10で割れば5が得れるし100で割れば45が得れる */
        int a = (n % mask) / pow(10, i); /* マスクしたものが45なら10で割って4を、345なら100で割って3を得る */
        char c = (char)(a + (unsigned short)'0'); /* '0'..'9'は連番なので、 N+'0'で'N'になる */
        str[len-i-1] = c; /* 下の桁から順番に代入する */
    }
}

char* concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

int main (int argc, char *argv[])
{
    long i = 12345;
    char str[] = "00000";
    int2str(i, str);
    puts(str);
    
    char* s = concat(str, "herp");
    puts(s);
    free(s); // deallocate the string
    return(0);
}
