#include <dirent.h>
#include <string.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
struct dirent ** dirEntry;
int n = scandir(".", &dirEntry, NULL, alphasort);
printf("n=%d\n", n);
int i = 0;

     while (i<n)
	{
        	puts ((*dirEntry)->d_name);
      		dirEntry++;
		i++;
	}


return 0;
}
