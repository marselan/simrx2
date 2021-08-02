#include <regex.h>
#include <stdio.h>


/*
 * Match string against the extended regular expression in
 * pattern, treating errors as no match.
 *
 * Return 1 for match, 0 for no match.
 */


int
match(const char *string, char *pattern)
{
    int    status;
    regex_t    re;


    if (regcomp(&re, pattern, REG_EXTENDED|REG_NOSUB) != 0) {
        return(0);      /* Report error. */
    }
    status = regexec(&re, string, (size_t) 0, NULL, 0);
    regfree(&re);
    if (status != 0) {
        return(0);      /* Report error. */
    }
    return(1);

}

int main(int argc, char *argv[])
{
int res = match(argv[1], argv[2]);
if(res) printf("MATCH\n");
else printf("NO MATCH\n");
return 0;
}
