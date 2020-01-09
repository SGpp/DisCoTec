//  Returns the timestamp of a file.  By khr, 05/2009

#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

int main(int argc, char **argv)
{

  struct stat buf;
  struct tm  *ts;
  int exists;
  char timestr[80];

  if (argc<2) {
    fprintf(stderr, "Argument <filename> is missing.\n");
    return(1);
  }

  exists = stat(argv[1], &buf);

  if (exists < 0) {

    fprintf(stderr, "%s not found\n", argv[1]);
    return(1);

  } else {

    ts = localtime(&buf.st_mtime);

    // Readable format:
    // strftime(timestr, sizeof(timestr), "%a %Y-%m-%d %H:%M:%S %Z", ts);

    // Compressed format:
    strftime(timestr, sizeof(timestr), "%Y%m%d%H%M%S", ts);

    printf("%s\n", timestr);

  }

  return(0);

}

