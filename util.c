
// System includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include "mpi.h"

// User includes
//#include "globals.h"
//#include "pprintf.h"

FILE *fp;

void cgetheader_(char *filename, int *width, int *height, int *depth)
{
  fp = fopen( filename, "r");
  if( !fp ) {
      printf( "Error: The file '%s' could not be opened.\n", filename );
      return;
  }
  char header[10];
  int bwidth, bheight, bdepth;
  int rv = fscanf( fp, "%6s\n%i %i\n%i\n", header, &bwidth, &bheight, &bdepth );
/*   printf("%s:%s %i %i %i\n", filename, header, bwidth, bheight, bdepth); */

  *width = bwidth;
  *height = bheight;
  *depth = bdepth;
}

void cgetheader(char *filename, int *width, int *height, int *depth)
{
  cgetheader_(filename, width, height, depth);
}

void cfgetc_(int *intval)
{
  *intval = fgetc(fp);
}

void cfgetc(int *intval)
{
  cfgetc_(intval);
}

void cfclose_()
{
  fclose(fp);
}

void cfclose()
{
  cfclose_();
}

