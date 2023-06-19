//
// File          _______ : io_fns.c
// Description   _______ : I/O functions for creating interpolation tables.
//                         These functions work only with float types and
//                         are intended to be use as a Fortran-C interface
//                         only by the Fortran programs in this directory.
// Project       _______ : DeSIRe
// Creation date _______ : 12/12/22
// Author        _______ : epm@iac.es
//


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static FILE **fpfile = NULL;
static int    nfiles = 0;


//____________________________________________________________________________
//
//  Method: openfile_()
//
/** Open a file.
 *  @param  ncodes size of 'codes'.
 *  @param  codes array with the ASCII codes for the filename.
 *  @return file descriptor index (-1 if error).
 */
//____________________________________________________________________________

int openfile_( int *ncodes, int *codes )
{
   char filename[*ncodes + 1];
   int  i;

   // Get the filename from the ASCII codes.
   for (i = 0; i < *ncodes; i++)
   {
      filename[i] = (char) codes[i];
   }
   filename[i] = '\0';
   
   // Open the file.
   fpfile = realloc(fpfile, (nfiles+1) * sizeof(FILE *));
   if ((fpfile[nfiles] = fopen(filename, "wb")) == NULL)
   {
      return(-1);
   }
   nfiles++;

   return(nfiles-1);
}


//____________________________________________________________________________
//
//  Method: writefile_()
//
/** Write the file.
 *  @param  ifile file descriptor index.
 *  @param  nvalues number of values to write.
 *  @param  values array with the values to write.
 *  @return number of items written (-1 if error).
 */
//____________________________________________________________________________

int writefile_( int *ifile, int *nvalues, float *values )
{
   int n;

   if (fpfile[*ifile] == NULL)
   {
      return(-1);
   }
   n = fwrite(values, sizeof(float), *nvalues, fpfile[*ifile]);

   return(n);
}


//____________________________________________________________________________
//
//  Method: closefile_()
//
/** Close the file.
 *  @param  ifile file descriptor index.
 *  @return error code (0 if successful).
 */
//____________________________________________________________________________

int closefile_( int *ifile )
{
   int ier;
   
   if (fpfile[*ifile] == NULL)
   {
      return(0);
   }
   ier = fclose(fpfile[*ifile]);
   fpfile[*ifile] = NULL;

   return(ier);
}


//____________________________________________________________________________
//
//  Method: main()
//
/** Main method.
 *  @param  argc number of arguments passed to the program.
 *  @param  argv array of strings with the arguments passed.
 *  @return exit code (0 if successful).
 *
 *  fopen() vs open()
 *  There are four main reasons to use fopen instead of open.
 *  1) fopen provides you with buffering I/O that may turn out to be a lot
 *     faster than what you're doing with open.
 *  2) fopen does line ending translation if the file is not opened in binary
 *     mode, which can be very helpful if your program is ever ported to a
 *     non-Unix environment.
 *  3) A FILE * gives you the ability to use fscanf and other stdio functions.
 *  4) Portability: use fopen() when you want to write portable code.
 *     The open() function is not defined in the C standard.
 *     If you're doing seeks (aka fsetpos or fseek the second of which is
 *     slightly trickier to use in a standards compliant way), the usefulness
 *     of buffering quickly goes down.
 */
//____________________________________________________________________________

int xmain( int argc, char *argv[] )
{
   char  *ptro;
   float  values[3][4] = {{100, 110, 120, 130},
                          {200, 210, 220, 230},
                          {300, 310, 320, 330}};
   float  info[6] = {4, 10, 1, 3, 20, 2}; // cols xini xdelta rows yini ydelta
   float  nboxes  = 1;
   float  xyscan  = 1;
   int    ncols   = 4;
   int    nrows   = 3;
   int    ninfo   = 6;
   int    n       = 1;
   int    codes[64], ncodes;
   int    ifile;
   int    i;

   // Command line.
   if (argc < 2)
   {
      fprintf(stderr, "Usage:> %s <filename>\n", argv[0]);
      return(-1);
   }

   // Filename in ASCII codes.
   ptro = argv[1];
   ncodes = strlen(argv[1]);
   for (i = 0; i < ncodes; i++, ptro++)
   {
      codes[i] = (int) *ptro;
   }
   codes[i] = '\0';

   // Create the table.
   if (ifile = openfile_(&ncodes, codes) == -1)
   {
      fprintf(stderr, "Error: Unable to create the file '%s'\n", argv[1]);
      return(-1);
   }

   // Write the table information.
   if (writefile_(&ifile, &n, &nboxes) != n)
   {
      fprintf(stderr, "Error: Unable to write the table on disk\n");
      return(-1);
   }
   if (writefile_(&ifile, &n, &xyscan) != n)
   {
      fprintf(stderr, "Error: Unable to write the table on disk\n");
      return(-1);
   }
   if (writefile_(&ifile, &ninfo, info) != ninfo)
   {
      fprintf(stderr, "Error: Unable to write the table on disk\n");
      return(-1);
   }

   // Write the table.
   for (i = 0; i < nrows; i++)
   {
      if (writefile_(&ifile, &ncols, values[i]) != ncols)
      {
         fprintf(stderr, "Error: Unable to write the table on disk\n");
         return(-1);
      }
   }

   // Close the file.
   if (closefile_(&ifile) != 0)
   {
      fprintf(stderr, "Error: Unable to close the table\n");
      return(-1);
   }

   return(0);
}


//____________________________________________________________________________
