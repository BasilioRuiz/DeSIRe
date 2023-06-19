//
// File          _______ : ptable.c
// Description   _______ : Program for printing interpolation tables.
// Project       _______ : DeSIRe
// Creation date _______ : 12/12/22
// Author        _______ : epm@iac.es
//


#include <stdio.h>
#include <stdlib.h>

#define NINFO 6



//____________________________________________________________________________
//
//  Method: main()
//
/** Main method.
 *  @param  argc number of arguments passed to the program.
 *  @param  argv array of strings with the arguments passed.
 *  @return exit code (0 if successful).
 */
//____________________________________________________________________________

int main( int argc, char *argv[] )
{
   FILE  *fp;
   float  (*info)[NINFO];  // bidimensional matrix with #columns fixed = NINFO
   float  num;
   float  xini, xdelta, yini, ydelta;
   int    ncols, nrows, col1, col2, row1, row2;
   int    ibox, nboxes, xyscan;  // 1 = scan in X; 2 = scan in Y
   int    i, j, k, n;
   long   offs;

   // Command line.
   if (argc < 2)
   {
      printf("\n  Usage:> %s <filename>\n\n", argv[0]);
      return(-1);
   }

   // Open and read the number of boxes and the scan type.
   if ((fp = fopen(argv[1], "rb")) == NULL)
   {
      printf("\n  ** Error: Unable to open the table '%s'\n\n", argv[1]);
      return(-1);
   }
   if (fread(&num, sizeof(float), 1, fp) != 1)
   {
      printf("\n  ** Error: Unable to read the table '%s'\n\n", argv[1]);
      return(-1);
   }
   nboxes = (int) num;
   if (fread(&num, sizeof(float), 1, fp) != 1)
   {
      printf("\n  ** Error: Unable to read the table '%s'\n\n", argv[1]);
      return(-1);
   }
   xyscan = (int) num;

   // Read the info of the boxes.
   info = malloc(nboxes * sizeof(*info));  // memory for 'nboxes' rows
   for (k = 0; k < nboxes; k++)
   {
      if (fread(info[k], sizeof(float), NINFO, fp) != NINFO)
      {
         printf("\n  ** Error: Unable to read the table '%s'\n\n", argv[1]);
         return(-1);
      }
   }

   printf("\n  Table make up with %d box/es", nboxes);
   if (nboxes == 1)
   {
      ibox = 1;
   }
   else
   {
      printf("\n  Which box do you want to print: ");
      if (scanf("%d", &ibox) != 1 || ibox < 1 || ibox > nboxes)
      {
         printf("\n  ** Error: Invalid box\n\n");
         return(-1);
      }
   }
   ncols  = (int) info[ibox-1][0];
   xini   =       info[ibox-1][1];
   xdelta =       info[ibox-1][2];
   nrows  = (int) info[ibox-1][3];
   yini   =       info[ibox-1][4];
   ydelta =       info[ibox-1][5];

   printf("\n  Box with %d columns (X) and %d rows (Y)", ncols, nrows);
   printf("\n  xini = %-10g  xfin = %-10g  xdelta = %-10g",
          xini, xini + (xdelta*(ncols-1)), xdelta);
   printf("\n  yini = %-10g  yfin = %-10g  ydelta = %-10g\n",
          yini, yini + (ydelta*(nrows-1)), ydelta);

   printf("\n  Please, introduce the range in columns and rows to print.");
   printf("\n  col-ini  col-fin  row-ini  row-fin ");
   printf("(max %d %d %d %d): ", 1, ncols, 1, nrows);
   n = scanf("%d %d %d %d", &col1, &col2, &row1, &row2);
   if ( n != 4 ||
        col1 < 1 || col1 > ncols || col2 < col1 || col2 > ncols ||
        row1 < 1 || row1 > nrows || row2 < row1 || row2 > nrows )
   {
      printf("\n  ** Error: Incorrect range (maximum range 1 %d 1 %d)\n\n",
             ncols, nrows);
      return(-1);
   }

   // Go to the beginning of the range.
   offs = 0;
   for (i = 0; i < (ibox-1); i++)
   {
      offs += (int) info[i][0] * (int) info[i][3];  // skip previous boxes
   }
   offs += ncols*(row1-1) + (col1-1);   // skip until the indicated beginning
   fseek(fp, offs * sizeof(float), SEEK_CUR);

   // Read the table values.
   printf("\n");
   for (i = row1; i <= row2; i++)
   {
      for (j = col1; j <= col2; j++)
      {
         fread(&num, sizeof(float), 1, fp);
         printf("  %10.3e", num);  // -1.123e+12 => 10 espacios
      }
      printf("\n");
      if (i < row2)
      {
         fseek(fp, (ncols - (col2-col1) - 1) * sizeof(float), SEEK_CUR);
      }
   }
   printf("\n");

   // Close the table and free the dynamic memory.
   fclose(fp);
   free(info);

   return(0);
}


//____________________________________________________________________________
