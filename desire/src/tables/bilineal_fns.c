//
// File          _______ : bilineal_fns.c
// Description   _______ : Bilineal interpolation functions.
// Project       _______ : DeSIRe
// Creation date _______ : 12/12/22
// Author        _______ : epm@iac.es
//


#include <stdio.h>
#include <stdlib.h>

#define MAX_TABLES  100  // maximum total number of tables
#define NINFO         6  // number of info values of each box

static FILE  *fp[MAX_TABLES];
static float *header[MAX_TABLES];
static int    ntables = 0;



//____________________________________________________________________________
//
//  Method: opentable_()
//
/** Open an interpolation table and read the header.
 *  @param  ncodes size of 'codes'.
 *  @param  codes array with the ASCII codes for the table filename.
 *  @param  vheader the header as a vector (output).
 *  @return table descriptor index (<0 if error).
 */
//____________________________________________________________________________

int opentable_( int *ncodes, int *codes, float *vheader )
{
   char  filename[*ncodes + 1];
   float fboxes;
   int   nboxes;
   int   nelem;
   int   i;

   // Table already opened.
   if (fp[ntables] != NULL)
   {
      return(ntables);
   }

   // Get the table filename from the ASCII codes.
   for (i = 0; i < *ncodes; i++)
   {
      filename[i] = (char) codes[i];
   }
   filename[i] = '\0';

   // Open the table.
   if ((fp[ntables] = fopen(filename, "rb")) == NULL)
   {
      return(-1);
   }

   // Read the number of boxes.
   if (fread(&fboxes, sizeof(float), 1, fp[ntables]) != 1)
   {
      return(-2);
   }
   nboxes = (int) fboxes;

   // Alloc memory to save the header.
   nelem = 2 + (nboxes * NINFO);
   if ((header[ntables] = malloc(nelem * sizeof(float))) == NULL)
   {
      return(-3);
   }
   header[ntables][0] = fboxes;

   // Read the remaining header.
   if (fread(&header[ntables][1], sizeof(float), nelem-1, fp[ntables])
       != nelem-1)
   {
      return(-2);
   }

   // Return the header as a vector.
   for (i = 0; i < nelem; i++)
   {
      vheader[i] = header[ntables][i];
   }
   ntables++;

   return(ntables-1);
}


//____________________________________________________________________________
//
//  Method: closetable_()
//
/** Close the table.
 *  @param  itable table descriptor index.
 *  @return error code (0 if successful).
 */
//____________________________________________________________________________

int closetable_( int *itable )
{
   int ier;
   
   if (fp[*itable] == NULL)
   {
      return(0);
   }
   ier = fclose(fp[*itable]);
   fp[*itable] = NULL;
   free(header[*itable]);

   return(ier);
}


//____________________________________________________________________________
//
//  Method: bilineal_()
//
/** Make a bilineal interpolation using a table on disk.
 *  @param  itable table descriptor index.
 *  @param  x x coordinate for the interpolation.
 *  @param  y y coordinate for the interpolation.
 *  @param  f resulting interpolated value.
 *  @return error code (0 if successful).
 */
//____________________________________________________________________________


int bilineal_( int *itable, float *x, float *y, float *f )
{
   float  q[4];
   float  Q11, Q21, Q12, Q22, R1, R2, P;
   float  x1, x2, y1, y2;
   float  xini = 0, yini = 0, xdelta = 0, ydelta = 0, xfin = 0, yfin = 0;
   int    ncols = 0, nrows = 0, nboxes = 0, nelem = 0, xyscan = 0;
   int    i, ibox, icol, irow, offs;

   if (fp[*itable] == NULL || header[*itable] == NULL)
   {
      return(-1);
   }
   nboxes = (int) header[*itable][0];
   xyscan = (int) header[*itable][1];  // 1 = scan in X; 2 = scan in Y
   nelem  = 2 + (nboxes * NINFO);      // number of elements of the header

   // Find out the box with the target coordinate.
   for (ibox = 0; ibox < nboxes; ibox++)
   {
      ncols  = (int) header[*itable][ibox * NINFO + 2];
      xini   =       header[*itable][ibox * NINFO + 3];
      xdelta =       header[*itable][ibox * NINFO + 4];
      nrows  = (int) header[*itable][ibox * NINFO + 5];
      yini   =       header[*itable][ibox * NINFO + 6];
      ydelta =       header[*itable][ibox * NINFO + 7];
      xfin   = xini + (xdelta * (ncols - 1));
      yfin   = yini + (ydelta * (nrows - 1));

      if (xyscan == 1 && *x < xfin)
      {
         if (*y >= yini && *y < yfin) break;  // point in the right box 
         else continue;                       // point overlapped in other box
      }
      if (xyscan == 2 && *y < yfin)
      {
         if (*x >= xini && *x < xfin) break;  // point in the right box
         else continue;                       // point overlapped in other box
      }
   }
   if (ibox == nboxes)
   {
      if (xyscan == 1) return(-2);
      if (xyscan != 1) return(-3);
   }

   // Check the point coordinates within the box.
   if (*x < xini || *x >= xfin) return(-2);
   if (*y < yini || *y >= yfin) return(-3);

   /**************************************************************************

      Bilineal interpolation in f(P) = f(x,y).

      Q11=(x1,y1), Q12=(x1,y2), Q21=(x2,y1), Q22=(x2,y2)


             x1             x                           x2
       _____________________________________________________
      |      :              :                            :
      |      :              :                            :
   y1 |.....Q11.............R1..........................Q21
      |      :              :                            :
      |      :              :                            :
      |      :              :                            :
      |      :              :                            :
    y |......:..............P............................:
      |      :              :                            :
      |      :              :                            :
      |      :              :                            :
   y2 |.....Q12.............R2..........................Q22
      |


              x2 - x              x - x1
      f(R1) = ------- * f(Q11) + ------- * f(Q21)    where R1 = (x,y1)
              x2 - x1            x2 - x1


              x2 - x              x - x1
      f(R2) = ------- * f(Q12) + ------- * f(Q22)    where R2 = (x,y2)
              x2 - x1            x2 - x1


              y2 - y              y - y1
      f(P)  = ------- * f(R1)  + ------- * f(R2)
              y2 - y1            y2 - y1

   **************************************************************************/

   // Find out the values in the table surrounding the point.
   offs = nelem;
   for (i = 0; i < ibox; i++)
   {
      offs += (int) (header[*itable][i*NINFO+2] * header[*itable][i*NINFO+5]);
   }
   irow = (int) ((*y - yini) / ydelta);   // row before
   icol = (int) ((*x - xini) / xdelta);   // column before
   offs += (ncols * irow) + icol;
   fseek(fp[*itable], offs * sizeof(float), SEEK_SET);
   if (fread(&q[0], sizeof(float), 2, fp[*itable]) != 2)
   {
      return(-4);
   }
   fseek(fp[*itable], (ncols - 2) * sizeof(float), SEEK_CUR);
   if (fread(&q[2], sizeof(float), 2, fp[*itable]) != 2)
   {
      return(-4);
   }
   Q11 = q[0];
   Q21 = q[1];
   Q12 = q[2];
   Q22 = q[3];

   // Find out the values in the axes surrounding the point.
   x1 = xini + (xdelta * icol);
   x2 = x1 + xdelta;
   y1 = yini + (ydelta * irow);
   y2 = y1 + ydelta;

   // Interpolate in x direction.
   R1 = ((x2 - *x) / (x2 - x1)) * Q11 + ((*x - x1) / (x2 - x1)) * Q21;
   R2 = ((x2 - *x) / (x2 - x1)) * Q12 + ((*x - x1) / (x2 - x1)) * Q22;
   // Interpolate in y direction.
   P  = ((y2 - *y) / (y2 - y1)) * R1  + ((*y - y1) / (y2 - y1)) * R2;
   *f = P;

   return(0);
}


//____________________________________________________________________________
