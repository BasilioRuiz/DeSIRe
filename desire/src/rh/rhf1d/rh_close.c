//
// File          _______ : rh_close.c
// Description   _______ : Close open writing RH files.
// Project       _______ : DeSIRe
// Creation date _______ : 09/09/19
// Author        _______ : epm@iac.es
//


#include <unistd.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "statistics.h"
#include "inputs.h"

void chrono_(unsigned long long *cpu_msec, unsigned long long *wall_msec);

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
extern ProgramStats stats;


//____________________________________________________________________________
//
//  Method: closefiles()
//
/** Close open writing RH files.
 */
//____________________________________________________________________________

void closefiles( void )
{
   Atom       *atom;
   AtomicLine *line;
   int         nact, i;
   unsigned long long cpu_msec1, cpu_msec2, wall_msec1, wall_msec2;

/*
   printf("\n");
   printf(" ... atmos.fd_background = %d\n", atmos.fd_background);
   printf(" ... atmos.fp_atmos      = %d\n", fileno(atmos.fp_atmos));
   printf(" ... commandline.logfile = %d\n", fileno(commandline.logfile));
   printf(" ... spectrum.fd_J       = %d\n", spectrum.fd_J);
   printf(" ... spectrum.fd_J20     = %d\n", spectrum.fd_J20);
   printf(" ... spectrum.fd_Imu     = %d\n", spectrum.fd_Imu);
   printf(" ... stats.fp_CPU        = %d\n", fileno(stats.fp_CPU));

   printf("     atmos.Natom         = %d\n", atmos.Natom);
   printf("     atmos.Nactiveatom   = %d\n", atmos.Nactiveatom);

   for (nact = 0; nact < atmos.Nactiveatom; nact++)
   {
     atom = atmos.activeatoms[nact];
     printf(" ... atom->fp_input[%d]   = %d\n", nact, fileno(atom->fp_input));
   }

   // To write 'line->fd_profile', LIMIT_MEMORY has to be TRUE in
   // "keyword.input". Otherwise, these descriptors have nonsense values.
   if (input.limit_memory)
   {
      for (nact = 0; nact < atmos.Nactiveatom; nact++)
      {
         atom = atmos.activeatoms[nact];
         printf("     atom[%d]->Nline = %d\n", nact, atom->Nline);
         for (i = 0; i < atom->Nline; i++)
         {
            line = &atom->line[i];
            printf(" ... fd_profile[%d]  = %d\n", i, line->fd_profile);
         }
      }
   }

   for (nact = 0; nact < atmos.Nactiveatom; nact++)
   {
     atom = atmos.activeatoms[nact];
     printf("     atom[%d]->Nline = %d\n", nact, atom->Nline);
     for (i = 0; i < atom->Nline; i++)
     {
       line = &atom->line[i];
       printf(" ... fp_GII[%d] = %d\n", i, line->fp_GII);
     }
   }
   printf("\n");
*/

   chrono_(&cpu_msec1, &wall_msec1);

   // Flushing files: About 40 ms per open file.
   // (Similar flushing individual files).

   //fflush(NULL); // flush every open output stream
   //sync();       // flush every open output file descriptor

   // Closing files: The same, about 40 ms per open file.
   // man close: A successful close does not guarantee that the data has been
   // successfully saved to disk, as the kernel defers writes. It is not
   // common for a file system to flush the buffers when the stream is closed.
   // If you need to be sure that the data is physically stored use fsync.

   //printf(" !!! FICHEROS QUE CIERRA closefiles():\n");
   if (atmos.fd_background > 0)                            // !! 1 see below
   {
      //printf(" !!! . atmos.fd_background (fd = %d)\n", atmos.fd_background);
      fsync(atmos.fd_background);
      close(atmos.fd_background);
      atmos.fd_background = 0;
   }
   if (spectrum.fd_J > 0)                                  // !! 4 see below
   {
      //printf(" !!! . spectrum.fd_J (fd = %d)\n", spectrum.fd_J);
      fsync(spectrum.fd_J);
      close(spectrum.fd_J);
      spectrum.fd_J = 0;
   }
   if (spectrum.fd_J20 > 0)                                // !! 4 see below
   {
      //printf(" !!! . spectrum.fd_J20 (fd = %d)\n", spectrum.fd_J20);
      fsync(spectrum.fd_J20);
      close(spectrum.fd_J20);
      spectrum.fd_J20 = 0;
   }
   if (spectrum.fd_Imu > 0)                                // !! 6 see below
   {
      //printf(" !!! . spectrum.fd_Imu (fd = %d)\n", spectrum.fd_Imu);
      fsync(spectrum.fd_Imu);
      close(spectrum.fd_Imu);
      spectrum.fd_Imu = 0;
   }

   for (nact = 0; nact < atmos.Nactiveatom; nact++)
   {
      atom = atmos.activeatoms[nact];
      for (i = 0;  i < atom->Nline; i++)
      {
         line = &atom->line[i];

         if (line->fd_profile > 0 && input.limit_memory)   // !! 9 see below
         {
            //printf(" !!! . line->fd_profile (fd = %d)\n", line->fd_profile);
            fsync(line->fd_profile);
            close(line->fd_profile);
            line->fd_profile = 0;
         }

         if (line->fp_GII != NULL)                         // !! 10 see below
         {
            //printf(" !!! . line->fp_GII (fd = %d)\n", fileno(line->fp_GII));
            fflush(line->fp_GII);
            fclose(line->fp_GII);
            line->fp_GII = NULL;
         }
      }
   }

   chrono_(&cpu_msec2, &wall_msec2);
   //printf(" !!! CONSUMO DE closefiles() = %llu ms\n",
   //       (wall_msec2 - wall_msec1));
}


//____________________________________________________________________________
//============================================================================
/*
                   GLOBAL DESCRIPTORS
                   ------------------

!! 1) atmos.h :: atmos.fd_background (global)
      - background.c     (open)
      - readj.c          (pread, pwrite)   !!!! force flush
      - rhf1d/solveray.c (open)

   2) atmos.h :: atmos->fp_atmos (local, argument pointing to a global)
      . backgrcontr.c call MULTIatmos(&atmos, &geometry), extern atmos
      . conversion.c  call MULTIatmos(&atmos, &geometry), extern atmos
      . rhf1d.c       call MULTIatmos(&atmos, &geometry), extern atmos
      . solveray.c    call MULTIatmos(&atmos, &geometry), extern atmos
      - rhf1d/multiatmos.c (fopen, getLine)   => OK, only readings

   3) inputs.h :: commandline.logfile (global)
      - error.c         (fprintf)
      - maxchange.c     (fprintf)
      - options.c       (fopen)   => OK, open once
      - rhf1d/solve1d.c (commandline.logfile = stderr)

!! 4) spectrum.h :: spectrum.fd_J|fd_J20 (global)
      - initial_xdr.c (open, close)
      - readj.c       (pread, pwrite)   !!!! just in case force flush

   5) spectrum.h :: spectrum->fd_J|fd_J20 (local: arg pointing to a global)
      . rhf1d.c call writeSpectrum(&spectrum), extern spectrum
      - writespect_xdr.c (open, pwrite, close)   => OK, close in situ

!! 6) spectrum.h :: spectrum.fd_Imu (global)
      - initial_xdr.c (open)
      - readj.c       (pread, pwrite)   !!!! force flush

   7) statistics.h :: stats.fp_CPU (global)
      - getcpu.c (fopen, fprintf)   => OK, open once

   8) atom.h :: atom->fp_input (local, argument pointing to a global)
      - readatom.c (fopen, getLine, fclose inside an 'if')
      - ltepops.c  (getLine through CollisionRate())   => OK, only readings

!! 9) atom.h :: line->fd_profile (local: argument pointing to a global)
      - readj.c   (pread, pwrite)   !!!! force flush
      - profile.c (open)

!!10) atom.h :: line->fp_GII (local: argument pointing to a global)
      - readatom.c (fclose)
      - scatter.c  (fopen, fread, fwrite)   => I am not sure

Be aware of local descriptors, although I think they are properly closed.

*/
//============================================================================
