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
#include "error.h"

void chrono_( unsigned long long *cpu_msec, unsigned long long *wall_msec );

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
extern ProgramStats stats;
extern char messageStr[];


//____________________________________________________________________________
//
//  Method: closefiles()
//
/** Close open writing RH files.
 */
//____________________________________________________________________________

void closefiles( void )
{
   const char routineName[] = "closefiles";

   Atom       *atom;
   AtomicLine *line;
   int         nact, i;
   unsigned long long cpu_msec1, cpu_msec2, wall_msec1, wall_msec2;

/*
   printf("\n");
   printf(" ... atmos.fd_background = %d\n", atmos.fd_background);
   if (atmos.fp_atmos)
   printf(" ... atmos.fp_atmos      = %d\n", fileno(atmos.fp_atmos));
   if (commandline.logfile)
   printf(" ... commandline.logfile = %d\n", fileno(commandline.logfile));
   printf(" ... spectrum.fd_J       = %d\n", spectrum.fd_J);
   printf(" ... spectrum.fd_J20     = %d\n", spectrum.fd_J20);
   printf(" ... spectrum.fd_Imu     = %d\n", spectrum.fd_Imu);
   if (stats.fp_CPU)
   printf(" ... stats.fp_CPU        = %d\n", fileno(stats.fp_CPU));

   printf("                           atmos.Natom       : %d\n",
          atmos.Natom);
   printf("                           atmos.Nactiveatom : %d\n",
          atmos.Nactiveatom);

   for (nact = 0; nact < atmos.Nactiveatom; nact++)
   {
     atom = atmos.activeatoms[nact];
     if (atom->fp_input)
     printf(" ... atom[%d]->fp_input = %d\n", nact, fileno(atom->fp_input));
   }

   // To write 'line->fd_profile', LIMIT_MEMORY has to be TRUE in
   // "keyword.input". Otherwise, these descriptors have nonsense values.
   if (input.limit_memory)
   {
      for (nact = 0; nact < atmos.Nactiveatom; nact++)
      {
         atom = atmos.activeatoms[nact];
         printf("                           atom[%d]->Nline : %d\n",
                nact, atom->Nline);
         for (i = 0; i < atom->Nline; i++)
         {
            line = &atom->line[i];
            printf(" ... line[%d]->fd_profile = %d\n", i, line->fd_profile);
         }
      }
   }

   for (nact = 0; nact < atmos.Nactiveatom; nact++)
   {
     atom = atmos.activeatoms[nact];
     printf("                           atom[%d]->Nline : %d\n",
            nact, atom->Nline);
     for (i = 0; i < atom->Nline; i++)
     {
       line = &atom->line[i];
       if (line->fp_GII)
       printf(" ... line[%d]->fp_GII = %d\n", i, fileno(line->fp_GII));
     }
   }
   printf("\n");
*/

   // chrono_(&cpu_msec1, &wall_msec1);

   // Flushing files: About 40 ms per open file.
   // (Similar flushing individual files).

   //fflush(NULL); // flush every open output stream
   //sync();       // flush every open output file descriptor

   // Closing files: The same, about 40 ms per open file.
   // man close: A successful close does not guarantee that the data has been
   // successfully saved to disk, as the kernel defers writes. It is not
   // common for a file system to flush the buffers when the stream is closed.
   // If you need to be sure that the data is physically stored use fsync.

   if (atmos.fd_background > 0)                            // !! 1 see below
   {
      sprintf(messageStr, "Closing the open file atmos.fd_background "
                          "(fd = %d) before going on", atmos.fd_background);
      Error(WARNING, routineName, messageStr);
      fsync(atmos.fd_background);
      close(atmos.fd_background);
      atmos.fd_background = 0;
   }
   if (spectrum.fd_J > 0)                                  // !! 4 see below
   {
      sprintf(messageStr, "Closing the open file spectrum.fd_J "
                          "(fd = %d) before going on", spectrum.fd_J);
      Error(WARNING, routineName, messageStr);
      fsync(spectrum.fd_J);
      close(spectrum.fd_J);
      spectrum.fd_J = 0;
   }
   if (spectrum.fd_J20 > 0)                                // !! 4 see below
   {
      sprintf(messageStr, "Closing the open file spectrum.fd_J20 "
                          "(fd = %d) before going on", spectrum.fd_J20);
      Error(WARNING, routineName, messageStr);
      fsync(spectrum.fd_J20);
      close(spectrum.fd_J20);
      spectrum.fd_J20 = 0;
   }
   if (spectrum.fd_Imu > 0)                                // !! 6 see below
   {
      sprintf(messageStr, "Closing the open file spectrum.fd_Imu "
                          "(fd = %d) before going on", spectrum.fd_Imu);
      Error(WARNING, routineName, messageStr);
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
            sprintf(messageStr,
                    "Closing the open file line[%d]->fd_profile "
                    "(fd = %d) before going on", i, line->fd_profile);
            Error(WARNING, routineName, messageStr);
            fsync(line->fd_profile);
            close(line->fd_profile);
            line->fd_profile = 0;
         }

         if (line->fp_GII != NULL)                         // !! 10 see below
         {
            sprintf(messageStr,
                    "Closing the open file line[%d]->fp_GII "
                    "(fd = %d) before going on", i, fileno(line->fp_GII));
            Error(WARNING, routineName, messageStr);
            fflush(line->fp_GII);
            fclose(line->fp_GII);
            line->fp_GII = NULL;
         }
      }
   }

   // chrono_(&cpu_msec2, &wall_msec2);
   // printf(" !!! CONSUMO DE closefiles() = %llu ms\n",
   //        (wall_msec2 - wall_msec1));
}


//____________________________________________________________________________
//============================================================================
/*
                   GLOBAL DESCRIPTORS
                   ------------------
Reopen an open file rewind the file. Be aware of that if you are reading
and writing the file since the writing might not be dumped into disk before
the next reading.

!! 1) atmos.fd_background (global int in atmos.h)
      - background.c     (open)
      - readj.c          (pread, pwrite)   !!!! force flush
      - rhf1d/solveray.c (open)
   => 04/04/20 epm: Save in memory rather than disk.

   2) atmos->fp_atmos (local, argument pointing to a global FILE in atmos.h)
      . backgrcontr.c calls MULTIatmos(&atmos, &geometry), extern atmos
      . conversion.c  calls MULTIatmos(&atmos, &geometry), extern atmos
      . rhf1d.c       calls MULTIatmos(&atmos, &geometry), extern atmos
      . solveray.c    calls MULTIatmos(&atmos, &geometry), extern atmos
      - rhf1d/multiatmos.c (fopen, getLine)   => OK, only readings

   3) commandline.logfile (global FILE in inputs.h)
      - error.c         (fprintf)
      - maxchange.c     (fprintf)
      - options.c       (fopen)   => OK, open once
      - rhf1d/solve1d.c (commandline.logfile = stderr)

!! 4) spectrum.fd_J|fd_J20 (global int in spectrum.h)
      - initial_xdr.c (open, close)
      - readj.c       (pread, pwrite)   !!!! just in case force flush
   => 04/04/20 epm: Save in memory rather than disk. And fd_J20 never used.

   5)  spectrum->fd_J|fd_J20 (local: arg pointing to a global in spectrum.h)
      . rhf1d.c calls writeSpectrum(&spectrum), extern spectrum
      - writespect_xdr.c (open, pwrite, close)   => OK, close in situ

!! 6) spectrum.fd_Imu (global int in spectrum.h)
      - initial_xdr.c (open)
      - readj.c       (pread, pwrite)   !!!! force flush
   => 04/04/20 epm: Never used.

   7) stats.fp_CPU (global FILE in statistics.h)
      - getcpu.c (fopen, fprintf)   => OK, open once

   8) atom->fp_input (local, argument pointing to a global FILE in atom.h)
      - readatom.c (fopen, getLine, fclose inside an 'if')
      - ltepops.c  (getLine through CollisionRate())   => OK, only readings

!! 9) line->fd_profile (local: argument pointing to a global int in atom.h)
      - readj.c   (pread, pwrite)   !!!! force flush
      - profile.c (open)
   => 04/04/20 epm: LIMIT_MEMORY = FALSE.

!!10) line->fp_GII (local: argument pointing to a global FILE in atom.h)
      - readatom.c (fclose)
      - scatter.c  (fopen, fread, fwrite)   => I am not sure
   => 04/04/20 epm: Never used.

Be aware of local descriptors, although I think they are properly closed.

*/
//============================================================================
