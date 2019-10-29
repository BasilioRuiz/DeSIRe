//
// File          _______ : rh_chrono.c
// Description   _______ : Chronometer functions.
// Project       _______ : DeSIRe
// Creation date _______ : 06/05/19
// Author        _______ : epm@iac.es
//


#include <stdio.h>    // printf()
#include <stdlib.h>   // EXIT_SUCCESS
#include <time.h>     // clock_gettime(), clock()
#include <unistd.h>   // usleep()

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

void cpuclock_( unsigned long long *zero );
void wallclock_( unsigned long long *zero );


//____________________________________________________________________________
//
//  Method: chrono_()
//
/** Measure times.
 *  @param  cpu_msec CPU time in milliseconds since the program was launched.
 *  @param  wall_msec wall time in milliseconds since the Epoch (01/01/1970).
 */
//____________________________________________________________________________

void chrono_( unsigned long long *cpu_msec, unsigned long long *wall_msec )
{
   cpuclock_(cpu_msec);
   wallclock_(wall_msec);
}


//____________________________________________________________________________
//
//  Method: cpuclock_()
//
/** Measure CPU time. CPU time (or process time) is the amount of time for
 *  which a central processing unit (CPU) was used for processing instructions
 *  of a computer program or operating system.
 *  @param  cpu_msec CPU time in milliseconds since the program was launched.
 */
//____________________________________________________________________________

void cpuclock_( unsigned long long *cpu_msec )
{
   *cpu_msec = clock() * 1000L / CLOCKS_PER_SEC;
}


//____________________________________________________________________________
//
//  Method: wallclock()
//
/** Measure elapsed real time. Elapsed real time (or simply real time, or
 *  wall-clock time) is the time taken from the start of a computer program
 *  until the end as measured by an ordinary clock. Elapsed real time includes
 *  I/O time, any multitasking delays, and all other types of waits incurred
 *  by the program.
 *  @param  wall_msec wall time in milliseconds since the Epoch (01/01/1970).
 */
//____________________________________________________________________________

void wallclock_( unsigned long long *wall_msec )
{
   struct timespec ts;

#ifdef __MACH__   // OS X does not have clock_gettime, use clock_get_time

   clock_serv_t cclock;
   mach_timespec_t mts;
   host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
   clock_get_time(cclock, &mts);
   mach_port_deallocate(mach_task_self(), cclock);
   ts.tv_sec = mts.tv_sec;
   ts.tv_nsec = mts.tv_nsec;

#else

   // All implementations support the system-wide realtime clock, which is
   // identified by CLOCK_REALTIME. Its time represents seconds and
   // nanoseconds since the Epoch. When its time is changed, timers for a
   // relative interval are unaffected, but timers for an absolute point in
   // time are affected.
   clock_gettime(CLOCK_REALTIME, &ts);

#endif

   // Don't return nanoseconds, that is, don't multiply the seconds by 1e9
   // because it can overflow even the long long 64 bits variable.
   // A millisecond precision is enough.
   *wall_msec = (ts.tv_sec * 1000L) + (ts.tv_nsec / 1000000L);
}


//____________________________________________________________________________
//
//  Method: main()
//
/** Main method for debugging.
 *  @param  argc number of arguments passed to the program.
 *  @param  argv array of strings with the arguments passed.
 *  @return exit code (0 if successful).
 */
//____________________________________________________________________________

int xmain( int argc, char *argv[] )
{
   int i;
   unsigned long long cpu_msec1, cpu_msec2, wall_msec1, wall_msec2;

   chrono_(&cpu_msec1, &wall_msec1);

   for (i = 0; i < 100000000; i++) { ; }
   usleep(2000000);   // 2 sec

   chrono_(&cpu_msec2, &wall_msec2);

   printf("cpuclock()  = %llu ms\n", (cpu_msec2 - cpu_msec1));
   printf("wallclock() = %llu ms\n", (wall_msec2 - wall_msec1));

   return(EXIT_SUCCESS);
}


//____________________________________________________________________________
//============================================================================
/*
                   BEST TIMING METHOD IN C
                   -----------------------

clock() measure the CPU time used by your process, not the wall-clock time.
When you have multiple threads running simultaneously, you can obviously burn
through CPU time much faster.

If you want to know the wall-clock execution time, you need to use an
appropriate function. The only one in ANSI C is time(), which typically only
has 1 second resolution.

However, if you're using POSIX, that means you can use clock_gettime(),
defined in "time.h". The CLOCK_MONOTONIC clock in particular is the best to
use for this.

If your OS doesn't provide CLOCK_MONOTONIC (which you can check at runtime
with sysconf(_SC_MONOTONIC_CLOCK)), then you can use CLOCK_REALTIME as a
fallback, but note that the latter has the disadvantage that it will generate
incorrect results if the system time is changed while your process is running.

CLOCK_REALTIME
   System-wide clock that measures real (i.e., wall-clock) time.  Setting
   this clock requires appropriate privileges. This clock is affected by
   discontinuous jumps in the system time (e.g. if the system administrator
   manually changes the clock) and by the incremental adjustments performed
   by adjtime(3) and NTP.

CLOCK_MONOTONIC
   Clock that cannot be set and represents monotonic time since some
   unspecified starting point.  This clock is not affected by discontinuous
   jumps in the system time (e.g., if the system administrator manually
   changes the clock), but is affected by the incremental adjustments
   performed by adjtime(3) and NTP.

CLOCK_BOOTTIME (since Linux 2.6.39; Linux-specific)
   Identical  to CLOCK_MONOTONIC, except it also includes any time that the
   system is suspended.  This allows applications to get a suspend-aware
   monotonic clock without having to deal with the complications of
CLOCK_REALTIME, which may have discontinuities if the time is changed
   using settimeofday(2).

CLOCK_PROCESS_CPUTIME_ID (since Linux 2.6.12)
   Per-process CPU-time clock (measures CPU time consumed by all threads
   in the process).

CLOCK_THREAD_CPUTIME_ID (since Linux 2.6.12)
   Thread-specific CPU-time clock.

*/
//============================================================================
/*
                   MEDIDA DEL TIEMPO TOTAL DE EJECUCION
                   ------------------------------------

El comando time de Linux permite obtener el tiempo total de ejecucion de un
programa. Su sintaxis es muy sencilla: 'time' seguido del comando cuya
ejecucion queremos medir (con los parametros correspondientes).

El resultado de time se escribe en la salida de error estandar (en linea de
comandos "2>"). Por ejemplo, si ejecutamos desde la linea de comandos:

> (time ls) 2> salida.txt

Podemos obtener en "salida.txt" un resultado del tipo:

    real 0m0.037s
    user 0m0.004s
    sys 0m0.008s

El valor "real" indica el tiempo total transcurrido en ejecutar el comando.
Por ejemplo, si hay otros procesos en el sistema, se contara tambien el tiempo
de los mismos. El valor "user" se refiere al tiempo de CPU del proceso en
cuestion, por lo tanto, se excluye el tiempo de otros procesos o de los
retardos del disco. El valor "sys" es el tiempo de CPU en las llamadas al
sistema del proceso. Idealmente, si no hubieran otros procesos y la lectura
de disco fuera inmediata, tendriamos que: user + sys = real.

Los parentesis no se pueden suprimir en el comando indicado, ya que en otro
caso se tomaria como:

    >> time (ls 2> salida.txt)

Lo cual tiene un resultado bien distinto.

*/
//============================================================================
