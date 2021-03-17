/* ------- file: -------------------------- Error.c -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Dec 23 14:47:16 2008 --

       --------------------------                      ----------RH-- */

/* --- Print warning, or print error string and exit when error level
       is above preset treshold

       Modifications:

       - 04/04/20 epm:
         Change to separate outputs on screen from logfile.

       --                                              -------------- */

#include <stdlib.h>
#include <errno.h>

#include "rh.h"
#include "error.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern CommandLine commandline;


/* ------- begin -------------------------- Error.c ----------------- */

void Error(enum errorlevel level, const char *routineName,
           const char *messageStr)
{
  char str[MAX_MESSAGE_LENGTH];
  enum errorlevel defaultLevel = ERROR_LEVEL_1;

  switch (level)
  {
    case MESSAGE:
    {
      sprintf(str, "%s", (messageStr) ? messageStr : "");

      if (!commandline.quiet) fprintf(stdout, "%s", str);
      if (commandline.logfile) fprintf(commandline.logfile, "%s", str);
    }
    return;

    case WARNING:
    {
      sprintf(str, "-WARNING in routine %s\n %s\n\n",
              routineName, (messageStr) ? messageStr : " ( Undocumented )");

      if (!commandline.quiet) fprintf(stdout, "%s", str);
      if (commandline.logfile) fprintf(commandline.logfile, "%s", str);
    }
    return;

    default:
    {
      if (level < defaultLevel)
      {
        sprintf(str, "\a-ERROR in routine %s\n %s\n %s\n\n",
                routineName, (messageStr) ? messageStr : " ( Undocumented )",
                "Trying to continue.....");

        if (!commandline.quiet) fprintf(stdout, "%s", str);
        if (commandline.logfile) fprintf(commandline.logfile, "%s", str);

        return;
      }
      else
      {
        sprintf(str, "\a\n-TERMINATING_ERROR in routine %s\n %s\n %s\n\n",
                routineName, (messageStr) ? messageStr : " ( Undocumented )",
                "Exiting.....");

        fprintf(stdout, "%s", str);
        if (commandline.logfile) fprintf(commandline.logfile, "%s", str);

        exit(level);
      }
    }
  }
}
/* ------- end ---------------------------- Error.c ----------------- */
