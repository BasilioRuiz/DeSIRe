/* ------- file: -------------------------- options.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Mar 25 14:46:56 2009 --

       --------------------------                      ----------RH-- */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "inputs.h"
#include "error.h"
#include "desire.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern CommandLine commandline;
extern char messageStr[];

// 04/04/20 epm: Structure to load some command line flags from SIR.
extern SIRflags sirflags;


/* ------- begin -------------------------- setOptions.c ------------ */

void setOptions(int argc, char *argv[])
{
  const  char routineName[] = "setOptions";
  static char logfileName[MAX_LINE_SIZE], wavetable[MAX_LINE_SIZE];

  int Noption;

  Option theOptions[] = {
    {"help", 1, FALSE, "", NULL, NULL, "Prints this message"},
    {"input", 1, TRUE, "keyword.input",
       commandline.keyword_input,
       setcharValue, "File name for input keywords"},
    {"logfile", 1, TRUE, "",
       logfileName,
       setcharValue, "File name log file"},
    {"quiet", 1, FALSE, "FALSE", &commandline.quiet, setboolValue,
       "Turns off warning messages"},
    {"showkeywords", 1, FALSE, "FALSE", &commandline.showkeywords,
       setboolValue,
       "Show keyword values with current keyword input file"}
  };
  Noption = sizeof(theOptions) / sizeof(Option);

  parse(argc, argv, Noption, theOptions);

  if (strlen(logfileName) > 0)
  {
    // 13/06/19 epm: After converting RH executables to functions,
    // the file remains open (there is not end of application to close it)
    // and, therefore, we should only open it once.
    if (commandline.logfile == NULL)
    {
      if ((commandline.logfile = fopen(logfileName, "w")) == NULL)
      {
        sprintf(messageStr, "Unable to open log file %s", logfileName);
        Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      setvbuf(commandline.logfile, NULL, _IOLBF, BUFSIZ);
    }
  }
  else
  {
    // 04/04/20 epm: Change to separate outputs on screen from logfile.
    // commandline.logfile = stderr;
  }
  commandline.wavetable = (strlen(wavetable) > 0) ? wavetable : NULL;

  // 04/04/20 epm: Overwrite with SIR command line verbose flag.
  if (sirflags.flagv == 1) commandline.quiet = FALSE;
  else commandline.quiet = TRUE;
}
/* ------- end ---------------------------- setOptions.c ------------ */
