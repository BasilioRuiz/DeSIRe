//
// File          _______ : rh_glob.c
// Description   _______ : Global variable declaration.
// Project       _______ : DeSIRe
// Creation date _______ : 23/05/19
// Author        _______ : epm@iac.es
//
// Remember that marking a C variable extern declares the variable without
// defining it. That is, no memory is allocated for it at that point in the
// program. Something elsewhere has to define the variable. And here they are.
// We can leave to the linker to find them, that is, we only need to add
// this module to the object list to link.
//


#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "statistics.h"
#include "inputs.h"
#include "desire.h"

enum Topology topology = ONE_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_MESSAGE_LENGTH];

// 02/07/19 epm: Boolean variables to control new RH executions.
bool_t new_background;   // for background.c :: Background()
bool_t new_hminus_ff;    // for hydrogen.c :: Hminus_ff()
bool_t new_h2minus_ff;   // for hydrogen.c :: H2minus_ff()
bool_t new_h2plus_ff;    // for hydrogen.c :: H2plus_ff()
bool_t new_passive_bb;   // for metal.c :: passive_bb()

// 10/10/19 epm: Structure to pass Barklem data from SIR to RH.
SIRBarklem sirbarklem;
// 04/04/20 epm: Structure to pass some command line flags from SIR to RH.
SIRflags sirflags;
// 07/07/20 epm: Structure to pass Kurucz data from SIR to RH.
SIRKurucz sirkurucz;
// 08/08/20 epm: Structure to pass abundance values from SIR to RH.
SIRabun sirabun = {0,
                   NULL,
                   {"H","HE","LI","BE","B","C","N","O","F","NE","NA","MG",
                    "AL","SI","P","S","CL","AR","K","CA","SC","TI","V","CR",
                    "MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR",
                    "KR","RB","SR","Y","ZR","NB","MO","TC","RU","RH","PD",
                    "AG","CD","IN","SN","SB","TE","I","XE","CS","BA","LA",
                    "CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER",
                    "TM","YB","LU","HF","TA","W","RE","OS","IR","PT","AU",
                    "HG","TL","PB","BI","PO","AT","RN","FR","RA","AC","TH",
                    "PA","U","NP","PU","AM","CM","BK","CF","ES"}
                  };
// 10/10/20 epm: Structure to pass the wavetable from SIR to RH.
SIRwave sirwave;
// 10/10/20 epm: Structure to pass some keywords from SIR to RH.
SIRkeywords sirkeywords;
// 11/11/20 epm: Structure to pass atmospheric models from SIR to RH.
SIRatmos siratmos;


//____________________________________________________________________________
