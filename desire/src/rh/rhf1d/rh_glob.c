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

// 10/10/19 epm: Structure to pass information from DeSIRe to RH.
DesireLines desirelines;


//____________________________________________________________________________
