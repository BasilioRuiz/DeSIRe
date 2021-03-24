//
// File          _______ : demo.c
// Description   _______ : Demo using the memory leak detection library.
// Project       _______ : DeSIRe
// Creation date _______ : 02/02/21
// Author        _______ : epm@iac.es
//


#include <stdio.h>
#include <stdlib.h>

#include "cMemDbg.h"  // after standard headers

struct contact
{
  char *name;
  int   age;
};


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
   char  *string1, *string2;
   struct contact *agenda;

   printf("sizeof(struct contact) = %d\n\n", sizeof(struct contact));

   // Allocations.

   string1 = calloc(1, 10);
   string2 = malloc(20);
   string2 = realloc(string2, 30);

   agenda = calloc(2, sizeof(struct contact));
   agenda[0].name = malloc(40);
   agenda[0].age  = 40;
   agenda[1].name = realloc(NULL, 50);
   agenda[1].age  = 50;

   printf("\n*(1)*\tPrintTotalAllocatedMemory() - ANTES DE LIBERAR\n");
   PrintTotalAllocatedMemory();

   printf("\n*(2)*\tPrintMemoryReservedByCMemDbgLibrary()\n");
   PrintMemoryReservedByCMemDbgLibrary();

   // Liberations.

   printf("\n");
   //free(string1);
   //free(string2);
   //free(agenda[0].name);
   //free(agenda[1].name);
   free(agenda);

   printf("\n*(3)*\tPrintTotalAllocatedMemory() - DESPUES DE LIBERAR\n");
   PrintTotalAllocatedMemory();

   printf("\n*(4)*\tPrintMemoryReservedByCMemDbgLibrary()\n");
   PrintMemoryReservedByCMemDbgLibrary();

   // Memory leaks detection.

   printf("\n*(5)*\tPrintMemoryLeakInfo()\n");
   PrintMemoryLeakInfo();

   printf("\n---------------------------------");
   printf("\n*(6)*\tFreeMemoryLeakNoRealloc()\n");
   FreeMemoryLeakNoRealloc();
   printf("---------------------------------\n");

   printf("\n*(7)*\tPrintTotalAllocatedMemory() - AL TERMINAR\n");
   PrintTotalAllocatedMemory();

   printf("\n*(8)*\tPrintMemoryReservedByCMemDbgLibrary()\n");
   PrintMemoryReservedByCMemDbgLibrary();

   printf("\n*(9)*\tPrintMemoryLeakInfo()\n");
   PrintMemoryLeakInfo();
   printf("\n");

   return(EXIT_SUCCESS);
}


//____________________________________________________________________________
