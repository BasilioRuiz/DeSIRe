
/*****************************************************************************
  This file is part of cMemDbg - Easy to use C memory leak detection library

        Copyright (C) 2009 Ezequiel Gastón Miravalles

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/


/*****************************************************************************
        Software: cMemDbg
        File: cMemDbg.c
        Version: 1.0
        Last modification: 8/4/2009
        Author: Ezequiel Gastón Miravalles
        Website: http://www.neoegm.com/software/cmemdbg/
        License: GNU GPL v3 (attached, read above)

        Notes: Read cMemDbg.h notes.

        Modification: 02/02/21
        Author: Esperanza Paez
        Email: epm@iac.es
        Change: Use of linked list instead of static array.
*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>

//Configuration Constants ***************
#define  PRINT_OPERATIONS  0
#define  PRINT_ERRORS      0
#define  MALLOC            0
#define  CALLOC            1
#define  REALLOC           2
//***************************************

struct DataAlloc
{
   void *pMem;
   int   nBytes;
   char *pszFile;
   int   nLine;
   int   type;
};
struct Node
{
   struct DataAlloc data;
   struct Node     *next;
};

int nAllocMem = 0;
int nNodes    = 0;
int nAllocs   = 0;
int nFrees    = 0;
struct Node *headList = NULL;


//============================================================================

// Add a new node to the linked list.

int AddAllocDataBytes(void *pPointer, size_t size, char *pszFile, int nLine,
                      int type)
{
   struct DataAlloc x;
   struct Node *node;

   if (!pPointer)
   {
      if (PRINT_ERRORS)
      {
         fprintf(stderr,
                 ">C_ERROR\tC function unable to allocate %ld bytes [%s:%d]\n",
                 size, pszFile, nLine);
      }
      return(0);
   }

   if ((node = malloc(sizeof(struct Node))) == NULL)
   {
      if (PRINT_ERRORS)
      {
         fprintf(stderr,
                 ">INTERNAL_ERROR\tError allocating memory for the node\n");
      }
      return(0);
   }

   x.pMem = pPointer;
   x.nBytes = size;
   x.pszFile = pszFile;
   x.nLine = nLine;
   x.type  = type;

   node->data = x;
   node->next = headList;   // add the node at the beggining
   headList = node;         // to use LIFO when removing nodes

   nAllocMem += size;
   nNodes++;
   nAllocs++;

   return(size);
}

//----------------------------------------------------------------------------

// Remove the node of the linked list matching a field.

int SubtractAllocDataBytes(void *pPointer, char **ppszFile, int *pnLine)
{
   struct DataAlloc x;
   struct Node *aux, *prev;

   if (!pPointer)
   {
      return(0);
   }

   for (prev=NULL, aux=headList;  aux != NULL;  prev=aux, aux=aux->next)
   {
      x = aux->data;

      if (x.pMem == pPointer)   // found!
      {
         if (prev == NULL)   // first node
         {
            headList = aux->next;
         }
         else   // following nodes
         {
            prev->next = aux->next;   // skip this node
         }
         free(aux);   // free the node

         if (ppszFile) *ppszFile = x.pszFile;
         if (pnLine) *pnLine = x.nLine;

         nAllocMem -= x.nBytes;
         nNodes--;
         nFrees++;

         return(x.nBytes);
      }
   }

   return(0);
}

//============================================================================

void PrintTotalAllocatedMemory()
{
   fprintf(stdout,
           ">INFO\tTotal allocated memory at this moment: %d bytes\n"
           ">INFO\t#allocs = %d   #frees = %d   #remainder = %d\n",
           nAllocMem, nAllocs, nFrees, nNodes);
}

//----------------------------------------------------------------------------

void PrintMemoryReservedByCMemDbgLibrary()
{
   fprintf(stdout,
           ">INFO\tMemory in use by the linked list: %ld bytes [%ld x %d]\n",
           sizeof(struct Node) * nNodes, sizeof(struct Node), nNodes);
}

//----------------------------------------------------------------------------

void PrintMemoryLeakBlocks()
{
   char   ctype[3] = {'M', 'C', 'R'};
   int    i = 1;
   struct DataAlloc x;
   struct Node *node = headList;

   while (node != NULL)
   {
      x = node->data;
      fprintf(stdout,
              ">NODE\t#%03d (%c)\t%p\t%d\t[%s:%d]\n",
              i, ctype[x.type], x.pMem, x.nBytes, x.pszFile, x.nLine);
      i++;
      node = node->next;
   }
}

//----------------------------------------------------------------------------

void PrintMemoryLeakInfo()
{
   if (nAllocMem == 0)
   {
      fprintf(stdout, ">INFO\tNo memory leaks detected\n");
   }
   else if (nAllocMem < 0)
   {
      fprintf(stdout,
              ">INFO\tPROBLEM: There was more freed memory than allocated"
              " (%d bytes) [This shouldn't happen, memory corruption?]\n",
              nAllocMem);
   }
   else
   {
      fprintf(stdout,
              ">INFO\tPROBLEM: Memory leak found (%d + %ld bytes)\n",
              nAllocMem, sizeof(struct Node) * nNodes);
   }

   PrintMemoryLeakBlocks();
}

//============================================================================

void FreeMemoryLeak()
{
   struct DataAlloc x;
   struct Node *node = headList;

   // To be consistent with bound allocations (like pointers as fields in a
   // struct), free the nodes under LIFO method: last in, first out.
   while (node != NULL)
   {
      x = node->data;
      headList = node->next;

      free(x.pMem);   // free the memory pointed by the node
      free(node);     // free the node
      node = headList;

      nAllocMem -= x.nBytes;
      nNodes--;
      nFrees++;
   }
}

//----------------------------------------------------------------------------

void FreeMemoryLeakNoRealloc()
{
   char   ctype[3] = {'M', 'C', 'R'};
   int    i = 1;
   struct DataAlloc x;
   struct Node *aux = headList, *prev = NULL;

   while (aux != NULL)
   {
      x = aux->data;

      if (x.type == REALLOC)   // do not free memory reserved with realloc()
      {
         prev = aux;
         aux = aux->next;

         if (PRINT_OPERATIONS)
         {
            fprintf(stdout,
                    ">NODE\t#%03d (%c)\t%p\t%d\t[%s:%d] => not freed\n",
                    i, ctype[x.type], x.pMem, x.nBytes, x.pszFile, x.nLine);
         }
      }
      else   // free memory reserved with malloc() or calloc()
      {
         if (aux == headList)         // first node:
         {
            headList = aux->next;     // first node will be the next
            free(x.pMem);             // free the memory pointed by the node
            free(aux);                // free the node
            aux = headList;           // reassign the auxiliar pointer
         }
         else                         // following nodes:
         {
            prev->next = aux->next;   // skip this node
            free(x.pMem);             // free the memory pointed by the node
            free(aux);                // free the node
            aux = prev->next;         // reassign the auxiliar pointer
         }

         nAllocMem -= x.nBytes;
         nNodes--;
         nFrees++;

         if (PRINT_OPERATIONS)
         {
            fprintf(stdout,
                    ">NODE\t#%03d (%c)\t%p\t%d\t[%s:%d] => freed\n",
                    i, ctype[x.type], x.pMem, x.nBytes, x.pszFile, x.nLine);
         }
      }
      i++;
   }
}

//============================================================================

void *mallocb(size_t size, char *pszFile, int nLine)
{
   void *pRet = malloc(size);

   AddAllocDataBytes(pRet, size, pszFile, nLine, MALLOC);

   if (PRINT_OPERATIONS)
   {
      fprintf(stdout,
              ">M\t%p\t%ld\t[%s:%d]\n", pRet, size, pszFile, nLine);
   }

   return(pRet);
}

//----------------------------------------------------------------------------

void *callocb(size_t num, size_t size, char *pszFile, int nLine)
{
   void *pRet = calloc(num, size);

   AddAllocDataBytes(pRet, num*size, pszFile, nLine, CALLOC);

   if (PRINT_OPERATIONS)
   {
      fprintf(stdout,
              ">C\t%p\t%ld\t(%ld x %ld)\t[%s:%d]\n",
              pRet, num*size, num, size, pszFile, nLine);
   }

   return(pRet);
}

//----------------------------------------------------------------------------

void *reallocb(void *memblock, size_t size, char *pszFile, int nLine)
{
   void *pRet;
   int   nBytes = SubtractAllocDataBytes(memblock, NULL, NULL);

   if (memblock && !nBytes)
   {
      if (PRINT_ERRORS)
      {
         fprintf(stderr,
                 ">ERROR\tTrying to free unallocated memory while"
                 " reallocating: %p [%s:%d]\n", memblock, pszFile, nLine);
      }
   }

   pRet = realloc(memblock, size);

   AddAllocDataBytes(pRet, size, pszFile, nLine, REALLOC);

   if (PRINT_OPERATIONS)
   {
      if (memblock)
      {
         fprintf(stdout,
                 ">R\t%p\t=>\t%p\t%d => %ld\t[%s:%d]\n",
                 memblock, pRet, nBytes, size, pszFile, nLine);
      }
      else
      {
         fprintf(stdout,
                 ">R(M)\t%p\t%ld\t[%s:%d]\n", pRet, size, pszFile, nLine);
      }
   }

   return(pRet);
}

//----------------------------------------------------------------------------

void freeb(void *memblock, char *pszFile, int nLine)
{
   char *pszOldFile = NULL;
   int   nOldLine = 0;
   int   nBytes;

   if (memblock)
   {
      nBytes = SubtractAllocDataBytes(memblock, &pszOldFile, &nOldLine);

      if (PRINT_OPERATIONS)
      {
         if (nBytes)
         {
            fprintf(stdout,
                    ">F\t%p\t%d\t(%s:%d)\t[%s:%d]\n",
                    memblock, nBytes, pszOldFile, nOldLine, pszFile, nLine);
         }
         else
         {
            fprintf(stdout,
                    ">F\t%p\t%d\t\t[%s:%d]\n",
                    memblock, nBytes, pszFile, nLine);
         }
      }

      if (!nBytes)
      {
         if (PRINT_ERRORS)
         {
            fprintf(stderr,
                    ">ERROR\tTrying to free unallocated memory: %p [%s:%d]\n",
                    memblock, pszFile, nLine);
         }
      }
   }

   free(memblock);
}

//----------------------------------------------------------------------------
