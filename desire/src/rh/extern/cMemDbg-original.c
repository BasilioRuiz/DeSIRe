
/*************************************************************************************************
	This file is part of cMemDbg - Easy to use C memory leak detection library

	Copyright (C) 2009 Ezequiel Gaston Miravalles

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
*************************************************************************************************/


/*************************************************************************************************
Software: cMemDbg
File: cMemDbg.c
Version: 1.0
Last modification: 8/4/2009
Author: Ezequiel Gaston Miravalles
Website: http://www.neoegm.com/software/cmemdbg/
License: GNU GPL v3 (attached, read above)

Notes: Read cMemDbg.h notes.
*************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>

//Configuration Constants ***************
#define	PRINT_OPERATIONS	1
#define MAX_ALLOC			256
#define PRINT_OUTPUT		stdout
//***************************************

int g_nAllocMem = 0;
struct DataAlloc
{
	void *pData;
	int nBytes;
	char *pszFile;
	int nLine;
} g_pAllocData[MAX_ALLOC] = {0};

void PrintTotalAllocatedMemory()
{
	fprintf(PRINT_OUTPUT, ">INFO\tTotal allocated memory at this moment: %d\n", g_nAllocMem);
}

void PrintMemoryReservedByCMemDbgLibrary()
{
	fprintf(PRINT_OUTPUT, ">INFO\tMemory reserved by cMemDbg Library (constant): %d bytes\n", sizeof(g_nAllocMem) + sizeof(g_pAllocData));
}

int FreeAllocDataBytes(void *pPointer, char **ppszFile, int *pnLine)
{
	int i;

	if (!pPointer)
		return 0;

	for (i = 0; i < MAX_ALLOC; i++)
		if (g_pAllocData[i].pData == pPointer)
		{
			g_pAllocData[i].pData = 0;		//Free position

			g_nAllocMem -= g_pAllocData[i].nBytes;

			if (ppszFile)
				*ppszFile = g_pAllocData[i].pszFile;

			if (pnLine)
				*pnLine = g_pAllocData[i].nLine;

			return g_pAllocData[i].nBytes;
		}

	return 0;
}

void AddAllocDataBytes(void *pPointer, size_t size, char *pszFile, int nLine)
{
	int i;

	if (!pPointer)
		fprintf(PRINT_OUTPUT, ">MEM_ERROR\tCould not allocate %d bytes [%s:%d]\n", size, pszFile, nLine);

	for (i = 0; i < MAX_ALLOC; i++)
		if (!g_pAllocData[i].pData)
		{
			g_pAllocData[i].pData = pPointer;
			g_pAllocData[i].nBytes = size;
			g_pAllocData[i].pszFile = pszFile;
			g_pAllocData[i].nLine = nLine;
			g_nAllocMem += size;
			return;
		}

	fprintf(PRINT_OUTPUT, ">INTERNAL_ERROR\tAllocation stack overflow, please increase MAX_ALLOC\n");
}

void DumpUnfreedBlocks()
{
	int i;

	for (i = 0; i < MAX_ALLOC; i++)
		if (g_pAllocData[i].pData)
			fprintf(PRINT_OUTPUT, ">INFO\tUnfreed block\t%p\t%d\t\t[%s:%d]\n", g_pAllocData[i].pData, g_pAllocData[i].nBytes, g_pAllocData[i].pszFile, g_pAllocData[i].nLine);
}

void * mallocb(size_t size, char *pszFile, int nLine)
{
	void * pRet = malloc(size);

	if (PRINT_OPERATIONS)
		fprintf(PRINT_OUTPUT, ">A\t%p\t%d\t[%s:%d]\n", pRet, size, pszFile, nLine);

	AddAllocDataBytes(pRet, size, pszFile, nLine);

	return pRet;
}

void * callocb(size_t num, size_t size, char *pszFile, int nLine)
{
	void * pRet = calloc(num, size);

	if (PRINT_OPERATIONS)
		fprintf(PRINT_OUTPUT, ">CA\t%p\t%d\t(%d x %d)\t[%s:%d]\n", pRet, num*size, num, size, pszFile, nLine);

	AddAllocDataBytes(pRet, num*size, pszFile, nLine);

	return pRet;
}

void freeb(void *memblock, char *pszFile, int nLine)
{
	int nBytes;
	char * pszOldFile = NULL;
	int nOldLine = 0;

	if (memblock)
	{
		nBytes = FreeAllocDataBytes(memblock, &pszOldFile, &nOldLine);

		if (PRINT_OPERATIONS)
			if (nBytes)
				fprintf(PRINT_OUTPUT, ">F\t%p\t%d\t(%s:%d)\t[%s:%d]\n", memblock, nBytes, pszOldFile, nOldLine, pszFile, nLine);
			else
				fprintf(PRINT_OUTPUT, ">F\t%p\t%d\t\t[%s:%d]\n", memblock, nBytes, pszFile, nLine);

		if (!nBytes)
			fprintf(PRINT_OUTPUT, ">ERROR\tTrying to free unallocated memory: %p [%s:%d]\n", memblock, pszFile, nLine);
	}

	free(memblock);
}

void * reallocb(void *memblock, size_t size, char *pszFile, int nLine)
{
	void * pRet;
	int nBytes = FreeAllocDataBytes(memblock, NULL, NULL);

	if (memblock && !nBytes)
		fprintf(PRINT_OUTPUT, ">ERROR\tTrying to free unallocated memory while reallocating: %p [%s:%d]\n", memblock, pszFile, nLine);

	pRet = realloc(memblock, size);

	AddAllocDataBytes(pRet, size, pszFile, nLine);

	if (PRINT_OPERATIONS)
		if (memblock)
			fprintf(PRINT_OUTPUT, ">R\t%p\t=>\t%p\t%d\t=>\t%d\t[%s:%d]\n", memblock, pRet, nBytes, size, pszFile, nLine);
		else
			fprintf(PRINT_OUTPUT, ">R(A)\t%p\t%d\t[%s:%d]\n", pRet, size, pszFile, nLine);

	return pRet;
}

void PrintMemoryLeakInfo()
{
	if (g_nAllocMem == 0)
		fprintf(PRINT_OUTPUT, ">INFO\tNo memory leaks detected\n");
	else if (g_nAllocMem < 0)
		fprintf(PRINT_OUTPUT, ">INFO\tPROBLEM: There was more freed memory than allocated (%d bytes) [This shouldn't happen, memory corruption?]\n", g_nAllocMem);
	else
		fprintf(PRINT_OUTPUT, ">INFO\tPROBLEM: Memory leak found (%d bytes)\n", g_nAllocMem);

	DumpUnfreedBlocks();
}
