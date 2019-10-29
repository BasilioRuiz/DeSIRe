#include <stdio.h>

int copy_( char *origin, char *destin )
{
   /* En linux tenemos dos grupos de funciones para lectura y escritura de
    * ficheros. Las funciones open(), write(), read() y close() son de algo
    * mas bajo nivel y especificas de linux, asi que no las usaremos en este
    * ejemplo. En su lugar usaremos fopen(), fwrite(), fread() y fclose(),
    * que son estandar de C y estan presentes en todos los compiladores.
    *
    * Para nuestro ejemplo, abriremos el fichero original como lectura y el
    * fichero en el que haremos la copia como escritura. Puesto que nos da
    * igual si el fichero es de texto o no, lo abriremos como binario, ya que
    * es mas general (un fichero de texto, con fines de linea, no es mas que
    * un caso particular de un fichero binario, que tiene bytes).
    *
    * http://www.chuidiang.org/clinux/ficheros/fichero-binario.php
    */

   FILE *f1, *f2;        // descriptores de los ficheros
   char  buffer[1024];   // buffer de lectura
   int   leidos;         // numero de items leidos o -1 en error
   int   n = 0;          // numero total de bytes leidos y copiados

   // Apertura del fichero original, para lectura en binario.
   if ((f1 = fopen(origin, "rb")) == NULL)
   {
      perror("Error in function copy_");
      return(-1);
   }
   // Apertura del fichero de destino, para escritura en binario.
   if ((f2 = fopen ("fichero2.dat", "wb")) == NULL)
   {
      perror("Error in function copy_");
      return(-1);
   }

   // Leemos bytes hasta 1024.
   leidos = fread(buffer, 1, sizeof(buffer), f1);

   // Bucle mientras hayamos leido algo (si leido == 0 es fin de fichero).
   while (leidos > 0)
   {
      // Lo escribimos en el fichero destino.
      fwrite(buffer, 1, leidos, f2);
      // Incrementamos el numero de bytes escritos.
      n += leidos;
      // Y leemos el siguiente bloque.
      leidos = fread(buffer, 1, sizeof(buffer), f1);
   }

   // Cerramos los ficheros.
   fclose(f1);
   fclose(f2);

   // Control de error.
   if (leidos == -1)
   {
      perror("Error in function copy_");
      return(-1);
   }

   // Retornamos el numero de bytes escritos.
   return(n);
}

int remove_( const char *filename )
{
   return(remove(filename));
}

int rename_( const char *old, const char *new )
{
   return(rename(old, new));
}