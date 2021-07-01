#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <ViennaRNA/mfe.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils/basic.h>
 
int main()
{
  /* The RNA sequence */
  char  *seq = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA";
 
  /* allocate memory for MFE structure (length + 1) */
  char  *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
 
  /* predict Minmum Free Energy and corresponding secondary structure */
  float mfe = vrna_fold(seq, structure);
 
  /* print sequence, structure and MFE */
  printf("%s\n%s [ %6.2f ]\n", seq, structure, mfe);
 
  /* cleanup memory */
  free(structure);
 
  return 0;
}


/*
After 4 hours, I got an example to work with 
gcc soft_constraints_up.c -o soft_constraints -pthread -I/usr/local/include -I/usr/local/include/ViennaRNA -L/usr/local/lib -fno-lto -Wl,-fno-lto -lRNA -fopenmp -lpthread -lstdc++ -lm -ldl
*/