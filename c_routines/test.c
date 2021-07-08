#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils/basic.h>

int main(){
  /* The RNA sequence */
  char      *seq = "CUCGCAGACUAACA";

  /* allocate memory for pairing propensity string (length + 1) */
  char      *propensity = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));

  /* pointers for storing and navigating through base pair probabilities */
  vrna_ep_t *ptr, *pair_probabilities = NULL;

  float     en = vrna_pf_fold(seq, propensity, &pair_probabilities);

  /* print sequence, pairing propensity string and ensemble free energy */
  double unpairdness[strlen(seq)];
  int a;
  int b;
  double c;

  a = 0;
  b =  0;
  c = 0;
  for(int i = 0; i<strlen(seq); i++){
      unpairdness[i] = 0.0;
  }  

  // for (ptr = pair_probabilities; ptr->i != 0; ptr++){
  //     printf("p(%d,%d) = %g\n", ptr->i,ptr->j,ptr->p);
  // }
  printf("starting at %g\n",unpairdness[9]);

  for (ptr = pair_probabilities; ptr->i != 0; ptr++){
    a = (*ptr).i;
    b = (*ptr).j;
    c = (*ptr).p;
    if(a == 10){
      printf("adding %g\n", (double)c);
    }
    if(b == 10){
      printf("adding %g\n", (double)c);
    }
    unpairdness[a-1] += (double)c;
    unpairdness[b-1] += (double)c;
  }

  printf("%g\n", unpairdness[9]);

  // for(int i = 0; i<strlen(seq); i++){
  //     unpairdness[i] = 1-unpairdness[i];
  //     printf("%g\n", unpairdness[i]);
  // }

  /* cleanup memory */

  free(pair_probabilities);
  free(propensity);

  return 0;
}


/*
After 4 hours, I got an example to work with 
gcc test.c -o test -pthread -I/usr/local/include -I/usr/local/include/ViennaRNA -L/usr/local/lib -fno-lto -Wl,-fno-lto -lRNA -fopenmp -lpthread -lstdc++ -lm -ldl
*/