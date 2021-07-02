#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils/basic.h>

int main(){
  /* The RNA sequence */
  char      *seq = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA";

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

  for (ptr = pair_probabilities; ptr->i != 0; ptr++){
    a = (*ptr).i;
    b = (*ptr).j;
    c = (*ptr).p;
    unpairdness[a-1] += (double)c;
    unpairdness[b-1] += (double)c;
    if(a==38 || b == 38){
      printf("%g\n", c);
    }
  }

  printf("%g\n", unpairdness[37]);
  /* cleanup memory */
  free(pair_probabilities);
  free(propensity);

  return 0;
}


/*
After 4 hours, I got an example to work with 
gcc soft_constraints_up.c -o soft_constraints -pthread -I/usr/local/include -I/usr/local/include/ViennaRNA -L/usr/local/lib -fno-lto -Wl,-fno-lto -lRNA -fopenmp -lpthread -lstdc++ -lm -ldl
*/