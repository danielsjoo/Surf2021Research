#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils/basic.h>

void sortArray(int arraySize, double *array)
{
    // int temp, i, j;
    // for (i = 0; i < arraySize; i++)
    //     for (j = i + 1; j < arraySize; j++)
    //         if (array[i] > array[j]) {
    //             temp = array[i];
    //             array[i] = array[j];
    //             array[j] = temp;
    //         }
    //key: 1:A 2:G 3:C 4:T
    
    int i;
    char sequence [arraySize];

    for (i = 0; i < arraySize; i++){
        if((int)array[i] == 1){
            sequence[i] = 'A';
        }else if((int)array[i] == 2){
            sequence[i] = 'G';
        }else if((int)array[i] == 3){
            sequence[i] = 'C';
        }else {
            sequence[i] = 'U';
        }
    }

    char *seq = sequence;

    char      *propensity = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));

    /* pointers for storing and navigating through base pair probabilities */
        vrna_ep_t *ptr, *pair_probabilities = NULL;

    float     en = vrna_pf_fold(seq, propensity, &pair_probabilities);

    /* print sequence, pairing propensity string and ensemble free energy */
    double unpairdness[arraySize];
    int a;
    int b;
    double c;

    for (ptr = pair_probabilities; ptr->i != 0; ptr++){
        a = (*ptr).i;
        b = (*ptr).j;
        c = (*ptr).p;
        if (c > 0.0001)
            unpairdness[a-1] += (double)c;
            unpairdness[b-1] += (double)c;
    }

    for(int i = 0; i<strlen(seq); i++){
        unpairdness[i] = 1-unpairdness[i];
    }

    for (i = 0; i < arraySize; i++)
        array[i] = unpairdness[i];

    /* cleanup memory */
    free(pair_probabilities);
    free(propensity);
}
