#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/constraints/basic.h>
#include <ViennaRNA/constraints/hard.h>

void deltaG(int arraySize, double *array)
{
    int i;
    char sequence [arraySize];
    char unconstrained[arraySize];
    char cstruc[arraySize+1];
    double deltaG[arraySize];
    for (i = 0; i < arraySize; i++){
        deltaG[i] = 0;
    }

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
        unconstrained[i] = '.';
        if(i <15){
            cstruc[i] = 'x';
        }else{
            cstruc[i] = '.';
        }
    }

    sequence[arraySize] = '\0';
    cstruc[arraySize] = '\0';

    char *seq = sequence;
    printf("%s",sequence);

    vrna_fold_compound_t  *fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);

    unsigned int constraint_options = VRNA_CONSTRAINT_DB_DEFAULT;

    char *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq)+1));
    vrna_constraints_add(fc, (const char *)unconstrained, constraint_options);
    float mfe = vrna_mfe(fc, structure);
    float constrainedMFE;

    for(i = 0; i<arraySize - 15; i++){
        cstruc[i+14]  = 'x';
        vrna_constraints_add(fc, (const char *)cstruc, constraint_options);
        constrainedMFE = vrna_mfe(fc, structure);
        deltaG[i] = constrainedMFE;
        vrna_hc_init(fc);
        cstruc[i] = '.';
    }

    deltaG[arraySize-1] = (double)mfe;

    for (i = 0; i < arraySize; i++)
        array[i] = deltaG[i];
    
    //char cstruc[21] = "...((((......))))..."; //gives -2.5 in rnafold, 
    //char cstruc[21] = ".....xxx............"; //returns same mfe, but different struc
    //char cstruc[21] = "xxxxxxxxxxxxxxxxxxxx"; // returns 0

    free(fc);
    free(structure);
}

// gcc -shared -fPIC -o deltaG.so deltaG.c -pthread -I/usr/local/include -I/usr/local/include/ViennaRNA -L/usr/local/lib -fno-lto -Wl,-fno-lto -lRNA -fopenmp -lpthread -lstdc++ -lm -ldl