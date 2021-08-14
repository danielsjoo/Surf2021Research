#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/constraints/basic.h>
#include <ViennaRNA/constraints/hard.h>

int main(
)
{
    int arraySize = 20;
    int i;
    //key: 1:A 2:G 3:C 4:T
    //sequence is UAGACGUAGCUAGACGUAGC
    //double is   4121324123
    //mfe is -2.5
    double deltaG [arraySize];
    double doubleSequence[20] = {4.0, 1.0, 2.0, 1.0, 3.0, 2.0, 4.0, 1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 1.0, 3.0, 2.0, 4.0, 1.0, 2.0, 3.0};
    double *array = doubleSequence;
    char sequence [arraySize+1];

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
    }

    printf("%s\n",sequence);

    char *seq = sequence;

    vrna_fold_compound_t  *fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);

    unsigned int constraint_options = VRNA_CONSTRAINT_DB_DEFAULT;

    char unconstrained[21] = "....................";
    vrna_constraints_add(fc, (const char *)unconstrained, constraint_options);
    char *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq)+1));
    float mfe = vrna_mfe(fc, structure);
    printf("%g\n", mfe);
    deltaG[arraySize] = mfe;
    float constrainedMFE;
    char cstruc[21] = "xxxxx...............";

    for(i = 0; i<arraySize - 5; i++){
        cstruc[i+5] = 'x';
        vrna_constraints_add(fc, (const char *)cstruc, constraint_options);
        constrainedMFE = vrna_mfe(fc, structure);
        vrna_hc_init(fc);
        array[i] = mfe - constrainedMFE;
        printf("%g : %s\n", constrainedMFE,(const char *)cstruc);
        cstruc[i] = '.';
    }

    char q[90];
    for (i = 0; i < 90; i++){

        if(i > 29 && i <45){
            q[i] = 'x';
        }else{
            q[i] = '.';
        }
    }

    printf("%s\n",q);
    //char cstruc[21] = "...((((......))))..."; //gives -2.6 in rnafold, 
    //char cstruc[21] = ".....xxx............"; //returns same mfe, but different struc
    //char cstruc[21] = "xxxxxxxxxxxxxxxxxxxx"; // returns 0

    free(fc);
    free(structure);

    return 1;
}

// gcc test2.c -o test2 -pthread -I/usr/local/include -I/usr/local/include/ViennaRNA -L/usr/local/lib -fno-lto -Wl,-fno-lto -lRNA -fopenmp -lpthread -lstdc++ -lm -ldl