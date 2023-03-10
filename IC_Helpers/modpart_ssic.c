/*
 ** modpart_ssic.c -- VH 22-05-2014
 ** =====
 ** Shifts all coordinates in pkdgrav_planet IC file by a random amount.
 */

#include "ssio.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main(int argc, char*argv[]){

    SSIO ssio_in;
    SSIO ssio_out;

    SSHEAD head;
    SSDATA data;

    char inputfilename[160];
    char outputfilename[160];

    int i;

    double factor;
    double shift_x;
    double shift_y;
    double shift_z;

    /* RNG seed */
    srand(pow(time(NULL), 2));

    sprintf(inputfilename, "ssic.ss");
    sprintf(outputfilename, "ssic_shifted.ss");

    printf("Input: %s\n", inputfilename);
    printf("Output: %s\n", outputfilename);

    /* Open input file */
    if (ssioOpen(inputfilename, &ssio_in, SSIO_READ)) {
        printf("Could not open input file: %s\n", inputfilename);
        exit(1);
    }

    /* Open output file */
    if (ssioOpen(outputfilename, &ssio_out, SSIO_WRITE)) {
        printf("Could not open output file: %s\n", outputfilename);
        exit(1);
    }

    /* Open header from input */
    if (ssioHead(&ssio_in, &head)) {
        printf("Could not read header of input file: %s\n", inputfilename);
        exit(1);
    }

    /* Write header for output */
    ssioHead(&ssio_out, &head);

    /* Print header */
    /* printf("%.20g %d %d %.20g %.20g\n\n", head.time, head.n_data, head.n_planets, head.dEcoll, head.dSunMass); */

    /* Loop data */
    for (i = 0; i < head.n_data; ++i){
        /* Load input data */
        if (ssioData(&ssio_in, &data)) {
            printf("Could not read data of input file: %s\n", inputfilename);
            exit(1);
        }

        /* Print data before shifting */
        /* printf("%.10g %d %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n", head.time, data.org_idx, data.mass, data.radius, data.pos[0], data.pos[1], data.pos[2], data.vel[0], data.vel[1], data.vel[2], data.spin[0], data.spin[1], data.spin[2]); */

        /* Randomly shift coordinates */
        /* With (double), the minimum factor is 1.0e-16; */
        /*
         * 
         * A number like 0.1xxx can store 16 significant digits.
         * A number like 1.xxxx can store 15 significant digits.
         *
         */
        factor = 1.0e-15;
        /*
        shift_x = ((((double)rand()/(double)(RAND_MAX)) - 0.5) * factor) + 1;
        shift_y = ((((double)rand()/(double)(RAND_MAX)) - 0.5) * factor) + 1;
        shift_z = ((((double)rand()/(double)(RAND_MAX)) - 0.5) * factor) + 1;
        data.pos[0] *= shift_x;
        data.pos[1] *= shift_y;
        data.pos[2] *= shift_z;
        */

        /* Randomly add or subtract factor */
        /* With (double), the minimum factor is 1.0e-16; */
        /* X */
        if (rand() > RAND_MAX/2.0) {
            shift_x = factor;
        } else {
            shift_x = -factor;
        }
        /* Y */
        if (rand() > RAND_MAX/2.0) {
            shift_y = factor;
        } else {
            shift_y = factor;
        }
        /* Z */
        if (rand() > RAND_MAX/2.0) {
            shift_z = factor;
        } else {
            shift_z = -factor;
        }
        data.pos[0] += shift_x;
        data.pos[1] += shift_y;
        data.pos[2] += shift_z;
        
        /* Print data after shifting */
        /* printf("%.10g %d %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n***\n", head.time, data.org_idx, data.mass, data.radius, data.pos[0], data.pos[1], data.pos[2], data.vel[0], data.vel[1], data.vel[2], data.spin[0], data.spin[1], data.spin[2]); */

        /* Write output data */
        (void) ssioData(&ssio_out, &data);

    }

    /* Close files */
    (void) ssioClose(&ssio_out);
    (void) ssioClose(&ssio_in);

    return 0;

} 
