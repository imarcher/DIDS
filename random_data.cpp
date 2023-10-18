#include "random_data.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <readline/readline.h>
#include <getopt.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <zconf.h>
#include <algorithm>

#define STD 1   // Standard deviation

static void z_normalize(float *ts, int size) {
    int i;
    float mean = 0;//gsl_stats_mean(ts, 1, size);
    float std = 0;//gsl_stats_sd(ts, 1, size);
    for (i=0; i<size; i++)
    {
        mean += ts[i];
    }
    mean /= size;

    for (i=0; i<size; i++)
    {
        std += (ts[i] - mean) * (ts[i] - mean);
    }
    std /= size;
    std = sqrt(std);
    for (i = 0; i < size; i++)
    {
        ts[i] = (ts[i] - mean) / std;
    }
}

static float * generate (float *ts, int size, gsl_rng * r, char normalize) {
    int i;
    float x = 0, dx;

    for (i = 0; i < size; i++)
    {
        dx = gsl_ran_gaussian (r, STD); // mean=0, std=STD
        x += dx;
        ts[i] = x;
    }

    if(normalize == 1)
    {
        z_normalize(ts, size);
    }
    return ts;
}

/**
    Generates a set of random time series.
**/
static void generate_random_timeseries(int length, int number_of_timeseries,
                                char normalize, const char * filename, int seed) {
    // Initialize random number generation
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
//    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    gsl_rng_default_seed = (seed);
    r = gsl_rng_alloc (T);
    FILE * data_file = fopen (filename,"w");


    size_t *p;


    float *ts = (float*)malloc(sizeof(float) * length);

    int i;
    for (i=1; i<=number_of_timeseries; i++)
    {
        generate(ts, length, r, normalize);
        fwrite(ts, sizeof(float), length,data_file);
        if(i % 1000 == 0) {
            fprintf(stderr,"\r\x1b[m>> Generating: \x1b[36m%2.2lf%%\x1b[0m",(float) ((float)i/(float)number_of_timeseries) * 100);
        }
    }
    fprintf(stderr, "\n");

    // Finalize random number generator
    fclose (data_file);
    gsl_rng_free (r);
}

void RandomData::generate_data(uint64_t ts_length, uint64_t ts_num, const std::string &filename_) {

    char normalize = 1;
    const char * filename = filename_.c_str();
    int seed = 0;


    fprintf(stderr, ">> generate random time series...\n");
    fprintf(stderr, ">> data filename: %s\n", filename);
    generate_random_timeseries(ts_length, ts_num, normalize, filename, seed);

    fprintf(stderr, ">> finished\n");
}