#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "initVeloFun.h"
const int MAXPART = 5001;
extern double v[MAXPART][3];
extern double vCM[3];
extern double vSqdSum;
extern double lambda;
extern double m;

//  Numerical recipes Gaussian distribution number generator
double gaussdist()
{
    static bool available = false;
    static double gset;
    double fac, rsq, v1, v2;
    if (!available)
    {
        do
        {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;

        return v2 * fac;
    }
    else
    {

        available = false;
        return gset;
    }
}

// parallelize the first nested loop: initialization of v
void *initGaussMat(void *arg){
    struct MD_initVelocities_task *s = (struct MD_initVelocities_task *)arg;
    int start = s->start_index;
    double stop = s->stop_index;
    for(int i=start;i<stop;i++){
        for (int j = 0; j < 3; j++)
        {
            v[i][j] = gaussdist();
        }
    }
    pthread_exit(NULL);
}

//parallelize the second nested loop
void *masCenterInit(void *arg){
    struct MD_initVelocities_task *s = (struct MD_initVelocities_task *)arg;
    int start = s->start_index;
    double stop = s->stop_index;
    for(int i=start;i<stop;i++) {
        for (int j = 0; j < 3; j++) {
            vCM[j] += m * v[i][j];
        }
    }
    pthread_exit(NULL);
}
//Third nested loop
void *nullifyCenter(void *arg){
    struct MD_initVelocities_task *s = (struct MD_initVelocities_task *)arg;
    int start = s->start_index;
    double stop = s->stop_index;
    for(int i=start;i<stop;i++) {
        for (int j = 0; j < 3; j++) {
            v[i][j] -= vCM[j];
        }
    }
    pthread_exit(NULL);
}

//Fourth nested loop
void *scaleAvgVeloc(void *arg){
    struct MD_initVelocities_task *s = (struct MD_initVelocities_task *)arg;
    int start = s->start_index;
    double stop = s->stop_index;
    for(int i=start;i<stop;i++) {
        for (int j = 0; j < 3; j++) {
            vSqdSum += v[i][j] * v[i][j];
        }
    }
    pthread_exit(NULL);
}

//Fifth nested loop
void *lambdaProduct(void *arg){
    struct MD_initVelocities_task *s = (struct MD_initVelocities_task *)arg;
    int start = s->start_index;
    double stop = s->stop_index;
    for(int i=start;i<stop;i++) {
        for (int j = 0; j < 3; j++) {
            v[i][j] *= lambda;
        }
    }
    pthread_exit(NULL);
}