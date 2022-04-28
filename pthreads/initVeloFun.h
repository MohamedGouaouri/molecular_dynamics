#ifndef H_PTHREAD_VELOFUN
#define H_PTHREAD_VELOFUN
#include <pthread.h>

struct MD_initVelocities_task {
    int start_index;
    int stop_index;
};

void *initGaussMat(void *arg);
void *masCenterInit(void *arg);
void *nullifyCenter(void *arg);
void *scaleAvgVeloc(void *arg);
void *lambdaProduct(void *arg);

#endif