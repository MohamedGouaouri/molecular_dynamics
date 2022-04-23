#ifndef H_PTHREAD_ROUTINES
#define H_PTHREAD_ROUTINES

#include <pthread.h>

struct MD_VelocityVerlet_task
{
    int start;
    int end;
    double dt;
    double computation;
};

void *updatePositionRoutine(void *arg);
void *updateVelocitiesRoutine(void *arg);
void *elasticWallsRoutine(void *arg);

#endif