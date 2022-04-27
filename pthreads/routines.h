#ifndef H_PTHREAD_ROUTINES
#define H_PTHREAD_ROUTINES

#include <pthread.h>
#define MAXPART 5001

struct MD_VelocityVerlet_task
{
    int start;
    int end;
    double dt;
    double computation;
};

struct MD_Kinetic_task
{
    int start;
    int end;
    double m;
};

void *updatePositionRoutine(void *arg);
void *updateVelocitiesRoutine(void *arg);
void *elasticWallsRoutine(void *arg);
struct MD_nullifyAccsTask
{
    int start;
    int end;
};

struct MD_calcNewAccsTask
{
    int start1;
    int end1;
    int start2;
    int end2;
};

void *nullifyAccsRoutine(void *arg);
void *calcNewAccsRoutine(void *arg);
// Kinetic
void *calculatePartialKinetic(void *arg);

#endif