#ifndef H_PTHREAD_ROUTINES
#define H_PTHREAD_ROUTINES

#include <pthread.h>
const int MAXPART = 5001;

struct MD_VelocityVerlet_task
{
    int start;
    int end;
    double dt;
    double computation;
};

struct MD_Kinetic_task{
    int start ;
    int end ;
    double m;
    double* velocity[MAXPART][3];
};

void *updatePositionRoutine(void *arg);
void *updateVelocitiesRoutine(void *arg);
void *elasticWallsRoutine(void *arg);
// Kinetic
void *calculatePartialKinetic(void *arg);

#endif