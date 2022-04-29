#ifndef H_PTHREAD_ROUTINES
#define H_PTHREAD_ROUTINES

#include <pthread.h>
#include <math.h>
#define MAXPART 5001

struct MD_Init_task
{
    int begin;
    int end;
};

void *ToDoInit(void *arg);

struct MD_initVelocities_task
{
    int start_index;
    int stop_index;
};

void *initGaussMat(void *arg);
void *masCenterInit(void *arg);
void *nullifyCenter(void *arg);
void *scaleAvgVeloc(void *arg);
void *lambdaProduct(void *arg);

struct MD_VelocityVerlet_task
{
    int start;
    int end;
    double dt;
    double computation;
};

struct MD_Potential_Task
{
    int start;
    int end;
    double Pot;
};
void *potentialRoutine(void *arg);

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