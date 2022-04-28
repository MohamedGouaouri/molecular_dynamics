#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "routines.h"

extern double r[MAXPART][3];
//  Velocity
extern double v[MAXPART][3];
//  Acceleration
extern double a[MAXPART][3];
//  Force
extern double F[MAXPART][3];

extern double m;

extern double L;

extern double vCM[3];
extern double vSqdSum;
extern double lambda;

extern double gaussdist();

void *ToDoInit(void *arg)
{

    int n, p, i, j, k;
    double pos;

    //  spacing between atoms along a given direction
    n = int(ceil(pow(N, 1.0 / 3)));
    pos = L / n;

    //  index for number of particles assigned positions
    p = 0;

    struct MD_Init_task *initFun = (struct MD_Init_task *)arg;

    for (i = initFun->begin; i <= initFun->end; i++)
    {

        for (j = initFun->begin; j <= initFun->end; j++)
        {

            for (k = initFun->begin; k <= initFun->end; k++)
            {

                if (p < N)
                {
                    r[p][0] = (i + 0.5) * pos;
                    r[p][1] = (j + 0.5) * pos;
                    r[p][2] = (k + 0.5) * pos;
                }
                p++;
            }
        }
    }
}

// parallelize the first nested loop: initialization of v
void *initGaussMat(void *arg)
{
    struct MD_initVelocities_task *s = (struct MD_initVelocities_task *)arg;
    int start = s->start_index;
    double stop = s->stop_index;
    for (int i = start; i < stop; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            v[i][j] = gaussdist();
        }
    }
    pthread_exit(NULL);
}

// parallelize the second nested loop
void *masCenterInit(void *arg)
{
    struct MD_initVelocities_task *s = (struct MD_initVelocities_task *)arg;
    int start = s->start_index;
    double stop = s->stop_index;
    for (int i = start; i < stop; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vCM[j] += m * v[i][j];
        }
    }
    pthread_exit(NULL);
}
// Third nested loop
void *nullifyCenter(void *arg)
{
    struct MD_initVelocities_task *s = (struct MD_initVelocities_task *)arg;
    int start = s->start_index;
    double stop = s->stop_index;
    for (int i = start; i < stop; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            v[i][j] -= vCM[j];
        }
    }
    pthread_exit(NULL);
}

// Fourth nested loop
void *scaleAvgVeloc(void *arg)
{
    struct MD_initVelocities_task *s = (struct MD_initVelocities_task *)arg;
    int start = s->start_index;
    double stop = s->stop_index;
    for (int i = start; i < stop; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vSqdSum += v[i][j] * v[i][j];
        }
    }
    pthread_exit(NULL);
}

// Fifth nested loop
void *lambdaProduct(void *arg)
{
    struct MD_initVelocities_task *s = (struct MD_initVelocities_task *)arg;
    int start = s->start_index;
    double stop = s->stop_index;
    for (int i = start; i < stop; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            v[i][j] *= lambda;
        }
    }
    pthread_exit(NULL);
}

void *updatePositionRoutine(void *arg)
{
    struct MD_VelocityVerlet_task *t = (struct MD_VelocityVerlet_task *)arg;

    int n = t->end;
    double dt = t->dt;
    for (int i = t->start; i < n; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            r[i][j] += v[i][j] * dt + 0.5 * a[i][j] * dt * dt;

            v[i][j] += 0.5 * a[i][j] * dt;
        }
        // printf("  %i  %6.4e   %6.4e   %6.4e\n",i,r[i][0],r[i][1],r[i][2]);
    }
    pthread_exit(NULL);
}

void *updateVelocitiesRoutine(void *arg)
{
    struct MD_VelocityVerlet_task *t = (struct MD_VelocityVerlet_task *)arg;

    int n = t->end;
    double dt = t->dt;
    for (int i = t->start; i < n; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            v[i][j] += 0.5 * a[i][j] * dt;
        }
    }
    pthread_exit(NULL);
}

void *elasticWallsRoutine(void *arg)
{
    struct MD_VelocityVerlet_task *t = (struct MD_VelocityVerlet_task *)arg;
    int n = t->end;
    double dt = t->dt;
    double psum;
    for (int i = t->start; i < n; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (r[i][j] < 0.)
            {
                v[i][j] *= -1.;                               //- elastic walls
                t->computation += 2 * m * fabs(v[i][j]) / dt; // contribution to pressure from "left" walls
            }
            if (r[i][j] >= L)
            {
                v[i][j] *= -1.;                               //- elastic walls
                t->computation += 2 * m * fabs(v[i][j]) / dt; // contribution to pressure from "right" walls
            }
        }
    }
    pthread_exit(NULL);
}
//  Position
double r[MAXPART][3];
//  Acceleration
double a[MAXPART][3];

extern int N;

void *nullifyAccsRoutine(void *arg)
{

    struct MD_nullifyAccsTask *mytask = (struct MD_nullifyAccsTask *)arg;

    for (int i = mytask->start; i <= mytask->end; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            a[i][k] = 0;
        }
    }

    pthread_exit(NULL);
}

void *calcNewAccsRoutine(void *arg)
{

    struct MD_calcNewAccsTask *myTask = (struct MD_calcNewAccsTask *)arg;
    int startMarks[] = {myTask->start1, myTask->start2};
    int endMarks[] = {myTask->end1, myTask->end2};
    int nbiter = 2;

    int i, j, k;
    double f, rSqd;
    double rij[3]; // position of i relative to j

    for (size_t iter = 0; iter < nbiter - 1; iter++)
    {

        for (i = startMarks[iter]; i < endMarks[iter]; i++)
        {
            // loop over all distinct pairs i,j
            for (j = i + 1; j < N; j++)
            {
                // initialize r^2 to zero
                rSqd = 0;

                for (k = 0; k < 3; k++)
                {
                    //  component-by-componenent position of i relative to j
                    rij[k] = r[i][k] - r[j][k];
                    //  sum of squares of the components
                    rSqd += rij[k] * rij[k];
                }

                //  From derivative of Lennard-Jones with sigma and epsilon set equal to 1 in natural units!
                f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
                for (k = 0; k < 3; k++)
                {
                    //  from F = ma, where m = 1 in natural units!
                    a[i][k] += rij[k] * f;
                    a[j][k] -= rij[k] * f;
                }
            }
        }
    }

    pthread_exit(NULL);
}

// kinetic
void *calculatePartialKinetic(void *arg)
{
    struct MD_Kinetic_task *t = (struct MD_Kinetic_task *)arg;
    double v2, partialKin;
    partialKin = 0.;
    for (int i = t->start; i < t->end; i++)
    {

        v2 = 0.;
        for (int j = 0; j < 3; j++)
        {
            v2 += v[i][j] * v[i][j];
        }
        partialKin += t->m * v2 / 2.;
    }
    pthread_exit(&partialKin);
}
