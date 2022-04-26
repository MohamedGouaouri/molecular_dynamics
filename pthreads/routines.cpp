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