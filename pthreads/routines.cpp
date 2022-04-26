#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "routines.h"

const int MAXPART = 5001;

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
