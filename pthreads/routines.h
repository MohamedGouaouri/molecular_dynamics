#ifndef H_PTHREAD_ROUTINES
#define H_PTHREAD_ROUTINES

#include <pthread.h>

struct MD_nullifyAccsTask {
    int start;
    int end;
};

struct MD_calcNewAccsTask {
    int start1;
    int end1;
    int start2;
    int end2;
};


void *nullifyAccsRoutine(void *arg);
void *calcNewAccsRoutine(void *arg);

#endif