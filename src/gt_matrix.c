#include <stdio.h>
#include <unistd.h>
#include <linux/unistd.h>
#include <sys/syscall.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sched.h>
#include <signal.h>
#include <setjmp.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include "gt_include.h"


#define ROWS 256
#define COLS ROWS
#define SIZE COLS

#define NUM_CPUS 2
#define NUM_GROUPS NUM_CPUS
#define PER_GROUP_COLS (SIZE/NUM_GROUPS)

#define NUM_THREADS 128
#define PER_THREAD_ROWS (SIZE/NUM_THREADS)


/* A[SIZE][SIZE] X B[SIZE][SIZE] = C[SIZE][SIZE]
 * Let T(g, t) be thread 't' in group 'g'. 
 * T(g, t) is responsible for multiplication : 
 * A(rows)[(t-1)*SIZE -> (t*SIZE - 1)] X B(cols)[(g-1)*SIZE -> (g*SIZE - 1)] */
int sched_fag;
unsigned long int run_time[128];
unsigned long int mean_run_time[16];
double sd_run_time[16];
unsigned long int exe_time[128];
unsigned long int mean_exe_time[128];
double sd_exe_time[128];
static flage = 0;
extern void gt_yield();
typedef struct matrix
{
	int m[SIZE][SIZE];

	int rows;
	int cols;
	unsigned int reserved[2];
} matrix_t;


typedef struct __uthread_arg
{
	matrix_t *_A, *_B, *_C;
	unsigned int reserved0;

	unsigned int tid;
	unsigned int gid;
	int start_row; /* start_row -> (start_row + PER_THREAD_ROWS) */
	int start_col; /* start_col -> (start_col + PER_GROUP_COLS) */
    int sz;
	int credits;
}uthread_arg_t;
	
struct timeval tv1;

static void generate_matrix(matrix_t *mat, int val, int sz)
{

	int i,j;
	mat->rows = SIZE;
	mat->cols = SIZE;
	for(i = 0; i < sz;i++)
		for( j = 0; j < sz; j++ )
		{
			mat->m[i][j] = val;
		}
	return;
}

static void print_matrix(matrix_t *mat, int sz)
{
	int i, j;

	for(i=0;i<sz;i++)
	{
		for(j=0;j<sz;j++)
			printf(" %d ",mat->m[i][j]);
		printf("\n");
	}

	return;
}

static void * uthread_mulmat(void *p) //changed
{
	int i, j, k;
	int start_row, end_row;
	int start_col, end_col;
	unsigned int cpuid;
	struct timeval tv2;

#define ptr ((uthread_arg_t *)p)

	i=0; j= 0; k=0;

	start_row = ptr->start_row;
	//end_row = (ptr->start_row + PER_THREAD_ROWS);

#ifdef GT_GROUP_SPLIT
	start_col = ptr->start_col;
	end_col = (ptr->start_col + PER_THREAD_ROWS);
#else
	start_col = 0;
	end_col = ptr->sz; //changed
#endif

#ifdef GT_THREADS
	cpuid = kthread_cpu_map[kthread_apic_id()]->cpuid;
	fprintf(stderr, "\nThread(id:%d, group:%d, cpu:%d) started",ptr->tid, ptr->gid, cpuid);
#else
	fprintf(stderr, "\nThread(id:%d, group:%d) started",ptr->tid, ptr->gid);
#endif

	for(i = 0; i < ptr->sz; i++)
		for(j = 0; j < ptr->sz; j++)
			for(k = 0; k < ptr->sz; k++) //changed
            {ptr->_C->m[i][j] += ptr->_A->m[i][k] * ptr->_B->m[k][j]; if(flage==0&& ptr->tid==50) {flage=1;gt_yield();}}

#ifdef GT_THREADS
	fprintf(stderr, "\nThread(id:%d, group:%d, cpu:%d) finished (TIME : %lu s and %lu us)",
			ptr->tid, ptr->gid, cpuid, (tv2.tv_sec - tv1.tv_sec), (tv2.tv_usec - tv1.tv_usec));
#else
	gettimeofday(&tv2,NULL);
	fprintf(stderr, "\nThread(id:%d, group:%d) finished (TIME : %lu s and %lu us)",
			ptr->tid, ptr->gid, (tv2.tv_sec - tv1.tv_sec), (tv2.tv_usec - tv1.tv_usec));
    run_time[ptr->tid]=((tv2.tv_sec - tv1.tv_sec)*1000L)+ ((tv2.tv_usec - tv1.tv_usec)/1000L);
#endif

#undef ptr
	return 0;
}

matrix_t A, B, C;

static void init_matrices(int sz)
{

	generate_matrix(&A, 1, sz);
	generate_matrix(&B, 1, sz);
	generate_matrix(&C, 0, sz);

	return;
}


uthread_arg_t uargs[NUM_THREADS];
uthread_t utids[NUM_THREADS];

int main(int argc, char* argv[])
{
	sched_fag = atoi(argv[1]);
	uthread_arg_t *uarg;

	int inx,c_i,j,sz,credits;
	gtthread_app_init();
    if (!sched_fag)
	init_matrices(SIZE);

	gettimeofday(&tv1,NULL);
    if (sched_fag) {
		for (inx = 0; inx < 4; inx++) {
			if (inx == 0) {
				sz = 32;
				init_matrices(sz);
			}
			if (inx == 1) {
				sz = 64;
				init_matrices(sz);
			}
			if (inx == 2) {
				sz = 128;
				init_matrices(sz);
			}
			if (inx == 3) {
				sz = 256;
				init_matrices(sz);
			}
			for (c_i = 0; c_i < 4; c_i++) {
				if (inx == 0) credits = 25;
				if (inx == 1) credits = 50;
				if (inx == 2) credits = 75;
				if (inx == 3) credits = 100;
				for (j = 0; j < 8; j++) {

					uarg = &uargs[(32 * inx) + (8 * c_i) + j];
					uarg->_A = &A;
					uarg->_B = &B;
					uarg->_C = &C;

					uarg->tid = (32 * inx) + (8 * c_i) + j;

					uarg->gid = (((32 * inx) + (8 * c_i) + j) % NUM_GROUPS);

					uarg->start_row = (((32 * inx) + (8 * c_i) + j) * PER_THREAD_ROWS);
#ifdef GT_GROUP_SPLIT
                    /* Wanted to split the columns by groups !!! */
                    uarg->start_col = (uarg->gid * PER_GROUP_COLS);
#endif
					uarg->start_row = 0;
					uarg->credits = credits;
					uarg->sz = sz;

					uthread_create(&utids[(32 * inx) + (8 * c_i) + j], uthread_mulmat, uarg, uarg->gid, credits);
                }
			}
		}

	}
    if (!sched_fag)
    {
        for(inx=0; inx<NUM_THREADS; inx++)

        {
            uarg = &uargs[inx];
            uarg->_A = &A;
            uarg->_B = &B;
            uarg->_C = &C;
            uarg->tid = inx;
            uarg->gid = (inx % NUM_GROUPS);
            uarg->sz=256;
            //uarg->start_row = (inx * PER_THREAD_ROWS);
            #ifdef GT_GROUP_SPLIT

            /* Wanted to split the columns by groups !!! */
            uarg->start_col = (uarg->gid * PER_GROUP_COLS);
            #endif
            uthread_create(&utids[inx], uthread_mulmat, uarg, uarg->gid,0);
        }
    }


	gtthread_app_exit();
    if(sched_fag){
    for(inx=0;inx<4;inx++)
    {
        for(c_i=0;c_i<4;c_i++)
        {
            for(j=0;j<8;j++)
            {
                mean_run_time[(inx*4)+c_i] = mean_run_time[(inx*4)+c_i]+run_time[(32*inx)+(8*c_i)+j];
                mean_exe_time[(inx*4)+c_i] = mean_exe_time[(inx*4)+c_i]+exe_time[(32*inx)+(8*c_i)+j];
            }
        }
    }
    for (j=0;j<16;j++) mean_run_time[j] = (float)mean_run_time[j]*(1/8.0);
        for(inx=0;inx<4;inx++)
        {
            for(c_i=0;c_i<4;c_i++)
            {
                for(j=0;j<8;j++)
                {
                    sd_run_time[(inx*4)+c_i] =  sd_run_time[(inx*4)+c_i] + pow(abs((float)run_time[(32*inx)+(8*c_i)+j]-(float)mean_run_time[(inx*4)+c_i]),2);
                    sd_exe_time[(inx*4)+c_i] =  sd_exe_time[(inx*4)+c_i] + pow(abs((float)exe_time[(32*inx)+(8*c_i)+j]-(float)mean_run_time[(inx*4)+c_i]),2);
                }
            }
        }
    for (j=0;j<16;j++)
        {
            //printf("The deviations are %g \n",sd_run_time[j]);

            sd_run_time[j] = sqrt(sd_run_time[j]/8.0);
            sd_exe_time[j] = sqrt(sd_exe_time[j]/8.0);
            //printf("SD is %f \n" ,sd_run_time[j]);
        }

        for (inx = 0; inx < 16; inx++)
        {

            printf("Mean Run Time for group %d are %f ms \n", inx,(float)mean_run_time[inx]);
            printf("Mean Execution Time for group %d are %f ms \n", inx,(float)mean_exe_time[inx]);
            printf("Standard Run Time Deviations for group %d are %f ms \n",inx,(float)sd_run_time[inx]);
            printf("Standard Execution Time Deviations for group %d are %f ms \n",inx,(float)sd_exe_time[inx]);
            printf("\n");
            printf("\n");

        }
    }
	 //print_matrix(&C,32);
	// fprintf(stderr, "********************************");
	return(0);
}
