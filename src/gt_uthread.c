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
#include <time.h>
#include "gt_include.h"
/**********************************************************************/
/** DECLARATIONS **/
/**********************************************************************/


/**********************************************************************/
/* kthread runqueue and env */

/* XXX: should be the apic-id */
#define KTHREAD_CUR_ID	0

/**********************************************************************/
/* uthread scheduling */
static void uthread_context_func(int);
static int uthread_init(uthread_struct_t *u_new);
static int flage=0;
static int flaga=0;
/**********************************************************************/
/* uthread creation */
#define UTHREAD_DEFAULT_SSIZE (16 * 1024)
void gt_yield();
extern int uthread_create(uthread_t *u_tid, int (*u_func)(void *), void *u_arg, uthread_group_t u_gid,int credit);
extern long int exe_time[128];
/**********************************************************************/
/** DEFNITIONS **/
/**********************************************************************/

/**********************************************************************/
/* uthread scheduling */

/* Assumes that the caller has disabled vtalrm and sigusr1 signals */
/* uthread_init will be using */
struct timeval tl1,tl2;
static int uthread_init(uthread_struct_t *u_new)
{
	stack_t oldstack;
	sigset_t set, oldset;
	struct sigaction act, oldact;

	gt_spin_lock(&(ksched_shared_info.uthread_init_lock));
	if(sched_fag==1)
	printf("My Credit is %d",u_new->credit);
	/* Register a signal(SIGUSR2) for alternate stack */
	act.sa_handler = uthread_context_func;
	act.sa_flags = (SA_ONSTACK | SA_RESTART);
	if(sigaction(SIGUSR2,&act,&oldact))
	{
		fprintf(stderr, "uthread sigusr2 install failed !!");
		return -1;
	}

	/* Install alternate signal stack (for SIGUSR2) */
	if(sigaltstack(&(u_new->uthread_stack), &oldstack))
	{
		fprintf(stderr, "uthread sigaltstack install failed.");
		return -1;
	}

	/* Unblock the signal(SIGUSR2) */
	sigemptyset(&set);
	sigaddset(&set, SIGUSR2);
	sigprocmask(SIG_UNBLOCK, &set, &oldset);


	/* SIGUSR2 handler expects kthread_runq->cur_uthread
	 * to point to the newly created thread. We will temporarily
	 * change cur_uthread, before entering the synchronous call
	 * to SIGUSR2. */

	/* kthread_runq is made to point to this new thread
	 * in the caller. Raise the signal(SIGUSR2) synchronously */
#if 0
	raise(SIGUSR2);
#endif
	syscall(__NR_tkill, kthread_cpu_map[kthread_apic_id()]->tid, SIGUSR2);

	/* Block the signal(SIGUSR2) */
	sigemptyset(&set);
	sigaddset(&set, SIGUSR2);
	sigprocmask(SIG_BLOCK, &set, &oldset);
	if(sigaction(SIGUSR2,&oldact,NULL))
	{
		fprintf(stderr, "uthread sigusr2 revert failed !!");
		return -1;
	}

	/* Disable the stack for signal(SIGUSR2) handling */
	u_new->uthread_stack.ss_flags = SS_DISABLE;

	/* Restore the old stack/signal handling */
	if(sigaltstack(&oldstack, NULL))
	{
		fprintf(stderr, "uthread sigaltstack revert failed.");
		return -1;
	}

	gt_spin_unlock(&(ksched_shared_info.uthread_init_lock));
	return 0;
}
extern void gt_yield()
{
    printf("\tI got removed\t ");

    uthread_schedule(&sched_find_best_uthread);
}
extern void uthread_schedule(uthread_struct_t * (*kthread_best_sched_uthread)(kthread_runqueue_t *))
{
	kthread_context_t *k_ctx;
	kthread_runqueue_t *kthread_runq;
	uthread_struct_t *u_obj;
    //struct timeval t1,t2;
     long int sub_cred;
     long int s,us,elapsed;
     int f;
	/* Signals used for cpu_thread scheduling */
	 //kthread_block_signal(SIGVTALRM);
	 //kthread_block_signal(SIGUSR1);

    //sub_cred = ((elapsed/100000)*25.0);


#if 0
	fprintf(stderr, "uthread_schedule invoked !!\n");
#endif
	k_ctx = kthread_cpu_map[kthread_apic_id()];
	kthread_runq = &(k_ctx->krunqueue);

	if((u_obj = kthread_runq->cur_uthread))
	{

        /*Go through the runq and schedule the next thread to run */
        gettimeofday(&(u_obj->t1),NULL);

        kthread_runq->cur_uthread = NULL;

        if(u_obj->uthread_state & (UTHREAD_DONE | UTHREAD_CANCELLED))
		{
			/* XXX: Inserting uthread into zombie queue is causing improper
			 * cleanup/exit of uthread (core dump) */
			uthread_head_t * kthread_zhead = &(kthread_runq->zombie_uthreads);
			gt_spin_lock(&(kthread_runq->kthread_runqlock));
			kthread_runq->kthread_runqlock.holder = 0x01;
			TAILQ_INSERT_TAIL(kthread_zhead, u_obj, uthread_runq);
			gt_spin_unlock(&(kthread_runq->kthread_runqlock));
		
			{
				ksched_shared_info_t *ksched_info = &ksched_shared_info;	
				gt_spin_lock(&ksched_info->ksched_lock);
				ksched_info->kthread_cur_uthreads--;
				gt_spin_unlock(&ksched_info->ksched_lock);
			}
		}
		else
		{
			/* XXX: Apply uthread_group_penalty before insertion */
			u_obj->uthread_state = UTHREAD_RUNNABLE;
            if (sched_fag==1)
            {

                elapsed =  (((u_obj->t1.tv_sec-u_obj->t2.tv_sec)*1000000L + (u_obj->t1.tv_usec-u_obj->t2.tv_usec))/100000.0L)*25;

                exe_time[u_obj->uthread_tid] = exe_time[u_obj->uthread_tid]+((u_obj->t1.tv_sec-u_obj->t2.tv_sec)*1000L + (u_obj->t1.tv_usec-u_obj->t2.tv_usec)/1000L) ;
                //if (u_obj->uthread_tid==50 && flage ==0){ printf("They kicked me out and burnt my credits %ld\n",elapsed); flage=1;}
                //printf("T1 is %ld\n",elapsed);
                //printf ("gfhgdfhgdfgchfh");
                //s = elapsed;
                //us = elapsed%1000000L;
                //sub_cred=25;
				//printf("My credits that are subtracted are %lu \n",elapsed);
                u_obj->credit=u_obj->credit-elapsed;
                    if(u_obj->credit<0)
                        add_to_runqueue(kthread_runq->expires_runq, &(kthread_runq->kthread_runqlock), u_obj);
                    if(u_obj->credit>0)
                        add_to_runqueue(kthread_runq->active_runq, &(kthread_runq->kthread_runqlock),u_obj);
            }
            if(!sched_fag) {add_to_runqueue(kthread_runq->expires_runq, &(kthread_runq->kthread_runqlock), u_obj);}
			/* XXX: Save the context (signal mask not saved) */
			if(sigsetjmp(u_obj->uthread_env, 0))
				return;
		}
	}

	/* kthread_best_sched_uthread acquires kthread_runqlock. Dont lock it up when calling the function. */
	if(!(u_obj = kthread_best_sched_uthread(kthread_runq)))
	{
		/* Done executing all uthreads. Return to main */
		/* XXX: We can actually get rid of KTHREAD_DONE flag */
		if(ksched_shared_info.kthread_tot_uthreads && !ksched_shared_info.kthread_cur_uthreads)
		{
			fprintf(stderr, "Quitting kthread (%d)\n", k_ctx->cpuid);
			k_ctx->kthread_flags |= KTHREAD_DONE;
		}

        siglongjmp(k_ctx->kthread_env, 1);
		return;
	}

	kthread_runq->cur_uthread = u_obj;
	if((u_obj->uthread_state == UTHREAD_INIT) && (uthread_init(u_obj)))
	{
		fprintf(stderr, "uthread_init failed on kthread(%d)\n", k_ctx->cpuid);
		exit(0);
	}

	u_obj->uthread_state = UTHREAD_RUNNING;
	
	/* Re-install the scheduling signal handlers */
	kthread_install_sighandler(SIGVTALRM, k_ctx->kthread_sched_timer);
	kthread_install_sighandler(SIGUSR1, k_ctx->kthread_sched_relay);
	/* Jump to the selected uthread context */
    gettimeofday(&(u_obj->t2),NULL);
    siglongjmp(u_obj->uthread_env, 1);

    return;
}


/* For uthreads, we obtain a seperate stack by registering an alternate
 * stack for SIGUSR2 signal. Once the context is saved, we turn this 
 * into a regular stack for uthread (by using SS_DISABLE). */
static void uthread_context_func(int signo)
{
    int f;
	uthread_struct_t *cur_uthread;
	kthread_runqueue_t *kthread_runq;

	kthread_runq = &(kthread_cpu_map[kthread_apic_id()]->krunqueue);

	printf("..... uthread_context_func .....\n");
	/* kthread->cur_uthread points to newly created uthread */
	if(!sigsetjmp(kthread_runq->cur_uthread->uthread_env,0))
	{
		/* In UTHREAD_INIT : saves the context and returns.
		 * Otherwise, continues execution. */
		/* DONT USE any locks here !! */
		assert(kthread_runq->cur_uthread->uthread_state == UTHREAD_INIT);
		kthread_runq->cur_uthread->uthread_state = UTHREAD_RUNNABLE;
		return;
	}

	/* UTHREAD_RUNNING : siglongjmp was executed. */

    cur_uthread = kthread_runq->cur_uthread;
	assert(cur_uthread->uthread_state == UTHREAD_RUNNING);
	/* Execute the uthread task */
	cur_uthread->uthread_func(cur_uthread->uthread_arg);
    if ((cur_uthread->uthread_tid == 50)&& flaga==0) {usleep(1000);flaga=1;}
	cur_uthread->uthread_state = UTHREAD_DONE;

    uthread_schedule(&sched_find_best_uthread);
	return;
}

/**********************************************************************/
/* uthread creation */

extern kthread_runqueue_t *ksched_find_target(uthread_struct_t *);

extern int uthread_create(uthread_t *u_tid, int (*u_func)(void *), void *u_arg, uthread_group_t u_gid ,int credit)
{
	kthread_runqueue_t *kthread_runq;
	uthread_struct_t *u_new;

	/* Signals used for cpu_thread scheduling */
	// kthread_block_signal(SIGVTALRM);
	// kthread_block_signal(SIGUSR1);

	/* create a new uthread structure and fill it */
	if(!(u_new = (uthread_struct_t *)MALLOCZ_SAFE(sizeof(uthread_struct_t))))
	{
		fprintf(stderr, "uthread mem alloc failure !!");
		exit(0);
	}

    if (sched_fag==1)u_new->credit=credit;
	u_new->uthread_state = UTHREAD_INIT;
	u_new->uthread_priority = DEFAULT_UTHREAD_PRIORITY;
	u_new->uthread_gid = u_gid;
	u_new->uthread_func = u_func;
	u_new->uthread_arg = u_arg;


	/* Allocate new stack for uthread */
	u_new->uthread_stack.ss_flags = 0; /* Stack enabled for signal handling */
	if(!(u_new->uthread_stack.ss_sp = (void *)MALLOC_SAFE(UTHREAD_DEFAULT_SSIZE)))
	{
		fprintf(stderr, "uthread stack mem alloc failure !!");
		return -1;
	}
	u_new->uthread_stack.ss_size = UTHREAD_DEFAULT_SSIZE;


	{
		ksched_shared_info_t *ksched_info = &ksched_shared_info;

		gt_spin_lock(&ksched_info->ksched_lock);
		u_new->uthread_tid = ksched_info->kthread_tot_uthreads++;
		ksched_info->kthread_cur_uthreads++;
		gt_spin_unlock(&ksched_info->ksched_lock);
	}

	/* XXX: ksched_find_target should be a function pointer */
	kthread_runq = ksched_find_target(u_new);

	*u_tid = u_new->uthread_tid;
	/* Queue the uthread for target-cpu. Let target-cpu take care of initialization. */
	add_to_runqueue(kthread_runq->active_runq, &(kthread_runq->kthread_runqlock), u_new);


	/* WARNING : DONOT USE u_new WITHOUT A LOCK, ONCE IT IS ENQUEUED. */

	/* Resume with the old thread (with all signals enabled) */
	// kthread_unblock_signal(SIGVTALRM);
	// kthread_unblock_signal(SIGUSR1);

	return 0;
}

#if 0
/**********************************************************************/
kthread_runqueue_t kthread_runqueue;
kthread_runqueue_t *kthread_runq = &kthread_runqueue;
sigjmp_buf kthread_env;

/* Main Test */
typedef struct uthread_arg
{
	int num1;
	int num2;
	int num3;
	int num4;	
} uthread_arg_t;

#define NUM_THREADS 10
static int func(void *arg);

int main()
{
	uthread_struct_t *uthread;
	uthread_t u_tid;
	uthread_arg_t *uarg;

	int inx;

	/* XXX: Put this lock in kthread_shared_info_t */
	gt_spinlock_init(&uthread_group_penalty_lock);

	/* spin locks are initialized internally */
	kthread_init_runqueue(kthread_runq);

	for(inx=0; inx<NUM_THREADS; inx++)
	{
		uarg = (uthread_arg_t *)MALLOC_SAFE(sizeof(uthread_arg_t));
		uarg->num1 = inx;
		uarg->num2 = 0x33;
		uarg->num3 = 0x55;
		uarg->num4 = 0x77;
		uthread_create(&u_tid, func, uarg, (inx % MAX_UTHREAD_GROUPS),0);
	}

	kthread_init_vtalrm_timeslice();
	kthread_install_sighandler(SIGVTALRM, kthread_sched_vtalrm_handler);
	if(sigsetjmp(kthread_env, 0) > 0)
	{
		/* XXX: (TODO) : uthread cleanup */
		exit(0);
	}
	
	uthread_schedule(&ksched_priority);
	return(0);
}

static int func(void *arg)
{
	unsigned int count;
#define u_info ((uthread_arg_t *)arg)
	printf("Thread %d created\n", u_info->num1);
	count = 0;
	while(count <= 0xffffff)
	{
		if(!(count % 5000000))
			printf("uthread(%d) => count : %d\n", u_info->num1, count);
		count++;
	}
#undef u_info
	return 0;
}
#endif
