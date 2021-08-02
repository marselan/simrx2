#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <signal.h>
#include "configparser.h"

extern "C"
{
    void simulate_(int *, int *);
    int getmax_(int*);
    void det2img_(int*, int*);
    void printparams_();
}

void timeoutCallback(int sig)
{
	int my_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	processcommon_.continueexecution = 0;
	printf("Process [%d]: Timer signal received. Shutting down process..\n", my_id);
}

int main(int argc, char* argv[])
{

     printf("INICIO ****\n");
	if( ! parseConfigFile("parameters.xml") )
	{
		std::cout<<"Error parsing parameters.xml file"<<std::endl;
		return -1;
	}	
	int ierr, num_procs, my_id, maxValue;
	time_t startTime, endTime;

	ierr = MPI_Init(NULL, NULL);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	std::cout<<"MPI rank:["<<my_id<<"]"<<std::endl;

	if( my_id == 0 )
	{
		std::cout<<"Number of projections: ["<<projectionCount<<"]"<<std::endl;
		std::cout<<"Number of processes: ["<<num_procs<<"]"<<std::endl;
		std::cout<<"Shutdown simrx after: ["<<processcommon_.shutdownafter<<"] seconds"<<std::endl;
		printparams_();
		time(&startTime);
	}

	processcommon_.continueexecution = 1;

	MPI_Barrier(MPI_COMM_WORLD);

	// set timer
	struct sigaction sact;
	struct itimerval timerval;

	sigemptyset( &sact.sa_mask );
	sact.sa_flags = 0;
	sact.sa_handler = timeoutCallback;
	sigaction( SIGALRM, &sact, NULL );

	timerval.it_interval.tv_sec = 0;
	timerval.it_interval.tv_usec = 0;
	timerval.it_value.tv_sec = processcommon_.shutdownafter;
	timerval.it_value.tv_usec = 0;

	setitimer(ITIMER_REAL, &timerval, NULL);

	// let's rock
	
	int procNum;
			
	for(int proyNum=1; proyNum <= projectionCount && processcommon_.continueexecution == 1; proyNum++)
	{
		procNum = proyNum % num_procs;
		if(procNum == my_id)
		{
			std::cout<<"Process: ["<<my_id<<"] Proy: ["<<proyNum<<"]"<<std::endl;
			simulate_(&proyNum, &projectionCount);
		}			
	} 


	MPI_Barrier(MPI_COMM_WORLD);

	if( my_id == 0 && processcommon_.continueexecution == 1 )
	{
		maxValue = getmax_(&projectionCount);
		std::cout<<"Global MaxValue: ["<<maxValue<<"]"<<std::endl;

		for(int proyNum=1; proyNum <= projectionCount; proyNum++)
		{
			det2img_(&maxValue, &proyNum);
		}
	
		time(&endTime);
		std::cout<<"Total execution time: ["<<difftime(endTime, startTime)<<"] seconds."<<std::endl; 
	}

	MPI_Finalize();
	return 0;

}
