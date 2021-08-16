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
    void simulate_(int *, int *, int*, int*);
    int getmax_(int*);
    void det2img_(int*, int*);
    void printparams_();
}

void timeoutCallback(int sig)
{
	int workerId;
	MPI_Comm_rank(MPI_COMM_WORLD, &workerId);
	processcommon_.continueexecution = 0;
	printf("Process [%d]: Timer signal received. Shutting down process..\n", workerId);
}

int main(int argc, char* argv[])
{

     printf("INICIO ****\n");
	if( ! parseConfigFile("parameters.xml") )
	{
		std::cout<<"Error parsing parameters.xml file"<<std::endl;
		return -1;
	}	
	int ierr, workerCount, workerId, maxValue;
	time_t startTime, endTime;

	ierr = MPI_Init(NULL, NULL);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &workerId);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &workerCount);

	std::cout<<"MPI rank:["<<workerId<<"]"<<std::endl;

	if( workerId == 0 )
	{
		std::cout<<"Number of projections: ["<<projectionCount<<"]"<<std::endl;
		std::cout<<"Number of processes: ["<<workerCount<<"]"<<std::endl;
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

	int npmax = (int)beamcommon_.npmax;
	int photonCount = npmax / workerCount;
	if (workerId == 0) {
		photonCount = npmax - photonCount * (workerCount - 1);
	}
			
	for(int proyNum=1; proyNum <= projectionCount && processcommon_.continueexecution == 1; proyNum++)
	{
		std::cout<<"Process: ["<<workerId<<"] Proy: ["<<proyNum<<"]"<<std::endl;
		simulate_(&projectionCount, &proyNum, &workerId, &photonCount);
		MPI_Barrier(MPI_COMM_WORLD);
	} 


	MPI_Barrier(MPI_COMM_WORLD);

	if( workerId == 0 && processcommon_.continueexecution == 1 )
	{
		//maxValue = getmax_(&projectionCount);
		//std::cout<<"Global MaxValue: ["<<maxValue<<"]"<<std::endl;

		//for(int proyNum=1; proyNum <= projectionCount; proyNum++)
		//{
		//	det2img_(&maxValue, &proyNum);
		//}
	
		time(&endTime);
		std::cout<<"Total execution time: ["<<difftime(endTime, startTime)<<"] seconds."<<std::endl; 
	}

	MPI_Finalize();
	return 0;

}
