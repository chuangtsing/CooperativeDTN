// CooperativeDTN.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "ProcessMIT.h"
#include "HeuristicAlgorithm.h"


extern void GenerateSyntheticNetwork(string inputfile, string outputfile, pair<double,double> dtnalpha, 
									 pair<double,double> dtnbeta, pair<double, double> dtnlambda, 
									 pair<double,double> cellularalpha, pair<double,double> cellularbeta,
									 pair<double,double> cellularlambda);

int _tmain(int argc, _TCHAR* argv[])
{

	GenerateSyntheticNetwork("SyntheticNetwork.dat", "100Nodes.dat", pair<double,double>(2,3), 
		pair<double,double>(6, 10), pair<double,double>(0.01, 0.1), pair<double,double>(2,3), 
		pair<double,double>(3, 4), pair<double,double>(0.001, 0.01));
	ProcessMIT MIT;
	//MIT.ExportNetwork();
	//MIT.StructNodesFromEdges();
	MIT.StructNodesFromFile("100Nodes.dat");
	
	FILE * fresult, *fstatresult;
	string strResult = "100NodesResult.dat";
	string strStatResult = "100NodesStatResult.dat";
	fresult = fopen(strResult.c_str(), "w+");
	if(fresult == NULL) return -1;
	fstatresult = fopen(strStatResult.c_str(), "w+");
	if(fstatresult == NULL) 
	{
		fclose(fresult);
		return -1;
	}

	
	//transmission parameters
	int datasize = 100, deadline = 400;
	int datastep = 1, deadlinestep = 400;
	int datainit = 1, deadlineinit = 400; 
	//simulation parameters
	double datapercent = 1;
	int runs = 1000;


	vector<vector<int>> vecData(datasize/datastep,vector<int>(deadline/deadlinestep, 0));  

	for(int i=2; i<3/*MIT.mapNodes.size()*/; i++) //Cellular is marked as the node with id=number of nodes
	{
		int helpedtran = 0;

		for(int j=datainit; j<=datasize; j=j+datastep)
		{
			for(int k=deadlineinit; k<=deadline; k=k+deadlinestep)
			{
				
				double cellularprob;
				double offloadprob;
				double cellularsimprob;
				double offloadsimprob;
				HeuristicAlgorithm heuristic(MIT.mapNodes);
				bool bOffload = heuristic.OffloadDecision(i, MIT.mapNodes.size(), j, k, cellularprob, offloadprob,
															true,datapercent,runs,cellularsimprob,offloadsimprob);
				if(bOffload) helpedtran++;
				fprintf(fresult, "%d,  %d,  %d,  %.8lf,  %.8lf,  %d,  %.8lf,  %.8lf\n", i, j, k, cellularprob, 
						offloadprob, bOffload, cellularsimprob, offloadsimprob);
				fflush(fresult);
				printf("%d  %d  %d  %d  Celluar %.8lf  Offload %.8lf  Cellular Sim %.8lf  Offload Sim %.8lf\n", 
					i, j, k, bOffload, cellularprob, offloadprob, cellularsimprob, offloadsimprob);
				fflush(stdout);

				if(bOffload)
					vecData[j/datastep-1][k/deadlinestep-1]+=1;
			}
		}
		fprintf(fstatresult, "%.8lf\n", (double)helpedtran/(datasize*deadline/datastep/deadlinestep));
		fflush(fstatresult);
	}

	for(int i=0; i<vecData.size(); i++)
	{
		for(int j=0; j<vecData[i].size(); j++)
			fprintf(fstatresult, "%d ", vecData[i][j]);
		fprintf(fstatresult, "\n");
	}

	fclose(fresult);
	fclose(fstatresult);

	return 0;
}

