#pragma once


class HeuristicAlgorithm
{
public:
	HeuristicAlgorithm(map<int, Node> mapNodes);
	~HeuristicAlgorithm(void);
	bool OffloadDecision(int source, int dest, double datasize, double deadline, double &cellprob, double &offloadprob, 
		bool bsimulate, double datapercen, int runs, double &cellularsimprob, double &offloadsimprob);

private:
	map<int, Node> mapNodesInside;
	bool PathPlanning(map<int, Node> & mapNodes, double datasize, double deadline, int source, int dest, double datapercent, 
		double &cellprob, double &offloadprob, bool bsimulate, int runs, double &cellularsimprob, double &offloadsimprob);
	
	//BOOL DijkstraPathCapacity(map<int, Node> & mapNodes, int source, int dest, double & pathcapacity, vector<int> & vecPath);
	//BOOL DijkstraContactProbability(map<int, Node> & mapNodes, int source, int dest, double deadline, double & probability, 
	//		double & pathcapacity, vector<int> & vecPath);
	BOOL DijkstraUploadingProbability(map<int, Node> & mapNodes, int source, int dest, double deadline, double threshold, 
		double & probability, double & pathcapacity, vector<int> & vecPath);
	
	//double PathCapacityPreferred(map<int, Node> & mapNodes, vector<vector<int>> & vecPaths, vector<double> & vecCapacities,
	//		vector<double> & vecProbs, double deadline, int source, int dest);
	//double ContactProbabilityPreferred(map<int, Node> & mapNodes, vector<vector<int>> & vecPaths, vector<double> & vecCapacities, 
	//		vector<double> & vecProbs, double datasize, double deadline, int source, int dest);
	double JointConsiderationPreferred(map<int, Node> & mapNodes, vector<vector<int>> & vecPaths, vector<double> & vecCapacities, 
			vector<double> & vecProbs, double datasize, double deadline, int source, int dest, double thresholdprob, bool bProportional = false);

	
	//double CalculateContactProbability(map<int, Node> & mapNodes, vector<int> vecPath, double deadline);
	double ApproximateContatcProbability(map<int, Node> & mapNodes, vector<int> vecPath, vector<int> vecContacts, double deadline, double dataassigment=0);
	double ApproximateContatcComplementProb(map<int, Node> & mapNodes, vector<int> vecPath, vector<int> vecContacts, int indication, double deadline, double dataassigment);
	double ApproximateUploadProbability(map<int, Node> & mapNodes, vector<int> vecPath, vector<int> vecContacts, double deadline, double dataassigment);

	

	void AllocateResourceInclusive(map<int,Node> & mapNodes, vector<int> vecPath, double capacity);
	void AllocateResourceExclusive(map<int,Node> & mapNodes, vector<int> vecPath);
	void ReleaseResourceInclusive(map<int,Node> & mapNodes, vector<int> vecPath, double capacity);
	void ReleaseResourceExclusive(map<int,Node> & mapNodes, vector<int> vecPath);

	//double AssignData2Paths(map<int, Node> & mapNodes, double dataallocated, vector<vector<int>> & vecPaths, 
	//		vector<double> & vecContactProbs, vector<double> & vecCapacities, vector<double> & vecProbs); 

	double ProcessiveAssignData2Paths(map<int, Node> & mapNodes, double dataallocated, double deadline, vector<vector<int>> & vecPaths, 
			vector<double> & vecCapacities, vector<double> & vecProbs);
	double ProportionalAssignData2Paths(map<int, Node> & mapNodes, double dataallocated, double deadline, vector<vector<int>> & vecPaths, 
			vector<double> & vecMetrics, vector<double> & vecCapacities, vector<double> & vecProbs);

	double TransSimulationIndividual(map<int, Node> & mapNodes, int source, int dest, 
		double datasize, double datapercent, double deadline, int number=1);
	double TransSimulationCooperative(map<int, Node> & mapNodes, vector<vector<int>> & vecPaths, vector<double> & vecCapacities, 
		int source, int dest, double datasize, double datapercent, double deadline, int number=1);

	
	//	void RecursiveCombination(vector<int> input, vector<int>& temp, vector<vector<int>> & output);
};

