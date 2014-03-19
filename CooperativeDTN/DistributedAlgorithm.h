#pragma once




class DistributedAlgorithm
{
public:
	DistributedAlgorithm(map<int, Node> Nodes);
	~DistributedAlgorithm(void);

	map<int, DistNode> mapDistNodes;
	
	void InitDataTransfer(int nSrc, int nDest, double datasize, double expiretime, double timeleft);
	void ExecuteNodeContact(int fromnode, int tonode, double contacttime, double contactduration);
	bool IsDataRecvAtDest(double datasize, double dbprecent);

private:
	map<int, Node> mapNodes;
	double dbRecvDataAtDest;
	
	void InitNode(int nID, int nSrc, int nDest, double datasize, double dataremain, double expiretime, double timeleft);
	void Find2HopPaths(int source, int dest, vector<vector<int>> & vecPaths);
	double CalculatePathContactProb(vector<int> & vecPath, double deadline);
	double CalculatePathCapacity(vector<int> & vecPath);
	double CalculatePathUploadProb(vector<int> & vecPath, double assigneddata, double deadline);
	void UpdateUploadProb(DistNode & node);
	double ReallocateData(DistNode & nodefrom, DistNode & nodeto);
	double ApproximateUploadProbability(vector<int> vecPath, vector<int> vecContacts, double deadline, double dataassigment=0.0);
	void CriterionAssignment(double datasize, double deadline, vector<vector<int>> & vecPaths, vector<double> &vecDataAssign);
	void AfterDataTransAtSender(DistNode & beforenode, double sentdata);
	void AfterDataTransAtReceiver(DistNode & afternode, double misseddata);
	
};

