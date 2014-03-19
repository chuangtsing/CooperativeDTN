#include "stdafx.h"
#include "DistributedAlgorithm.h"
#include "mathfunctions.h"
#include "Alglib/specialfunctions.h"


using namespace alglib;


extern void RecursiveCombination(vector<int> input, vector<int>& temp, vector<vector<int>> & output);


DistributedAlgorithm::DistributedAlgorithm(map<int, Node> Nodes)
{
	mapNodes = Nodes;
	dbRecvDataAtDest = 0.0;
}


DistributedAlgorithm::~DistributedAlgorithm(void)
{

}


void DistributedAlgorithm::InitDataTransfer(int nSrc, int nDest, double datasize, double expiretime, double timeleft)
{
	InitNode(nSrc, nSrc, nDest, datasize, datasize, expiretime, timeleft);
}


void DistributedAlgorithm::InitNode(int nID, int nSrc, int nDest, double datasize, double dataremain, double expiretime, double timeleft)
{
	DistNode node;
	
	node.nID = nID;
	node.nSource = nSrc;
	node.nDest = nDest;
	node.dbDataSize = datasize;
	node.dbTimeExpired = expiretime;
	node.dbTimeLeft = timeleft;
	node.dbDataRemain = dataremain;

	if(nSrc == nID)
	{
		Find2HopPaths(nSrc, nDest, node.vecPaths);
		CriterionAssignment(node.dbDataSize, node.dbTimeLeft, node.vecPaths, node.vecDataAssginments);
	}

	mapDistNodes.insert(pair<int, DistNode>(nID, node));
}


void DistributedAlgorithm::Find2HopPaths(int from, int to, vector<vector<int>> & vecPaths)
{
	vecPaths.clear();

	Node src = mapNodes[from];

	for(map<int, Parameters>::iterator i= src.NeighborContactInfo.begin(); 
		i!=src.NeighborContactInfo.end(); i++)
	{
		if(i->first == to)
		{
			vector<int> path;
			path.push_back(from);
			path.push_back(to);
			vecPaths.push_back(path);
		}
		else
		{
			
			for(map<int, Parameters>::iterator j=mapNodes[i->first].NeighborContactInfo.begin();
				j!=mapNodes[i->first].NeighborContactInfo.end(); j++)
			{
				if(j->first == to)
				{
					vector<int> path;
					path.push_back(from);
					path.push_back(i->first);
					path.push_back(to);
					vecPaths.push_back(path);
				}
			}
		}
	}
}


double DistributedAlgorithm::CalculatePathContactProb(vector<int> & vecPath, double deadline)
{
	if(vecPath.size() <= 1) return 0.0;
	int K = vecPath.size()-1;
	vector<double> lambda;

	for(int i=0, j=1; i<K; i++,j++)
	{
		lambda.push_back(mapNodes[vecPath[i]].NeighborContactInfo[vecPath[j]].lambda);
	}

	double prob = 0.0;

	if(K==1)
	{
		prob = 1-exp(-lambda[0] * deadline);
	}
	else
	{
		for(int i=0; i<K; i++)
		{
			double part1 = 1-exp(-deadline * lambda[i]);
			double part2 = 1.0;
			for(int j=0; j<K; j++)
			{
				if(j==i) continue;
				part2 *= lambda[j]/(lambda[j]-lambda[i]);
			}
			prob += part1*part2;
		}
	}
	return prob;
}


double DistributedAlgorithm::CalculatePathCapacity(vector<int> & vecPath)
{
	double capacity = -1.0;

	for(int i=0; i<vecPath.size()-1; i++)
	{
		if(capacity > mapNodes[vecPath[i]].NeighborContactInfo[vecPath[i+1]].beta
			|| capacity < 0)
			capacity = mapNodes[vecPath[i]].NeighborContactInfo[vecPath[i+1]].beta;
	}

	return capacity;
}


double DistributedAlgorithm::CalculatePathUploadProb(vector<int> & vecPath, double assigneddata, double deadline)
{
	if(assigneddata == 0.0) return 1.;

	//double prob = CalculatePathContactProb(vecPath, deadline);

	//for(size_t i=0; i<vecPath.size()-1; i++)
	//{
	//	double alpha = mapNodes[vecPath[i]].NeighborContactInfo[vecPath[i+1]].alpha;
	//	double beta = mapNodes[vecPath[i]].NeighborContactInfo[vecPath[i+1]].beta;
	//	if(assigneddata > beta)
	//		prob *= pow(beta/assigneddata, alpha); 
	//}
	//return prob;


	double contactprob=0.0; 
	vector<int> vecContacts, vecTemp;
	vector<vector<int>> vecCombinations;

	for(int i=0; i<vecPath.size()-1; i++)
	{
		double beta = mapNodes[vecPath[i]].NeighborContactInfo[vecPath[i+1]].beta;
		vecContacts.push_back((int)ceil(assigneddata/beta));
	}

	RecursiveCombination(vecContacts, vecTemp, vecCombinations);

	for(int i=0; i<vecCombinations.size(); i++)
	{
		double dbprob = ApproximateUploadProbability(vecPath, vecCombinations[i], deadline, assigneddata);
		if(dbprob > contactprob) contactprob = dbprob;
	}

	return contactprob;
}


void DistributedAlgorithm::UpdateUploadProb(DistNode & node)
{
	for(size_t i=0; i<node.vecPaths.size(); i++)
	{
		node.vecDataUploadProbs[i] = CalculatePathUploadProb(node.vecPaths[i], 
			node.vecDataAssginments[i], node.dbTimeLeft);
	}
}


double DistributedAlgorithm::ApproximateUploadProbability(vector<int> vecPath, vector<int> vecContacts, double deadline, double dataassigment)
{
	if(vecPath.size() - vecContacts.size() != 1) return 0.0;
	if(vecPath.size() <= 1) return 0.0;
	int K = vecPath.size()-1;
	vector<double> alpha, beta, lambda;
	
	for(int i=0, j=1; i<K; i++,j++)
	{
		alpha.push_back(mapNodes[vecPath[i]].NeighborContactInfo[vecPath[j]].alpha);
		beta.push_back(mapNodes[vecPath[i]].NeighborContactInfo[vecPath[j]].beta);
		lambda.push_back(mapNodes[vecPath[i]].NeighborContactInfo[vecPath[j]].lambda);
	}

	
	double contactprob, dataprob=1.;
	if(K==1)
	{
		contactprob = incompletegamma(vecContacts[0], lambda[0]*deadline);
		if(dataassigment > beta[0])
			dataprob = pow(beta[0]/dataassigment, alpha[0]);
		return contactprob*dataprob;
	}

	double mu=0.0, sqr = 0.0;

	for(int i=0; i<K; i++)
	{
		mu += vecContacts[i]/lambda[i];
		sqr += vecContacts[i]/(lambda[i]*lambda[i]);
		if(dataassigment > beta[i])
			dataprob *= 1.0-pow(1.0-pow(beta[i]/dataassigment, alpha[i]),vecContacts[i]);
	}

	double p, theta;

	p = mu*mu/sqr;
	theta = sqr/mu;

	contactprob = incompletegamma(p, deadline/theta);
		
	return contactprob * dataprob;
}


void DistributedAlgorithm::CriterionAssignment(double datasize, double deadline, vector<vector<int>> & vecPaths, vector<double> &vecDataAssign)
{
	vecDataAssign.clear();

	vector<double> vecMetrics;
	double totalMetrics = 0.;
	for(size_t i=0; i<vecPaths.size(); i++)
	{
		double prob = CalculatePathContactProb(vecPaths[i], deadline);
		double capa = CalculatePathCapacity(vecPaths[i]);
		vecMetrics.push_back(prob*capa);
		totalMetrics += prob*capa;
	}

	double alreadyassgined = 0.;
	for(size_t i=0; i<vecPaths.size(); i++)
	{
		if(i = vecPaths.size()-1)
		{
			vecDataAssign.push_back(datasize - alreadyassgined);
		}
		else
		{
			double assigned = datasize * vecMetrics[i]/totalMetrics;
			alreadyassgined += assigned;
			vecDataAssign.push_back(assigned);
		}
	}
}


void DistributedAlgorithm::AfterDataTransAtSender(DistNode & beforenode, double sentdata)
{
	beforenode.dbDataRemain = beforenode.dbDataRemain - sentdata <= 0.0 ? 0.0 : beforenode.dbDataRemain - sentdata;
	if(beforenode.dbDataRemain == 0.0)
	{
		for(size_t i=0; i<beforenode.vecPaths.size(); i++)
		{
			beforenode.vecDataAssginments[i] = 0.0;  
		}
	}
	else
	{
		//same as criterion assignment
		//CriterionAssignment(node.dbDataRemain, node.dbTimeLeft, node.vecPaths, node.vecDataAssginments);
		//deduct the transfered data from current assginment
		while(sentdata != 0.0)
		{
			double minprob = 1.0;
			int index = -1;
			for(size_t i=0; i<beforenode.vecPaths.size(); i++)
			{
				if(beforenode.vecDataUploadProbs[i] < minprob)
				{
					minprob = beforenode.vecDataUploadProbs[i];
					index = -1;
				}
			}
			if(index == -1) break;
			if(sentdata > beforenode.vecDataAssginments[index])
			{
				sentdata -= beforenode.vecDataAssginments[index];
				beforenode.vecDataAssginments[index] = 0;
			}
			else
			{
				beforenode.vecDataAssginments[index] -= sentdata;
				sentdata = 0.0;
			}			
		}
	}
}


void DistributedAlgorithm::AfterDataTransAtReceiver(DistNode & afternode, double misseddata)
{
	while(misseddata != 0.0)
	{
		double minprob = 1.0;
		int index = -1;
		for(size_t i=0; i<afternode.vecPaths.size(); i++)
		{
			if(afternode.vecDataUploadProbs[i] < minprob)
			{
				minprob = afternode.vecDataUploadProbs[i];
				index = -1;
			}
		}
		if(index == -1) break;
		if(misseddata > afternode.vecDataAssginments[index])
		{
			misseddata -= afternode.vecDataAssginments[index];
			afternode.vecDataAssginments[index] = 0;
		}
		else
		{
			afternode.vecDataAssginments[index] -= misseddata;
			misseddata = 0.0;
		}			
	}

	double dbtotaldata = 0.0;

	for(size_t i=0; i<afternode.vecDataAssginments.size(); i++)
	{
		dbtotaldata += afternode.vecDataAssginments[i];
	}

	afternode.dbDataRemain = dbtotaldata;
}


double DistributedAlgorithm::ReallocateData(DistNode & nodefrom, DistNode & nodeto)
{
	double totalreallocatedData = 0.0;

	double dataAtnNoData = 0.0;
	for (size_t i=0; i<nodefrom.vecPaths.size(); i++)
	{
		if(nodefrom.vecPaths[i][1] == nodeto.nID)
		{
			dataAtnNoData = nodefrom.vecDataAssginments[i];
			nodefrom.vecDataAssginments[i] = 0.0;
			break;
		}
	}

	for (size_t i=0; i<nodeto.vecPaths.size(); i++)
	{
		if(nodeto.vecPaths[i][1] == nodeto.nDest)
		{
			nodeto.vecDataAssginments[i] = dataAtnNoData;
			break;
		}
	}
	totalreallocatedData += dataAtnNoData;

	UpdateUploadProb(nodefrom);
	UpdateUploadProb(nodeto);

	while(1)
	{
		double minprob = 1.;
		double reallocatedData = 0.;
		int candidate = -1;
		for(size_t i=0; i<nodefrom.vecDataUploadProbs.size(); i++)
		{
			if(minprob > nodefrom.vecDataUploadProbs[i])
			{
				minprob = nodefrom.vecDataUploadProbs[i];
				candidate = i;
				reallocatedData = nodefrom.vecDataAssginments[i];
			}
		}

		if(candidate = -1) break;

		bool bFirstFlag = true;
		double combinedprob; 
		int curpath = -1;
		double curprob;
		for(size_t i=0; i<nodeto.vecPaths.size(); i++)
		{
			double preprob = CalculatePathUploadProb(nodeto.vecPaths[i], nodeto.vecDataAssginments[i],nodeto.dbTimeLeft);
			double combined = nodeto.vecDataAssginments[i] + reallocatedData;
			double afterprob = CalculatePathUploadProb(nodeto.vecPaths[i], combined, nodeto.dbTimeLeft);
			if(preprob * minprob < afterprob)
			{
				if(bFirstFlag == true)
				{
					bFirstFlag = false;
					combinedprob = afterprob;
					curprob = preprob;
					curpath = i;
				}
				else
				{
					if(combinedprob * preprob < curprob * afterprob)
					{
						combinedprob = afterprob;
						curprob = preprob;
						curpath = i;
					}
				}
			}
		}
		if(curpath == -1) break;
		nodefrom.vecDataAssginments[candidate] = 0.0;
		nodeto.vecDataAssginments[curpath] += reallocatedData;
		totalreallocatedData += reallocatedData;
	}

	return totalreallocatedData;
}


void DistributedAlgorithm::ExecuteNodeContact(int fromnode, int tonode, double contacttime, double contactduration)
{
	double maxdata = 0.;

	bool bFrom = mapDistNodes.count(fromnode);
	bool bTo = mapDistNodes.count(tonode);

	//neithor has data to upload
	if(!bFrom && !bTo) return;

	//one of them has the data
	if((!bFrom && bTo) || (bFrom && !bTo))
	{
		int nHasData = bFrom == true ? fromnode : tonode;
		int nNoData = bFrom == true ? tonode : fromnode;
		DistNode nodewithdata = mapDistNodes[nHasData];
		nodewithdata.dbTimeLeft = nodewithdata.dbTimeExpired - contacttime;

		if(nNoData == nodewithdata.nDest)
		{
			if(contactduration >= nodewithdata.dbDataRemain)
			{
				dbRecvDataAtDest += nodewithdata.dbDataRemain;
				AfterDataTransAtSender(nodewithdata, nodewithdata.dbDataRemain);
			}
			else
			{
				dbRecvDataAtDest += contactduration;
				AfterDataTransAtSender(nodewithdata, contactduration);
			}
			mapDistNodes[nHasData] = nodewithdata;
			return;
		}

		InitNode(nNoData, nodewithdata.nSource, nodewithdata.nDest, nodewithdata.dbDataSize, 0.0, 
			nodewithdata.dbTimeExpired, nodewithdata.dbTimeLeft);
		DistNode nodenodata = mapDistNodes[nNoData];

		//find the 2-hop paths from the encountered node to the dest;
		//if found no paths, return
		Find2HopPaths(nNoData, nodenodata.nDest, nodenodata.vecPaths);
		if(nodenodata.vecPaths.size() == 0) return;
		//exclude the path from current node to the dest;
		for (size_t i = nodenodata.vecPaths.size()-1; i<0; i--)
		{
			if(nodenodata.vecPaths[i][1] == nodewithdata.nID || 
				nodenodata.vecPaths[i][1] == nodewithdata.nSource)
				nodenodata.vecPaths.erase(nodenodata.vecPaths.begin()+i);
		}
		
		//assign data at nNoData;
		for (size_t i=0; i<nodenodata.vecPaths.size(); i++)
		{
			nodenodata.vecDataAssginments.push_back(0.0);
		}

		DistNode nodebefore = nodewithdata;
		maxdata = ReallocateData(nodewithdata, nodenodata);


		AfterDataTransAtSender(nodebefore, contactduration > maxdata ? maxdata : contactduration);
		AfterDataTransAtReceiver(nodenodata, contactduration > maxdata ? 0 : maxdata-contactduration);

		if(maxdata > 0.0 && contactduration > 0.0)
		{
			nodenodata.setPreNodes.insert(nodebefore.nID);
		}

		mapDistNodes[nHasData] = nodebefore;
		mapDistNodes[nNoData] = nodenodata;

		return;
	}

	//both of them have the data
	if(bFrom && bTo)
	{
		bool fromto = mapDistNodes[tonode].setPreNodes.count(fromnode);
		bool tofrom = mapDistNodes[fromnode].setPreNodes.count(tonode);

		if(fromto && tofrom) return;
		
		if(fromto || tofrom)
		{
			int nHasData = fromto == true ? fromnode : tonode;
			int nNoData = fromto == true ? tonode : fromnode;
			DistNode nodewithdata = mapDistNodes[nHasData];
			DistNode nodenodata = mapDistNodes[nNoData];
			nodewithdata.dbTimeLeft = nodewithdata.dbTimeExpired - contacttime;
			nodenodata.dbTimeLeft = nodenodata.dbTimeExpired - contacttime;
			DistNode nodebefore = nodewithdata;
			maxdata = ReallocateData(nodewithdata, nodenodata);
			AfterDataTransAtSender(nodebefore, contactduration > maxdata ? maxdata : contactduration);
			AfterDataTransAtReceiver(nodenodata, contactduration > maxdata ? 0 : maxdata-contactduration);

			if(maxdata > 0.0 && contactduration > 0.0)
			{
				nodenodata.setPreNodes.insert(nodebefore.nID);
			}

			mapDistNodes[nHasData] = nodebefore;
			mapDistNodes[nNoData] = nodenodata;

			return;
		}

		if(!fromto && !tofrom)
		{
			DistNode nodefrom = mapDistNodes[fromnode];
			DistNode nodeto = mapDistNodes[tonode];
			nodefrom.dbTimeLeft = nodefrom.dbTimeExpired - contacttime;
			nodeto.dbTimeLeft = nodeto.dbTimeExpired - contacttime;

			DistNode nodefromcopy = nodefrom;
			DistNode nodetocopy = nodeto;
			double dbfromto = ReallocateData(nodefromcopy, nodetocopy);
			nodefromcopy = nodefrom;
			nodetocopy = nodeto;
			double dbtofrom = ReallocateData(nodeto, nodefrom);

			if(dbfromto > dbtofrom)
			{
				DistNode nodebefore = nodefrom;
				maxdata = ReallocateData(nodefrom, nodeto);
				AfterDataTransAtSender(nodebefore, contactduration > maxdata ? maxdata : contactduration);
				AfterDataTransAtReceiver(nodeto, contactduration > maxdata ? 0 : maxdata-contactduration);

				if(maxdata > 0.0 && contactduration > 0.0)
				{
					nodeto.setPreNodes.insert(nodebefore.nID);
				}

				mapDistNodes[fromnode] = nodebefore;
				mapDistNodes[tonode] = nodeto;
			}
			else
			{
				DistNode nodebefore = nodeto;
				maxdata = ReallocateData(nodeto, nodefrom);
				AfterDataTransAtSender(nodebefore, contactduration > maxdata ? maxdata : contactduration);
				AfterDataTransAtReceiver(nodefrom, contactduration > maxdata ? 0 : maxdata-contactduration);

				if(maxdata > 0.0 && contactduration > 0.0)
				{
					nodefrom.setPreNodes.insert(nodebefore.nID);
				}

				mapDistNodes[fromnode] = nodefrom;
				mapDistNodes[tonode] = nodebefore;
			}
			return;
		}
	}
}


bool DistributedAlgorithm::IsDataRecvAtDest(double datasize, double dbprecent)
{
	if(datasize * dbprecent <= dbRecvDataAtDest)
		return true;
	return false;
}