#include "stdafx.h"
#include "HeuristicAlgorithm.h"
#include "mathfunctions.h"
#include "Alglib/specialfunctions.h"


using namespace alglib;


void RecursiveCombination(vector<int> input, vector<int>& temp, vector<vector<int>> & output)
{
	for(int i=1; i<=input[0]; i++)
	{
		if(input.size() > 1)
		{
			temp.push_back(i);
			vector<int> subinput(input);
			subinput.erase(subinput.begin());
			RecursiveCombination(subinput, temp, output);
		}
		else
		{
			vector<int> temptemp(temp);
			temptemp.push_back(i);
			output.push_back(temptemp);
		}
	}
	temp.clear();
}


HeuristicAlgorithm::HeuristicAlgorithm(map<int, Node> mapNodes, double datasize, double deadline, int source, int dest)
{
	mapNodesInside = mapNodes;
	PathPlanning(mapNodesInside, datasize, deadline, source, dest);
}


HeuristicAlgorithm::~HeuristicAlgorithm(void)
{


}


bool HeuristicAlgorithm::PathPlanning(map<int, Node> & mapNodes, double datasize, double deadline, int source, int dest)
{
	double pathcapacity, totalcapacity=0.;
	vector<int> vecPath;
	vector<vector<int>> vecPaths;
	vector<double> vecCapacities;
	vector<double> vecProbs;
	BOOL bEnoughCapacity = FALSE;


	double cellularprob=0.0;
	double offloadprob=0.0;
	double thresholdprob = 0.0;

	//Calculate the prob of uploading by itself

	if(mapNodes[source].NeighborContactInfo.count(dest) != 0)
	{
		double alpha = mapNodes[source].NeighborContactInfo[dest].alpha;
		double beta = mapNodes[source].NeighborContactInfo[dest].beta;
		double lambda = mapNodes[source].NeighborContactInfo[dest].lambda;
		double k = ceil(datasize/beta);

		for(int i=1; i<=k; i++)
		{
			//calculate the contact prob according to gamma distribiton
			double grammaprob = incompletegamma(i,lambda*deadline);
			double paretoprob = 1.0-pow(1.0-pow(beta/datasize, alpha),i);
			cellularprob = paretoprob * grammaprob > cellularprob ? paretoprob * grammaprob : cellularprob;
		}

		thresholdprob = incompletegamma(1,lambda*deadline)*(alpha*beta/(alpha-1));
	}

	
	// calculate the prob of offloading

	offloadprob = JointConsiderationPreferred(mapNodes,vecPaths,vecCapacities, vecProbs, datasize, deadline, source, dest, thresholdprob, true);
	return offloadprob > cellularprob;

	//while(DijkstraPathCapacity(mapNodes, source, dest, pathcapacity, vecPath))
	//{
	//	//Allocated resource of path
	//	AllocateResourceInclusive(mapNodes,vecPath, pathcapacity);

	//	totalcapacity += pathcapacity;
	//	vecCapacities.push_back(pathcapacity);
	//	vecPaths.push_back(vecPath);
	//	if(totalcapacity >= datasize) 
	//	{
	//		bEnoughCapacity = TRUE;
	//		break;
	//	}
	//}

	//if(bEnoughCapacity)
	//{
	//	//Path capacity preferred algorithm
	//	offloadprob = PathCapacityPreferred(mapNodes, vecPaths, vecCapacities, vecProbs, deadline, source, dest);
	//}
	//else
	//{
	//	//Release resource
	//	for(int i=0; i<vecPaths.size(); i++)
	//	{
	//		ReleaseResourceInclusive(mapNodes, vecPaths[i], vecCapacities[i]);
	//	}
	//	vecPaths.clear();
	//	vecCapacities.clear();

	//	//Contact probability preferred algorithm
	//	offloadprob = ContactProbabilityPreferred(mapNodes, vecPaths, vecCapacities, vecProbs, datasize, deadline, source, dest);
	//}

	
}


BOOL HeuristicAlgorithm::DijkstraPathCapacity(map<int, Node> & mapNodes, int source, int dest, double & pathcapacity, vector<int> & vecPath)
{
	pathcapacity = 0.0;
	vecPath.clear();

	map<int, double> mapUnvisitedNode;
	map<int, int> mapPreNode;

	for(map<int,Node>::iterator i = mapNodes.begin(); i != mapNodes.end(); i++)
	{
		if(i->first != source)
			mapUnvisitedNode.insert(pair<int,double>(i->first, -1.0));
		else
			mapUnvisitedNode.insert(pair<int,double>(i->first, 0.0));
		mapPreNode.insert(pair<int,double>(i->first, 0));
	}

	while(mapUnvisitedNode.size())
	{

		//find the node in unvisited set with the maximum path capacity.
		int next = -1;
		double max = -1.;
		for(map<int,double>::iterator i = mapUnvisitedNode.begin(); i!=mapUnvisitedNode.end(); i++)
		{
			if(i->second > max)
			{
				max = i->second;
				next = i->first;
			}
		}

		//disconnected from source to dest.
		if(next == -1) 
			return FALSE;
		
		//reach dest, output path and capacity.
		if(next == dest)
		{
			pathcapacity = max;
			
			while(1)
			{
				vecPath.push_back(next);
				int preNode = mapPreNode[next];
				if(preNode == source) 
				{
					vecPath.push_back(preNode);
					break;
				}
				next = mapPreNode[next];
			}
			return TRUE;
		}


		//compute path capacity and update the value of unvisited node set
		for(map<int,Parameters>::iterator i = mapNodes[next].NeighborContactInfo.begin(); 
				i != mapNodes[next].NeighborContactInfo.end(); i++)
		{
			if(mapUnvisitedNode.count(i->first) == 0) continue;

			if(i->second.beta >= max && max != 0.0)
			{
				if(mapUnvisitedNode[i->first] < max)
				{
					mapUnvisitedNode[i->first] = max;
					mapPreNode[i->first] = next;
				}
			}
			else if(i->second.beta > max && max == 0.0)
			{
				if(mapUnvisitedNode[i->first] < i->second.beta)
				{
					mapUnvisitedNode[i->first] = i->second.beta;
					mapPreNode[i->first] = next;
				}
			}
			else if(i->second.beta < max && i->second.beta != 0.0)
			{
				if(mapUnvisitedNode[i->first] < i->second.beta)
				{
					mapUnvisitedNode[i->first] = i->second.beta;
					mapPreNode[i->first] = next;
				}
			}
		}

		//remove visited node
		mapUnvisitedNode.erase(next);
	}
}


BOOL HeuristicAlgorithm::DijkstraContactProbability(map<int, Node> & mapNodes, int source, int dest, double deadline, 
		double & probability, double & pathcapacity, vector<int> & vecPath)
{
	pathcapacity = 0.0;
	probability = 0.0;
	vecPath.clear();

	map<int, double> mapUnvisitedNode;
	map<int, int> mapPreNode;

	for(map<int,Node>::iterator i = mapNodes.begin(); i != mapNodes.end(); i++)
	{
		if(i->first != source)
			mapUnvisitedNode.insert(pair<int,double>(i->first, -1.0));
		else
			mapUnvisitedNode.insert(pair<int,double>(i->first, 0));
		mapPreNode.insert(pair<int,double>(i->first, 0));
	}

	while(mapUnvisitedNode.size())
	{

		//find the node in unvisited set with the maximum contact probability.
		int next = -1;
		double max = -1.;
		for(map<int,double>::iterator i = mapUnvisitedNode.begin(); i!=mapUnvisitedNode.end(); i++)
		{
			if(i->second > max)
			{
				max = i->second;
				next = i->first;
			}
		}

		//disconnected from source to dest.
		if(next == -1) 
			return FALSE;
		
		//reach dest, output path capacity and path, allocate resources.
		if(next == dest)
		{
			probability = max;

			//find the path with maximum probability
			while(1)
			{
				vecPath.push_back(next);
				int preNode = mapPreNode[next];
				if(mapNodes[next].NeighborContactInfo[preNode].beta < pathcapacity || pathcapacity==0.0)
				{
					pathcapacity = mapNodes[next].NeighborContactInfo[preNode].beta;
				}
				if(preNode == source) 
				{
					vecPath.push_back(preNode);
					break;
				}
				next = mapPreNode[next];
			}

			return TRUE;
		}


		//compute probability and update the value of unvisited node set
		for(map<int,Parameters>::iterator i = mapNodes[next].NeighborContactInfo.begin(); 
				i != mapNodes[next].NeighborContactInfo.end(); i++)
		{
			//skip the disconnected node
			if(i->second.removed == true) continue;
			//skip the visited node
			if(mapUnvisitedNode.count(i->first) == 0) continue;
			//retrieve the path from source to current node
			if(i->second.beta <= 0.0) continue;
			vector<int> vecTempPath;
			vecTempPath.push_back(i->first);
			vecTempPath.push_back(next);
			for(int j=mapPreNode[next]; j!=0; j=mapPreNode[j])
			{
				vecTempPath.push_back(j);
			}

			//calculate and update probability 
			double tempProb = CalculateContactProbability(mapNodes,vecTempPath,deadline);

		
			if(mapUnvisitedNode[i->first] < tempProb)
			{
				mapUnvisitedNode[i->first] = tempProb;
				mapPreNode[i->first] = next;
			}
		}

		//remove visited node
		mapUnvisitedNode.erase(next);
	}
}


BOOL HeuristicAlgorithm::DijkstraUploadingProbability(map<int, Node> & mapNodes, int source, int dest, double deadline, 
		double threshold, double & probability, double & pathcapacity, vector<int> & vecPath)
{
	pathcapacity = 0.0;
	probability = 0.0;
	vecPath.clear();

	map<int, double> mapUnvisitedNode;
	map<int, int> mapPreNode;

	for(map<int,Node>::iterator i = mapNodes.begin(); i != mapNodes.end(); i++)
	{
		if(i->first != source)
			mapUnvisitedNode.insert(pair<int,double>(i->first, -1.0));
		else
			mapUnvisitedNode.insert(pair<int,double>(i->first, 0));
		mapPreNode.insert(pair<int,double>(i->first, 0));
	}

	while(mapUnvisitedNode.size())
	{

		//find the node in unvisited set with the maximum contact probability.
		int next = -1;
		double max = -1.;
		for(map<int,double>::iterator i = mapUnvisitedNode.begin(); i!=mapUnvisitedNode.end(); i++)
		{
			if(i->second > max)
			{
				max = i->second;
				next = i->first;
			}
		}

		//disconnected from source to dest.
		if(next == -1) 
			return FALSE;
		
		//reach dest, output path capacity and path, allocate resources.
		if(next == dest)
		{
			probability = max;
			if(probability < threshold) return FALSE;

			//find the path with maximum probability
			while(1)
			{
				vecPath.push_back(next);
				int preNode = mapPreNode[next];
				double beta = mapNodes[next].NeighborContactInfo[preNode].beta;
				double alpha = mapNodes[next].NeighborContactInfo[preNode].alpha;
				if(alpha*beta/(alpha-1) < pathcapacity || pathcapacity==0.0)
				{
					pathcapacity = alpha*beta/(alpha-1);
				}
				if(preNode == source) 
				{
					vecPath.push_back(preNode);
					break;
				}
				next = mapPreNode[next];
			}
			return TRUE;
		}


		//compute probability and update the value of unvisited node set
		for(map<int,Parameters>::iterator i = mapNodes[next].NeighborContactInfo.begin(); 
				i != mapNodes[next].NeighborContactInfo.end(); i++)
		{
			//skip the disconnected node
			if(i->second.removed == true) continue;
			//skip the visited node
			if(mapUnvisitedNode.count(i->first) == 0) continue;
			//retrieve the path from source to current node
			if(i->second.beta <= 0.0) continue;

			vector<int> vecTempPath;
			vector<int> vecContact;
			vecTempPath.push_back(i->first);
			vecTempPath.push_back(next);
			vecContact.push_back(1);
			for(int j=mapPreNode[next]; j!=0; j=mapPreNode[j])
			{
				vecTempPath.push_back(j);
				vecContact.push_back(1);
			}

			//calculate and update probability 
			double dbContactProb = ApproximateContatcProbability(mapNodes,vecTempPath,vecContact,deadline);

			//find the capacity (minimum mean) of the path
			double dbCapacity = -1.0;
			for(int j=0; j<vecTempPath.size()-1; j++)
			{
				double beta = mapNodes[vecTempPath[j]].NeighborContactInfo[vecTempPath[j+1]].beta;
				double alpha = mapNodes[vecTempPath[j]].NeighborContactInfo[vecTempPath[j+1]].alpha;
				double mu = alpha*beta/(alpha-1);
				if(mu < dbCapacity || dbCapacity == -1.0)
					dbCapacity = mu;
			}

		
			if(mapUnvisitedNode[i->first] < dbContactProb * dbCapacity)
			{
				mapUnvisitedNode[i->first] = dbContactProb * dbCapacity;
				mapPreNode[i->first] = next;
			}
		}

		//remove visited node
		mapUnvisitedNode.erase(next);
	}
}



double HeuristicAlgorithm::PathCapacityPreferred(map<int, Node> & mapNodes, vector<vector<int>> & vecPaths, vector<double> & vecCapacities, 
		vector<double> & vecProbs, double deadline, int source, int dest)
{
	vecProbs.clear();
 	
	for(vector<vector<int>>::iterator i=vecPaths.begin(); i!=vecPaths.end(); i++)
	{
		vector<int> vecContacts;
		for(int j=0; j<i->size()-1; j++)
		{
			vecContacts.push_back(1);
		}

		double prob = CalculateContactProbability(mapNodes, *i, deadline);
		double appprob = ApproximateContatcProbability(mapNodes, *i, vecContacts, deadline);
		vecProbs.push_back(prob);
	}

	while(1)
	{
		//find the path with minimum contact probability
		double prob=1.0;
		int index;
		for(int i=0; i<vecProbs.size(); i++)
		{
			if(vecProbs[i]<=prob)
			{
				prob = vecProbs[i];
				index = i;
			}
		}

		double pathcapacity = vecCapacities[index];
		vector<int> vecPath = vecPaths[index];

		//Release the resource of found path
		ReleaseResourceInclusive(mapNodes, vecPath, pathcapacity);


		//Dijstra with maximizing contact probability iteratively
		double newprob, newtotalprob=1.0, newcapacity, newtotalcapacity=0.0;
		vector<int> vecNewpath;
		vector<vector<int>> vecNewPaths;
	    vector<double> vecNewCapacities;
		vector<double> vecNewProbabilities;

		BOOL bReplace = FALSE;

		while(DijkstraContactProbability(mapNodes, source, dest, deadline, newprob, newcapacity, vecNewpath))
		{
			//exclude the path capacity allocated by the found path
			AllocateResourceInclusive(mapNodes, vecNewpath, newcapacity);

			newtotalcapacity += newcapacity;
			newtotalprob *= newprob;
			vecNewCapacities.push_back(newcapacity);
			vecNewPaths.push_back(vecNewpath);
			vecNewProbabilities.push_back(newprob);

			if(newtotalcapacity >= pathcapacity && newtotalprob >= prob)
			{
				bReplace = TRUE;
				break;
			}
			else if(newtotalcapacity < pathcapacity && newtotalprob >= prob)
			{
				continue;
			}
			else if(newtotalprob < prob)
			{
				bReplace = FALSE;
				break;
			}
		}

		if(bReplace)
		{
			//delete the previous path
			vecCapacities.erase(vecCapacities.begin()+index);
			vecPaths.erase(vecPaths.begin()+index);
			vecProbs.erase(vecProbs.begin()+index);

			//include the found paths
			for(size_t i=0; i<vecNewPaths.size(); i++)
			{
				vecPaths.push_back(vecNewPaths[i]);
				vecProbs.push_back(vecNewProbabilities[i]);
				vecCapacities.push_back(vecNewCapacities[i]);
			}
		}
		else
		{
			for(int i=0; i<vecNewPaths.size(); i++)
			{
				ReleaseResourceInclusive(mapNodes, vecNewPaths[i], vecNewCapacities[i]);
			}
			AllocateResourceInclusive(mapNodes, vecPaths[index], vecCapacities[index]); 
			break;
		}
	}

	double totalprob = 1.0;
	for(int i=0; i<vecProbs.size(); i++)
		totalprob *= vecProbs[i];
	return totalprob;
}


double HeuristicAlgorithm::ContactProbabilityPreferred(map<int, Node> & mapNodes, vector<vector<int>> & vecPaths, 
		vector<double> & vecCapacities, vector<double> & vecProbs, double datasize, double deadline, int source, int dest)
{
	vecPaths.clear();
	vecCapacities.clear();
	vecProbs.clear();

	vector<int> vecPath;
	double pathcapacity;
	double probability;
	double totalcapacity=0.0;
	double totalprobability = 1.0;

	while(DijkstraContactProbability(mapNodes, source, dest, deadline, probability, pathcapacity, vecPath)
		&& totalcapacity < datasize)
	{
		vecPaths.push_back(vecPath);
		vecCapacities.push_back(pathcapacity);
		vecProbs.push_back(probability);
		totalcapacity += pathcapacity;
		totalprobability *= probability;
		AllocateResourceExclusive(mapNodes, vecPath);
	}

	if(totalcapacity >= datasize) return totalprobability;
	vector<double> vecContactProbs = vecProbs;

	AssignData2Paths(mapNodes, datasize-totalcapacity, vecPaths,vecContactProbs,vecCapacities,vecProbs);

	/*relocate the data assigned at the path with minimum uploading probability to 
	  other paths to achieve better performance*/
	while(1)
	{
		//find the path with minimum uploading probability
		double minprob=1.0;
		int index=-1;
		double totalprob=1.0;
		for(int i=0; i<vecProbs.size(); i++)
		{
			totalprob *= vecProbs[i];
			if(vecProbs[i] <= minprob)
			{
				minprob = vecProbs[i];
				index = i;
			}
		}
		//release the resource taken by the path
		ReleaseResourceExclusive(mapNodes, vecPaths[index]);

		double relocatedcapacity = vecCapacities[index];
		vector<double> vecNewCapacities = vecCapacities;
		vector<double> vecNewProbs = vecProbs;
		vector<vector<int>> vecNewPaths = vecPaths;
		vector<double> vecNewContactProbs = vecContactProbs;
		vecNewCapacities.erase(vecNewCapacities.begin()+index);
		vecNewProbs.erase(vecNewProbs.begin()+index);
		vecNewPaths.erase(vecNewPaths.begin()+index);
		vecNewContactProbs.erase(vecNewContactProbs.begin()+index);

		double newtotalprob = AssignData2Paths(mapNodes,relocatedcapacity,vecNewPaths,vecNewContactProbs,vecNewCapacities,vecNewProbs);

		if(newtotalprob > totalprob)
		{
			vecPaths = vecNewPaths;
			vecContactProbs = vecNewContactProbs;
			vecCapacities = vecNewCapacities;
			vecProbs = vecNewProbs;
		}
		else
			return totalprob;
	}
}


double HeuristicAlgorithm::JointConsiderationPreferred(map<int, Node> & mapNodes, vector<vector<int>> & vecPaths, vector<double> & vecDataAssignment, 
			vector<double> & vecProbs, double datasize, double deadline, int source, int dest, double thresholdprob, bool bProportional /*=false*/)
{
	vecPaths.clear();
	vecDataAssignment.clear();
	vecProbs.clear();
	
	vector<double> vecCapacities;
	vector<double> vecPathMean;
	vector<int> vecPath;
	vector<double> vecMetrics;
	double pathmean;
	double probability;
	double pathcapacity=0.0;
	double totalcapacity=0.0;
	
	//find the paths that have higher uploading probability than source
	while(DijkstraUploadingProbability(mapNodes, source, dest, deadline, thresholdprob, probability, pathmean, vecPath))
	{	
		vecPaths.push_back(vecPath);
		vecPathMean.push_back(pathmean);
		vecProbs.push_back(probability/pathmean);
		vecMetrics.push_back(probability);
		AllocateResourceExclusive(mapNodes, vecPath);

		for(int i=0; i<vecPath.size()-1; i++)
		{
			double beta = mapNodes[vecPath[i]].NeighborContactInfo[vecPath[i+1]].beta;
			if(beta < pathcapacity || pathcapacity == 0.0) pathcapacity = beta;
		}
		if(totalcapacity + pathcapacity > datasize)
		{
			pathcapacity = datasize - totalcapacity;
			vecCapacities.push_back(pathcapacity);
			totalcapacity = datasize;
			break;
		}
		else
		{
			vecCapacities.push_back(pathcapacity);
			totalcapacity += pathcapacity;
		}
	}

	if(totalcapacity < datasize && bProportional == false) 
		ProcessiveAssignData2Paths(mapNodes, deadline, datasize-totalcapacity, vecPaths, vecCapacities, vecProbs);
	if(totalcapacity < datasize && bProportional == true)
		ProportionalAssignData2Paths(mapNodes, deadline, datasize, vecPaths, vecMetrics, vecCapacities, vecProbs);

	/*relocate the data assigned at the path with minimum uploading probability to 
	  other paths to achieve better performance*/
	while(1)
	{


		//find the path with minimum uploading probability
		double minprob=1.0;
		int index=-1;
		double totalprob=1.0;
		for(int i=0; i<vecProbs.size(); i++)
		{
			totalprob *= vecProbs[i];
			if(vecProbs[i] <= minprob)
			{
				minprob = vecProbs[i];
				index = i;
			}
		}

		if(vecPaths.size() == 1) return totalprob;

		//release the resource taken by the path
		ReleaseResourceExclusive(mapNodes, vecPaths[index]);

		double relocatedcapacity = vecCapacities[index];
		vector<double> vecNewCapacities = vecCapacities;
		vector<double> vecNewProbs = vecProbs;
		vector<vector<int>> vecNewPaths = vecPaths;
	
		vecNewCapacities.erase(vecNewCapacities.begin()+index);
		vecNewProbs.erase(vecNewProbs.begin()+index);
		vecNewPaths.erase(vecNewPaths.begin()+index);

		double newtotalprob = ProcessiveAssignData2Paths(mapNodes,deadline, relocatedcapacity,vecNewPaths, vecNewCapacities,vecNewProbs);

		if(newtotalprob > totalprob)
		{
			vecPaths = vecNewPaths;
			vecCapacities = vecNewCapacities;
			vecProbs = vecNewProbs;
		}
		else
			return totalprob;
	}
}


double HeuristicAlgorithm::CalculateContactProbability(map<int, Node> & mapNodes, vector<int> vecPath, double deadline /*Hours*/)
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


double HeuristicAlgorithm::ApproximateContatcProbability(map<int, Node> & mapNodes, vector<int> vecPath, vector<int> vecContacts, double deadline, double dataassigment /*=0*/)
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


double HeuristicAlgorithm::AssignData2Paths(map<int, Node> & mapNodes, double dataallocated, vector<vector<int>> & vecPaths, 
		vector<double> &vecContactProbs, vector<double> & vecCapacities, vector<double> &vecProbs)
{
	//assign the remaining data
	while(1)
	{
		//find the path with maximum contact probability
		double maxprob=0.0;
		int index=-1;
		for(int i=0; i<vecProbs.size(); i++)
		{
			if(vecProbs[i] > maxprob)
			{
				maxprob = vecProbs[i];
				index = i;
			}
		}

		//calculate the amount of data assignment
		double oldcapacity = vecCapacities[index];
		double newcapacity = oldcapacity;
		set<double> setRank;
		for(int i=0; i<vecPaths[index].size()-1; i++)
		{
			setRank.insert(mapNodes[vecPaths[index][i]].NeighborContactInfo[vecPaths[index][i+1]].beta);
		}

		for(set<double>::iterator i=setRank.begin(); i!=setRank.end(); i++)
		{
			if(*i > oldcapacity)
			{
				if(*i - oldcapacity >= dataallocated)
				{
					newcapacity += dataallocated;
					dataallocated = 0.0;
				}
				else
				{
					newcapacity = *i;
					dataallocated -= *i - oldcapacity;
				}
				break;
			}
		}
		if(oldcapacity == newcapacity)
		{
			if(dataallocated <= *setRank.begin())
			{
				newcapacity += dataallocated;
				dataallocated = 0.0;
			}
			else
			{
				newcapacity += *setRank.begin();
				dataallocated -= *setRank.begin();
			}
		}

		double contactprob = vecContactProbs[index];
		for(int i=0; i<vecPaths[index].size()-1; i++)
		{
			double alpha = mapNodes[vecPaths[index][i]].NeighborContactInfo[vecPaths[index][i+1]].alpha;
			double beta = mapNodes[vecPaths[index][i]].NeighborContactInfo[vecPaths[index][i+1]].beta;
			if(beta > newcapacity)
				contactprob *= 1.0;
			else
				contactprob *= pow(beta/newcapacity,alpha);
		}

		vecProbs[index]=contactprob;
		vecCapacities[index]=newcapacity;

		if(dataallocated <= 0.0) 
			break;
	}

	double totalprob = 1.0;
	for(int i=0; i<vecProbs.size(); i++)
	{
		totalprob *= vecProbs[i];
	}
	return totalprob;
}


double HeuristicAlgorithm::ProcessiveAssignData2Paths(map<int, Node> & mapNodes, double deadline, 
		double dataallocated, vector<vector<int>> & vecPaths, vector<double> & vecCapacities, vector<double> &vecProbs)
{
	//assign the remaining data
	while(1)
	{
		//find the path with maximum contact probability
		double maxprob=0.0;
		int index=-1;
		for(int i=0; i<vecProbs.size(); i++)
		{
			if(vecProbs[i] > maxprob)
			{
				maxprob = vecProbs[i];
				index = i;
			}
		}

		//calculate the amount of data assignment
		double oldcapacity = vecCapacities[index];
		double newcapacity = oldcapacity;
		set<double> setRank;
		for(int i=0; i<vecPaths[index].size()-1; i++)
		{
			setRank.insert(mapNodes[vecPaths[index][i]].NeighborContactInfo[vecPaths[index][i+1]].beta);
		}

		for(set<double>::iterator i=setRank.begin(); i!=setRank.end(); i++)
		{
			if(*i > oldcapacity)
			{
				if(*i - oldcapacity >= dataallocated)
				{
					newcapacity += dataallocated;
					dataallocated = 0.0;
				}
				else
				{
					newcapacity = *i;
					dataallocated -= *i - oldcapacity;
				}
				break;
			}
		}
		if(oldcapacity == newcapacity)
		{
			if(dataallocated <= *setRank.begin())
			{
				newcapacity += dataallocated;
				dataallocated = 0.0;
			}
			else
			{
				newcapacity += *setRank.begin();
				dataallocated -= *setRank.begin();
			}
		}

		vecCapacities[index]=newcapacity;

		double contactprob=0.0; 
		vector<int> vecContacts, vecTemp;
		vector<vector<int>> vecCombinations;

		for(int i=0; i<vecPaths[index].size()-1; i++)
		{
			double beta = mapNodes[vecPaths[index][i]].NeighborContactInfo[vecPaths[index][i+1]].beta;
			vecContacts.push_back((int)ceil(newcapacity/beta));
		}

		RecursiveCombination(vecContacts, vecTemp, vecCombinations);

		for(int i=0; i<vecCombinations.size(); i++)
		{
			double dbprob = ApproximateContatcProbability(mapNodes, vecPaths[index],vecCombinations[i], deadline, newcapacity);
			if(dbprob > contactprob) contactprob = dbprob;
		}

		vecProbs[index]=contactprob;
		
		if(dataallocated <= 0.0) 
			break;
	}

	double totalprob = 1.0;
	for(int i=0; i<vecProbs.size(); i++)
	{
		totalprob *= vecProbs[i];
	}
	return totalprob;
}


double HeuristicAlgorithm::ProportionalAssignData2Paths(map<int, Node> & mapNodes, double deadline, 
		double dataallocated, vector<vector<int>> & vecPaths, vector<double> &vecMetrics, vector<double> & vecCapacities, vector<double> &vecProbs)
{
	//assign the data
	double totalmetrics = 0.0;

	for(int i=0; i<vecPaths.size(); i++)
	{
		totalmetrics += vecMetrics[i];
	}

	for(int i=0; i<vecCapacities.size(); i++)
	{
		vecCapacities[i] = dataallocated*vecMetrics[i]/totalmetrics;

		double contactprob=0.0; 
		vector<int> vecContacts, vecTemp;
		vector<vector<int>> vecCombinations;

		for(int j=0; j<vecPaths[i].size()-1; j++)
		{
			double beta = mapNodes[vecPaths[i][j]].NeighborContactInfo[vecPaths[i][j+1]].beta;
			vecContacts.push_back((int)ceil(vecCapacities[i]/beta));
		}

		RecursiveCombination(vecContacts, vecTemp, vecCombinations);

		for(int j=0; j<vecCombinations.size(); j++)
		{
			double dbprob = ApproximateContatcProbability(mapNodes, vecPaths[i], vecCombinations[j], deadline, vecCapacities[i]);
			if(dbprob > contactprob) contactprob = dbprob;
		}

		vecProbs[i]=contactprob;
	}


	double totalprob = 1.0;
	for(int i=0; i<vecProbs.size(); i++)
	{
		totalprob *= vecProbs[i];
	}
	return totalprob;
}


void HeuristicAlgorithm::ReleaseResourceInclusive(map<int,Node> & mapNodes, vector<int> vecPath, double pathcapacity)
{
	for(int i=0,j=1; i<vecPath.size()-1; i++,j++)
	{
		mapNodes[vecPath[i]].NeighborContactInfo[vecPath[j]].beta += pathcapacity;
		mapNodes[vecPath[j]].NeighborContactInfo[vecPath[i]].beta += pathcapacity;
	}
}


void HeuristicAlgorithm::ReleaseResourceExclusive(map<int,Node> & mapNodes, vector<int> vecPath)
{
	for(int i=0,j=1; i<vecPath.size()-1; i++,j++)
	{
		mapNodes[vecPath[i]].NeighborContactInfo[vecPath[j]].removed = false;
		mapNodes[vecPath[j]].NeighborContactInfo[vecPath[i]].removed = false;
	}
}


void HeuristicAlgorithm::AllocateResourceInclusive(map<int,Node> & mapNodes, vector<int> vecPath, double pathcapacity)
{
	for(int i=0,j=1; i<vecPath.size()-1; i++,j++)
	{
		mapNodes[vecPath[i]].NeighborContactInfo[vecPath[j]].beta -= pathcapacity;
		mapNodes[vecPath[j]].NeighborContactInfo[vecPath[i]].beta -= pathcapacity;
	}
}


void HeuristicAlgorithm::AllocateResourceExclusive(map<int,Node> & mapNodes, vector<int> vecPath)
{
	for(int i=0, j=1; i<vecPath.size()-1; i++,j++)
	{
//		mapNodes[vecPath[i]].NeighborContactInfo.erase(vecPath[j]);
//		mapNodes[vecPath[j]].NeighborContactInfo.erase(vecPath[i]);
		mapNodes[vecPath[i]].NeighborContactInfo[vecPath[j]].removed = true;
		mapNodes[vecPath[j]].NeighborContactInfo[vecPath[i]].removed = true;

	}
}

