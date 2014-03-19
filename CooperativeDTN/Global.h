#pragma once

#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <set>

using namespace std;

typedef struct{
	int source;
	int dest;
	time_t starttime;
	time_t endtime;
	time_t duration;
}Contact;

typedef struct{
	int source;
	int dest;
	int counter;
	time_t sumofduration;
	time_t sumofinterarrival;
	time_t lasttimemeet;
	double lambda;
	double alpha;
	double beta;
}Edge;

typedef struct{
	double lambda;
	double alpha;
	double beta;
	bool   removed;
}Parameters;

typedef struct Node{
	map<int, Parameters> NeighborContactInfo;
}Node;

typedef struct DistNode{
	int nID;
	int nSource;
	int nDest;
	double dbDataSize;
	double dbDataRemain;
	double dbTimeLeft;
	double dbTimeExpired;
	set<int> setPreNodes;
	vector<vector<int>> vecPaths;
	vector<double> vecDataAssginments;
	vector<double> vecDataUploadProbs;
}DistNode;

