// CooperativeDTN.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "ProcessMIT.h"
#include "HeuristicAlgorithm.h"

int _tmain(int argc, _TCHAR* argv[])
{

	ProcessMIT MIT;
	MIT.ExportNetwork();
	MIT.StructNodesFromEdges();

	HeuristicAlgorithm heuristic(MIT.mapNodes, 50, 24, 6, 29);

	return 0;
}

