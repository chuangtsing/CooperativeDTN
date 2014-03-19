#include "stdafx.h"
#include "ProcessMIT.h"


ProcessMIT::ProcessMIT(void)
{
	vecContacts.clear();
	mapEdges.clear();
	FILE * fr = fopen("MITContact.csv", "r");
	FILE * fw = fopen("MITTrace.dat", "w");
	if(fr == NULL) return;
	if(fw == NULL) return;
	CStdioFile file(fr);
	CString strLine = _T("");
	file.ReadString(strLine);
	CTime tTraceStart (2004, 7, 10, 15, 57, 0);
	CTime tTraceEnd(2005, 5, 5, 3, 23, 0);
	//columns explanation 
	//oid, endtime, starttime, person_oid, device_oid 
	while(!feof(fr))
	{
		file.ReadString(strLine);
		strLine.Trim();
		if(strLine == _T("")) continue;
		int oid, EndMonth,EndDate, EndYear, EndHour, EndMinite, StartMonth, StartDate, StartYear, StartHour, StartMinite, Person_oid, Device_oid;
		int result = swscanf_s(strLine.GetBuffer(strLine.GetLength()), L"%d,%d/%d/%d %d:%d,%d/%d/%d %d:%d,%d,%d", 
			&oid, &EndMonth,&EndDate, &EndYear, &EndHour, &EndMinite, &StartMonth, &StartDate, &StartYear, &StartHour, &StartMinite, &Person_oid, &Device_oid);
		if(result != 13 || Device_oid > 97) continue;
		if(Device_oid == Person_oid) continue;
		CTime EndTime(EndYear, EndMonth, EndDate, EndHour, EndMinite, 0);
		CTime StarTime(StartYear, StartMonth, StartDate, StartHour, StartMinite, 0);
		Contact tempContact;
		tempContact.source = Person_oid;
		tempContact.dest = Device_oid;
		tempContact.starttime = StarTime.GetTime();
		tempContact.endtime = EndTime.GetTime();
		tempContact.duration = EndTime.GetTime()-StarTime.GetTime();
		if(tempContact.duration <= 0) continue;
		if(StarTime.GetTime() > (tTraceEnd.GetTime()-tTraceStart.GetTime())/2+tTraceStart.GetTime())
		{
			vecContacts.push_back(tempContact);
			fprintf(fw, "%-6d %-6d %-6I64d %-6d\n", tempContact.source, tempContact.dest, tempContact.starttime, tempContact.duration);
			continue;
		}

		BOOL bFound = FALSE;
		
		if(mapEdges.count(Person_oid+Device_oid) > 0)
		{
			multimap<int, Edge>::iterator it;
			pair<multimap<int, Edge>::iterator, multimap<int, Edge>::iterator> ret;
			ret = mapEdges.equal_range(Person_oid+Device_oid);
			for(it=ret.first; it!=ret.second; it++)
			{
				if((it->second.source == Person_oid && it->second.dest == Device_oid) 
					|| (it->second.source == Device_oid && it->second.dest == Person_oid))
				{
					it->second.counter++;
					it->second.sumofduration+=tempContact.duration;
					it->second.sumofinterarrival += tempContact.starttime-it->second.lasttimemeet;
					it->second.lasttimemeet=tempContact.starttime;
					if(it->second.beta > tempContact.duration)
						it->second.beta = tempContact.duration;
					bFound = TRUE;
					break;
				}
			}
		}
		if(!bFound)
		{
			Edge tempEdge;
			tempEdge.source = Person_oid;
			tempEdge.dest = Device_oid;
			tempEdge.counter = 1;
			tempEdge.sumofduration = tempContact.duration;
			tempEdge.sumofinterarrival = 0;
			tempEdge.lasttimemeet = tempContact.starttime;
			tempEdge.beta = tempContact.duration;
			mapEdges.insert(pair<int, Edge>(Person_oid+Device_oid, tempEdge));
		}
		
	}
	fclose(fr);
	fclose(fw);

	multimap<int,Edge>::iterator it;
	for(it=mapEdges.begin(); it!=mapEdges.end(); )
	{
		if(it->second.counter < 2 || it->second.sumofinterarrival <= 300) //training the trace
		{
			mapEdges.erase(it++);
		}
		else
		{
			it->second.lambda = 1.0/(it->second.sumofinterarrival/(3600.0*(it->second.counter-1))); //interarrival in Hours
			it->second.beta /= 60.0;
			double mu = it->second.sumofduration/60.0/it->second.counter;
			it->second.alpha = mu/(mu-it->second.beta);  //contact duration in Minutes
			++it;
		}
	}

}


ProcessMIT::~ProcessMIT(void)
{
}


void ProcessMIT::ExportNetwork()
{
	FILE * file = fopen("MITNetwork.dat", "w");
	if(file == NULL) return;
	multimap<int, Edge>::iterator it;
	for (it=mapEdges.begin(); it!=mapEdges.end(); it++)
	{
		fprintf(file, "%-6d %-6d %-6f %-6f %-6f\n", (*it).second.source, (*it).second.dest, (*it).second.alpha, (*it).second.beta, (*it).second.lambda);
	}
	fclose(file);
}


void ProcessMIT::StructNodesFromEdges()
{
	mapNodes.clear();

	for(multimap<int,Edge>::iterator i=mapEdges.begin(); i!=mapEdges.end(); i++)
	{
		Parameters para;
		para.alpha = i->second.alpha;
		para.beta = i->second.beta;
		para.lambda = i->second.lambda;
		para.removed = 0;
		
		if(mapNodes.count(i->second.source)==0)
		{
			Node node;
			node.NeighborContactInfo.insert(pair<int,Parameters>(i->second.dest,para));
			mapNodes.insert(pair<int, Node>(i->second.source,node));
		}
		else
		{
			mapNodes.find(i->second.source)->second.NeighborContactInfo.insert(pair<int,Parameters>(i->second.dest,para));
		}

		if(mapNodes.count(i->second.dest)==0)
		{
			Node node;
			node.NeighborContactInfo.insert(pair<int,Parameters>(i->second.source,para));
			mapNodes.insert(pair<int, Node>(i->second.dest,node));
		}
		else
		{
			mapNodes.find(i->second.dest)->second.NeighborContactInfo.insert(pair<int,Parameters>(i->second.source,para));
		}
	}
}


BOOL ProcessMIT::StructNodesFromFile()
{
	FILE * fp = fopen("MITNetwork.dat", "r");
	if(fp == NULL) return FALSE;
	CStdioFile file(fp);
	CString strLine = _T("");

	mapNodes.clear();

	while(!feof(fp))
	{
		file.ReadString(strLine);
		strLine.Trim();
		if (strLine == "") continue;
		int source, dest;
		double alpha, beta, lambda;
		int result = swscanf_s(strLine.GetBuffer(strLine.GetLength()), L"%d %d %lf %lf %lf", &source, &dest, &alpha, &beta, &lambda);
		if(result != 5) continue;

		Parameters para;
		para.alpha = alpha;
		para.beta = beta;
		para.lambda = lambda;
		para.removed = false;

		if(mapNodes.count(source)==0)
		{
			Node node;
			node.NeighborContactInfo.insert(pair<int,Parameters>(dest,para));
			mapNodes.insert(pair<int, Node>(source,node));
		}
		else
		{
			mapNodes.find(source)->second.NeighborContactInfo.insert(pair<int,Parameters>(dest,para));
		}

		if(mapNodes.count(dest)==0)
		{
			Node node;
			node.NeighborContactInfo.insert(pair<int,Parameters>(source,para));
			mapNodes.insert(pair<int, Node>(dest,node));
		}
		else
		{
			mapNodes.find(dest)->second.NeighborContactInfo.insert(pair<int,Parameters>(source,para));
		}
	}
	fclose(fp);
	return TRUE;
}