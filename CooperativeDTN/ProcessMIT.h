#pragma once

class ProcessMIT
{
public:
	ProcessMIT(void);
	~ProcessMIT(void);
	
	vector<Contact> vecContacts;
	multimap<int, Edge> mapEdges;
	map<int, Node> mapNodes;

	void ExportNetwork();
	void StructNodesFromEdges();
	BOOL StructNodesFromFile();
};

