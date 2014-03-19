#include "stdafx.h"
#include "mathfunctions.h"


int factorial(int i)
{
	if(i==0) return 1;

	int value = 1;
	for(int j=1; j<=i; j++)
		value*=j;

	return value;
}