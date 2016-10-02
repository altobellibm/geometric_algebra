// AlgebraGemetrica.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Multivetor.h"

int main(int argc, char* argv[])
{
	Multivetor<double> A = 1.0*e(1) + 1.0*e(2);
	Multivetor<double> B = 1.0*e(3) + 1.0*e(2);
	Multivetor<double> C = A^B;

	return 0;
	
}

