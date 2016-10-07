// AlgebraGemetrica.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Multivector.h"

int main(int argc, char* argv[])
{
	Multivector<double> A = 1.0*e(1) + 1.0*e(2);
	Multivector<double> B = 1.0*e(3) + 1.0*e(2);
	Multivector<double> C = A^B;

	

	return 0;
	
}

