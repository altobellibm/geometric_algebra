// AlgebraGemetrica.cpp : Defines the entry point for the console application.
//

#include "Multivector.h"

int main(int argc, char* argv[])
{


	Orthonormal ortonormal(4);

	Multivector<double> A = (1.0*e(1) ^ 1.0*e(2)) + (1.0*e(1) ^ 1.0*e(3)) + (1.0*e(2) ^ 1.0*e(3));
	Multivector<double> B = (1.0*e(1) ^ 1.0*e(2)) + (1.0*e(1) ^ 1.0*e(4)) + (1.0*e(2) ^ 1.0*e(4));
	Multivector<double> C = GP(A,B,ortonormal);

	//Multivector<double> C = RP(A, B, 4);

	return 0;
	
}

