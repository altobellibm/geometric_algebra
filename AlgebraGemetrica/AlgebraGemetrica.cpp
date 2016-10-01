// AlgebraGemetrica.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Multivetor.h"

static int wordbits[65536] = { /* bitcounts of integers 0 through 65535, inclusive */ };
static int popcount(int i)
{
	return (wordbits[i & 0xFFFF] + wordbits[i >> 16]);
}

int NumberOfSetBits(unsigned int i)
{
	i = i - ((i >> 1) & 0x55555555);
	i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

int main(int argc, char* argv[])
{
	Multivetor<double> A = 1.0*e(1) + 1.0*e(2);
	Multivetor<double> B = 1.0*e(1) + 1.0*e(2);
	Multivetor<double> C = A^B;

	return 0;

	/*int qtd =  popcount(14);
	qtd = NumberOfSetBits(14);

	int a = 14;
	int b = 3;

	int c = a & b;

	unsigned int t = 1;

	t = t << 2;
	unsigned int f = t << 2;


	int x = 1; // 0000 0001

	//int x0 = (x << 3) | (); // 0000 010 Não deslocado
	*/

	
}

