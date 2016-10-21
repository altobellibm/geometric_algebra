// AlgebraGemetrica.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <tchar.h>
#include "Multivector.h"

//_tsetlocale(LC_ALL, _T("portuguese"));

int main(int , char** )
{
	    // Exercício 1

		std::cout << "Questão 1." << std::endl;

        std::cout << "Item a) " << std::endl;
		std::cout << ((e(1)+e(2))^(e(3)+e(2))) << std::endl;

		std::cout << "Item b) " << std::endl;
        std::cout << ((e(2) - e(1)) ^ (e(1) - 2.0*e(3)))  << std::endl;

		std::cout << "Item c) " << std::endl;
        std::cout << ((4.0*e(1) + e(2) + e(3)) ^ (3.0*e(1))) << std::endl;
		
		std::cout << "Item d) " << std::endl;
        std::cout << ((e(2) + e(3)) ^ ((1.0/2.0)*e(1) + e(2) + (3.0/2.0)*e(3))) << std::endl;

		std::cout << "Item e) " << std::endl;
        std::cout << ((e(1) + e(2)) ^ ((e(2)^e(1)) + (e(3)^e(2)))) << std::endl;

	// Exercício 2 

		std::cout << std::endl;
		std::cout << "Questão 2." << std::endl;

		std::cout << "Item a) " << std::endl;
		std::cout << (e(1) ^ (e(1) ^ (e(2) + 2 * e(3)) ^ e(4))) << std::endl;

		std::cout << "Item b) " << std::endl;
		std::cout << ((e(1) - 3 * e(4)) ^ (e(1) ^ (e(2) + 2 * e(3)) ^ e(4))) << std::endl;

		std::cout << "Item c) " << std::endl;
		std::cout << ((e(2) + e(3)) ^ (e(1) ^ (e(2) + 2 * e(3)) ^ e(4))) << std::endl;


	// Exercício 3 

		std::cout << std::endl;
		std::cout << "Questão 3." << std::endl;

		std::cout << ( (2*e(2) + e(3)) ^ (e(2) - 1*e(3)) ) << std::endl;


	// Exercício 4 

		Orthonormal ortonormal(3);

		std::cout << std::endl;
		std::cout << "Questão 4." << std::endl;

		std::cout << "Item a) " << std::endl;
		std::cout << SCP((e(1) + e(3)), (e(1) + e(2)), ortonormal) << std::endl;

	return 0;
	
}

