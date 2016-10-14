// AlgebraGemetrica.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include "Multivector.h"

int main(int , char** )
{
	Orthonormal ortonormal(4);

    // Exercício 1
        // Item a)
        std::cout << ((e(1)+e(2))^(e(3)+e(2))) << std::endl;

        // Item b)
        std::cout << ((e(2) - e(1)) ^ (e(1) - 2.0*e(3)))  << std::endl;

        // Item c)
        std::cout << ((4.0*e(1) + e(2) + e(3)) ^ (3.0*e(1))) << std::endl;

        // Item d)
        std::cout << ((e(2) + e(3)) ^ ((1.0/2.0)*e(1) + e(2) + (3.0/2.0)*e(3))) << std::endl;

        // Item e)
        std::cout << ((e(1) + e(2)) ^ ((e(2)^e(1)) + (e(3)^e(2)))) << std::endl;



	return 0;
	
}

