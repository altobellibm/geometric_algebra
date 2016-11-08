// AlgebraGemetrica.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include "Multivector.h"

typedef std::vector<Multivector<int>> BASE_BLADE;
typedef std::vector<std::vector<Multivector<int>>> CLIFFORD_TABLE;

enum class OPERATION { GP, SCP, LCONST, RCONST, OUTP };

void createCliffordTable(const BASE_BLADE& _bladesBase, OPERATION _type, CLIFFORD_TABLE& _table, const Orthonormal* ortonormal = NULL ){

	if ((OPERATION::GP == _type || OPERATION::SCP == _type || OPERATION::LCONST == _type || OPERATION::RCONST == _type) && (ortonormal == NULL)) {
		std::cout << "Must pass metric matrix!" << std::endl;
		return;
	}

	_table.resize(_bladesBase.size(), std::vector<Multivector<int>>(_bladesBase.size()));
 
   for(BASE_BLADE::size_type i = 0; i < _bladesBase.size(); i++)
      for(BASE_BLADE::size_type j = 0; j < _bladesBase.size(); j++){
          switch (_type) {
          case OPERATION::GP: _table[i][j] = GP(_bladesBase[i], _bladesBase[j],*ortonormal); break;
          case OPERATION::SCP: _table[i][j] = SCP(_bladesBase[i], _bladesBase[j],*ortonormal); break;
          case OPERATION::LCONST: _table[i][j] = LConst(_bladesBase[i], _bladesBase[j],*ortonormal); break;
          case OPERATION::RCONST: _table[i][j] = RConst(_bladesBase[i], _bladesBase[j],*ortonormal); break;
		  case OPERATION::OUTP:	_table[i][j] = _bladesBase[i] ^ _bladesBase[j]; break;
          default:
              break;
          }
      }
}

void printTableBladeBase(const std::string& _title, const BASE_BLADE& _bladesBase, const CLIFFORD_TABLE& _table) {

	if (_table.size() != _bladesBase.size())
		return;

	std::cout << _title << std::endl << std::endl;
	for (BASE_BLADE::size_type i = 0; i < _bladesBase.size(); i++){
		std::cout << "Blade reference = " << _bladesBase[i] << std::endl << std::endl;
		for (BASE_BLADE::size_type j = 0; j < _bladesBase.size(); j++)
			std::cout << _bladesBase[i]  << " " << _bladesBase[j] << " = " << _table[i][j] << std::endl;		
		std::cout << std::endl;
	}
}

int main(int , char** )
{
	    // Exercício 1

		std::cout << "Question 1." << std::endl;

        std::cout << "Item a) " << std::endl;
		std::cout << ((e(1)+e(2))^(e(3)+e(2))) << std::endl << std::endl;

		std::cout << "Item b) " << std::endl;
        std::cout << ((e(2) - e(1)) ^ (e(1) - 2.0*e(3)))  << std::endl << std::endl;

		std::cout << "Item c) " << std::endl;
        std::cout << ((4.0*e(1) + e(2) + e(3)) ^ (3.0*e(1))) << std::endl << std::endl;
		
		std::cout << "Item d) " << std::endl;
        std::cout << ((e(2) + e(3)) ^ ((1.0/2.0)*e(1) + e(2) + (3.0/2.0)*e(3))) << std::endl << std::endl;

		std::cout << "Item e) " << std::endl;
        std::cout << ((e(1) + e(2)) ^ ((e(2)^e(1)) + (e(3)^e(2)))) << std::endl << std::endl;

	// Exercício 2 

		std::cout << std::endl;
		std::cout << "Question 2." << std::endl;

		std::cout << "Item a) " << std::endl;
		std::cout << (e(1) ^ (e(1) ^ (e(2) + 2 * e(3)) ^ e(4))) << std::endl << std::endl;

		std::cout << "Item b) " << std::endl;
		std::cout << ((e(1) - 3 * e(4)) ^ (e(1) ^ (e(2) + 2 * e(3)) ^ e(4))) << std::endl << std::endl;

		std::cout << "Item c) " << std::endl;
		std::cout << ((e(2) + e(3)) ^ (e(1) ^ (e(2) + 2 * e(3)) ^ e(4))) << std::endl << std::endl;


	// Exercício 3 

		std::cout << std::endl;
		std::cout << "Question 3." << std::endl;

		std::cout << ( (2*e(2) + e(3)) ^ (e(2) - 1*e(3)) ) << std::endl;


	// Exercício 4 

		Orthonormal ortonormal(3);

        Multivector<int> a = (e(1) + e(3));
        Multivector<int> b = (e(1) + e(2));

		std::cout << std::endl;
		std::cout << "Question 4." << std::endl;

		std::cout << "Item a) " << std::endl;
        std::cout << SCP(a, b, ortonormal) << std::endl << std::endl;

        std::cout << "Item b) " << std::endl;
        std::cout << LConst((e(3)),b, ortonormal) << std::endl << std::endl;

        std::cout << "Item c) " << std::endl;
        std::cout << LConst((e(3)),(a^b),ortonormal) << std::endl << std::endl;

        std::cout << "Item d) " << std::endl;
        std::cout << LConst((a^b),(e(1)),ortonormal) << std::endl << std::endl;

        std::cout << "Item e) " << std::endl;
        std::cout << SCP((2*a + b),(a+b),ortonormal) << std::endl << std::endl;

        std::cout << "Item f) " << std::endl;
        std::cout << RConst((e(1)^e(2)^e(3)),b,ortonormal) << std::endl << std::endl;

     // Exercício 5
            // dúvida

     // Exercício 6
            // dúvida

     // Exercício 7

		std::cout << std::endl;
		std::cout << "Question 7." << std::endl;

		BASE_BLADE bladesBase = {(e(0)),(e(1)),(e(2)),(e(3)), (e(1)^e(2)), (e(1)^e(3)), (e(2)^e(3)),(e(1)^e(2)^e(3))};

		CLIFFORD_TABLE tableGP;
		createCliffordTable(bladesBase, OPERATION::GP, tableGP, &ortonormal);
		std::cout << std::endl;
		printTableBladeBase("Geometric product of blade!", bladesBase ,tableGP);

		CLIFFORD_TABLE tableSCP;
		createCliffordTable(bladesBase, OPERATION::SCP, tableSCP, &ortonormal);
		std::cout << std::endl;
		printTableBladeBase("Scalar product of blade!", bladesBase, tableSCP);

		CLIFFORD_TABLE tableLCONST;
		createCliffordTable(bladesBase, OPERATION::LCONST, tableLCONST, &ortonormal );
		std::cout << std::endl;
		printTableBladeBase("Left contraction of blade!", bladesBase, tableLCONST);

		CLIFFORD_TABLE tableOUTP;
		createCliffordTable(bladesBase, OPERATION::OUTP, tableOUTP);
		std::cout << std::endl;
		printTableBladeBase("Outer product!", bladesBase, tableOUTP);


    // Exercício 8

      // Não tem como executar esse exercício.





	return 0;
	
}

