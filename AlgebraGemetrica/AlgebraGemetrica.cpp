// AlgebraGemetrica.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include "Multivector.h"

typedef std::vector<Multivector<int>> BASE_BLADE;
typedef std::vector<std::vector<Multivector<int>>> CLIFFORD_TABLE;

enum class OPERATION { GP, SCP, LCONST, RCONST };

void createCliffordTable(const BASE_BLADE& _bladesBase, const Orthonormal& ortonormal, OPERATION _type, CLIFFORD_TABLE& _table){

   _table.resize(_bladesBase.size(), std::vector<Multivector<int>>(_bladesBase.size()));

   for(BASE_BLADE::size_type i = 0; i < _bladesBase.size(); i++)
      for(BASE_BLADE::size_type j = 0; j < _bladesBase.size(); j++){
          switch (_type) {
          case OPERATION::GP: _table[i][j] = GP(_bladesBase[i], _bladesBase[j],ortonormal); break;
          case OPERATION::SCP: _table[i][j] = SCP(_bladesBase[i], _bladesBase[j],ortonormal); break;
          case OPERATION::LCONST: _table[i][j] = LConst(_bladesBase[i], _bladesBase[j],ortonormal); break;
          case OPERATION::RCONST: _table[i][j] = RConst(_bladesBase[i], _bladesBase[j],ortonormal); break;
          default:
              break;
          }
      }
}

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

        Multivector<int> a = (e(1) + e(3));
        Multivector<int> b = (e(1) + e(2));

		std::cout << std::endl;
		std::cout << "Questão 4." << std::endl;

		std::cout << "Item a) " << std::endl;
        std::cout << SCP(a, b, ortonormal) << std::endl;

        std::cout << "Item b) " << std::endl;
        std::cout << LConst((e(3)),b, ortonormal) << std::endl;

        std::cout << "Item c) " << std::endl;
        std::cout << LConst((e(3)),(a^b),ortonormal) << std::endl;

        std::cout << "Item d) " << std::endl;
        std::cout << LConst((a^b),(e(1)),ortonormal) << std::endl;

        std::cout << "Item e) " << std::endl;
        std::cout << SCP((2*a + b),(a+b),ortonormal) << std::endl;

        std::cout << "Item f) " << std::endl;
        std::cout << RConst((e(1)^e(2)^e(3)),b,ortonormal) << std::endl;


     // Exercício 5
            // dúvida

     // Exercício 6
            // dúvida

     // Exercício 7

      BASE_BLADE bladesBase = {(e(0)),(e(1)),(e(2)),(e(3)), (e(1)^e(2)), (e(1)^e(3)), (e(2)^e(3)),(e(1)^e(2)^e(3))};

	  CLIFFORD_TABLE tableGP;
      createCliffordTable(bladesBase, ortonormal, OPERATION::GP, tableGP);

      CLIFFORD_TABLE tableSCP;
      createCliffordTable(bladesBase, ortonormal, OPERATION::SCP, tableSCP);

      CLIFFORD_TABLE tableLCONST;
      createCliffordTable(bladesBase, ortonormal, OPERATION::LCONST, tableLCONST);

      CLIFFORD_TABLE tableRCONST;
      createCliffordTable(bladesBase, ortonormal, OPERATION::RCONST, tableRCONST);
	  
    // Exercício 8

      // Não tem como executar esse exercício.





	return 0;
	
}

