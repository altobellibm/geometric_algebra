
#include <map>
#include <type_traits>
#include <bitset>
#include <iostream>

#include "Utils.h"
#include "Metric.h"


typedef unsigned int type;
template <typename T>
class Multivector{
public:
	Multivector(){}
	
	friend Multivector<int> e(type i);

	template <typename T1, typename T2>
	friend Multivector<typename std::common_type<T1, T2>::type> operator+(const Multivector<T1>& A, const Multivector<T2>& B);

	template <typename T1, typename T2>
	friend Multivector<typename std::common_type<T1, T2>::type> operator-(const Multivector<T1>& A, const Multivector<T2>& B);

	template <typename T1, typename T2>
	friend Multivector<typename std::common_type<T1, T2>::type> operator*(const T1& l, const Multivector<T2>& A);

	template <typename T1, typename T2>
	friend Multivector<typename std::common_type<T1, T2>::type> operator*(const Multivector<T1>& A, const T2& l);

	template <typename T1, typename T2>
	friend Multivector<typename std::common_type<T1, T2>::type> operator^(const Multivector<T1>& A, const Multivector<T2>& B);

	template <typename T1, typename T2>
	friend Multivector<typename std::common_type<T1, T2>::type> RP(const Multivector<T1>& A, const Multivector<T2>& B, int dimention);

	template<typename T1, typename T2, typename T3>
	friend Multivector<typename std::common_type<T1, T2>::type> GP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort);

	template<typename T1, typename T2, typename T3>
	Multivector<typename std::common_type<T1, T2>::type> LConst(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort);


private:
	std::map<type, T> m;

};

Multivector<int> e(type i){
	Multivector<int> a;
	a.m[1 << (i - 1)] = 1;
	return a;
}

/*
std::bitset<8> a();
std::bitset<8> b();
std::bitset<8> c1();
std::bitset<8> c2();
std::cout << a << std::endl;
std::cout << b << std::endl;
std::cout << c1 << std::endl;
std::cout << c2 << std::endl;
*/

int canonical_order(type masc1, type masc2){

    type trocas = 0;
    masc1 = masc1 >> 1;
    while (masc1 != 0){
        trocas += Utils::HammingWeight(masc1 & masc2);
        masc1 = masc1 >> 1;
    }

    if ((trocas & 1) == 0)
        return +1;
    else
        return -1;
}

template<typename T>
T metric_factor(type masc, const Orthogonal<T>& orth){
	
	int indice = 0;
	T product_value = 1;
	type mascr = 1;

	do {

		if ((masc & mascr) != 0)
			product_value *= orth.eval(indice);

		masc = masc >> 1;
		indice++;

	} while (masc != 0);
		
	return product_value;
}

int take_grade(type masc){
    return Utils::HammingWeight(masc);
}

template<typename T> 
Multivector<T> grade_extraction(Multivector<T> mul, type){
	return mul;
}

 template <typename T1, typename T2>
 Multivector<typename std::common_type<T1, T2>::type> operator+(const Multivector<T1>& A, const Multivector<T2>& B){
	 Multivector<typename std::common_type<T1, T2>::type> C;

	 auto itA = A.m.begin();
	 auto itB = B.m.begin();

	 while (itA != A.m.end() && itB != B.m.end()){
		 if (itA->first < itB->first){
			 C.m[itA->first] = itA->second;
			 itA++;
		 }
		 else
			 if (itA->first > itB->first){
			 C.m[itB->first] = itB->second;
			 itB++;
			 }
			 else{
				 C.m[itB->first] = itA->second + itB->second;
				 itA++;
				 itB++;
			 }
	 }

	 while (itA != A.m.end()){
		 C.m[itA->first] = itA->second;
		 itA++;
	 }

	 while (itB != B.m.end()){
		 C.m[itB->first] = itB->second;
		 itB++;
	 }

	 return C;
 }

 template <typename T1, typename T2>
 Multivector<typename std::common_type<T1, T2>::type> operator-(const Multivector<T1>& A, const Multivector<T2>& B){
	 return A + uminus(B);
 }


template <typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> operator*(const T1& l, const Multivector<T2>& A){
	
	Multivector<typename std::common_type<T1, T2>::type> C;

    auto itA = A.m.begin();
	while (itA != A.m.end()){
        C.m[itA->first] = itA->second *l;
		itA++;
	}

	return C;
 }
 
template <typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> operator*(const Multivector<T1>& A, const T2& l){

	Multivector<typename std::common_type<T1, T2>::type> C;

	auto itA = A.m.begin();
    while(itA != A.m.end()){
        C.m[itA->first] = itA->second * l;
	}

	return C;
}

template <typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> operator^(const Multivector<T1>& A, const Multivector<T2>& B){

	Multivector<typename std::common_type<T1, T2>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
		for (auto itB = B.m.begin(); itB != B.m.end(); itB++){
            unsigned int zero = 0;
			if ((itA->first & itB->first) == zero){

				Multivector<typename std::common_type<T1, T2>::type> D;
				D.m[(itA->first | itB->first)] = canonical_order(itA->first, itB->first)*itA->second*itB->second;
				C = C + D;
			}
		}

	return C;
}


template <typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> RP(const Multivector<T1>& A, const Multivector<T2>& B, int dimention){

	Multivector<typename std::common_type<T1, T2>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
		for (auto itB = B.m.begin(); itB != B.m.end(); itB++){
				Multivector<typename std::common_type<T1, T2>::type> D;

				type mascr = itA->first & itB->first;
				if ((take_grade(itA->first) + take_grade(itB->first) - take_grade(mascr)) == dimention){					
					D.m[mascr] = canonical_order(itA->first^mascr, itB->first^mascr)*itA->second*itB->second;
				}
				C = C + D;
		}

	return C;
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2>::type> GP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
	Multivector<typename std::common_type<T1, T2>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
		for (auto itB = B.m.begin(); itB != B.m.end(); itB++){
			Multivector<typename std::common_type<T1, T2>::type> D;
			D.m[itA->first^itB->first] = canonical_order(itA->first, itB->first)*metric_factor(itA->first & itB->first, ort)*itA->second*itB->second;
			C = C + D;
		}

	return C;
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2>::type> LConst(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
	Multivector<typename std::common_type<T1, T2>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
		for (auto itB = B.m.begin(); itB != B.m.end(); itB++){
			Multivector<typename std::common_type<T1, T2>::type> D;
			D.m[itA->first^itB->first] = canonical_order(itA->first, itB->first)*metric_factor(itA->first & itB->first, ort)*itA->second*itB->second;
			C = C + grade_extraction(D, take_grade(itB->first) - take_grade(itA->first));
		}

	return C;
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2>::type> RConst(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
	Multivector<typename std::common_type<T1, T2>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
	for (auto itB = B.m.begin(); itB != B.m.end(); itB++){
		Multivector<typename std::common_type<T1, T2>::type> D;
		D.m[itA->first^itB->first] = canonical_order(itA->first, itB->first)*metric_factor(itA->first & itB->first, ort)*itA->second*itB->second;
		C = C + grade_extraction(D, take_grade(itA->first) - take_grade(itB->first));
	}

	return C;
}



template <typename T>
Multivector<T> uminus(const Multivector<T>& A){

	Multivector<T> B;

	auto itA = A.m.begin();
    while(itA != A.m.end()){
        B.m[itA->first] = -itA->second;
	}

	return B;
}
