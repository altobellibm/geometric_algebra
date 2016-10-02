
#include <map>
#include <type_traits>
#include <bitset>
#include <iostream>
#include "Utils.h"

typedef unsigned int type;
template <typename T>
class Multivetor{
public:
	Multivetor(){}
	
	friend Multivetor<int> e(int i);

	template <typename T1, typename T2>
	friend Multivetor<typename std::common_type<T1, T2>::type> operator+(const Multivetor<T1>& A, const Multivetor<T2>& B);

	template <typename T1, typename T2>
	friend Multivetor<typename std::common_type<T1, T2>::type> operator-(const Multivetor<T1>& A, const Multivetor<T2>& B);

	template <typename T1, typename T2>
	friend Multivetor<typename std::common_type<T1, T2>::type> operator*(const T1& l, const Multivetor<T2>& A);

	template <typename T1, typename T2>
	friend Multivetor<typename std::common_type<T1, T2>::type> operator*(const Multivetor<T1>& A, const T2& l);

	template <typename T1, typename T2>
	friend Multivetor<typename std::common_type<T1, T2>::type> operator^(const Multivetor<T1>& A, const Multivetor<T2>& B);


private:
	std::map<type, T> m;

};

 Multivetor<int> e(int i){
	Multivetor<int> a;
	
	a.m[1 << i] = 1;
	return a;
}

 template <typename T1, typename T2>
 Multivetor<typename std::common_type<T1, T2>::type> operator+(const Multivetor<T1>& A, const Multivetor<T2>& B){
	 Multivetor<typename std::common_type<T1, T2>::type> C;

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
 friend Multivetor<typename std::common_type<T1, T2>::type> operator-(const Multivetor<T1>& A, const Multivetor<T2>& B){
	 return A + uminus(B);
 }


template <typename T1, typename T2>
Multivetor<typename std::common_type<T1, T2>::type> operator*(const T1& l, const Multivetor<T2>& A){
	
	Multivetor<typename std::common_type<T1, T2>::type> C;

	Multivetor<T1> B;
	B.m[0] = l;

	auto itA = A.m.begin();
	while (itA != A.m.end()){
		C.m[itA->first] = itA->second * B.m[0];
		itA++;
	}

	return C;
 }
 
template <typename T1, typename T2>
Multivetor<typename std::common_type<T1, T2>::type> operator*(const Multivetor<T1>& A, const T2& l){

	Multivetor<typename std::common_type<T1, T2>::type> C;

	Multivetor<T1> B;
	B.m[0] = l;

	auto itA = A.m.begin();
	While(itA != A.m.end()){
		C.m[itA->first] = itA->second * B.m[0];
	}

	return C;
}

template <typename T1, typename T2>
Multivetor<typename std::common_type<T1, T2>::type> operator^(const Multivetor<T1>& A, const Multivetor<T2>& B){

	Multivetor<typename std::common_type<T1, T2>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++){
		for (auto itB = B.m.begin(); itB != B.m.end(); itB++){
			int zero = 0;
			if ((itA->first & itB->first) == zero){
				Multivetor<typename std::common_type<T1, T2>::type> D;
				std::bitset<8> a(itA->first);
				std::bitset<8> b(itB->first);
				std::bitset<8> c1(itA->first | itB->first);
				std::bitset<8> c2(a | b);
				std::cout << a << std::endl;
				std::cout << b << std::endl;
				std::cout << c1 << std::endl;
				std::cout << c2 << std::endl;
				D.m[(itA->first | itB->first) >> 1] = canonical_order(itA->first, itB->first)*itA->second*itB->second;
				C = C + D;
			}
		}
	}

	return C;
}

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

template <typename T>
Multivetor<T> uminus(const Multivetor<T1>& A){

	Multivetor<T> B;

	auto itA = A.m.begin();
	While(itA != A.m.end()){
		C.m[itA->first] = - itA->second;
	}
	
	return B;
}