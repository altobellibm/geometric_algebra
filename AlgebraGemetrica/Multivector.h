
#include <map>
#include <type_traits>
#include <bitset>
#include <iostream>
#include <assert.h>


#include "Utils.h"
#include "Metric.h"


typedef unsigned int type;
template <typename T>
class Multivector{
public:
	Multivector(){}
	
    template<typename T1>
    friend Multivector<T1> e(type i);

    template<typename T1>
    friend bool is_blade(const Multivector<T1>& A);

    template<typename T1>
    friend Multivector<T1> grade_extraction(const Multivector<T1>& M, type);

    template<typename T1>
    friend std::ostream& operator<<(std::ostream& out, const Multivector<T1>& A);

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
    friend Multivector<typename std::common_type<T1, T2, T3>::type> GP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort);

	template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type> LConst(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort);

    template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type> RConst(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort);

    template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type> SCP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort);

    template<typename T1>
    friend Multivector<T1> Reverse(const Multivector<T1>& A);

    template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type> SQR_Norm_Reverse(const Multivector<T1>& A, const Orthogonal<T2>& ort);

    template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type> INV(const Multivector<T1>& A, const Orthogonal<T2>& ort);

    template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type> IGP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort);

private:
	typedef std::map<type, T> MAP;
	MAP m;
};

template<typename T1 = int>
Multivector<T1> e(type i){

    Multivector<T1> a;
    a.m[1 << (i - 1)] = 1;
	return a;
}

template<typename T1>
std::ostream& operator<<(std::ostream& out, const Multivector<T1>& A){

    std::string strMultivector;
    for(auto it = A.m.begin(); it != A.m.end(); it++){

        it->second > 0 ? strMultivector+="+" : strMultivector;

        strMultivector+= std::to_string(it->second) + "*";

        type masc = it->first;

        if(masc == 0)
            strMultivector+= "e(0)";
        else
        {
            int indice = 1;
            type mascr = 1;
            strMultivector+= "( ";

            while (masc != 0){

                if ( (masc & mascr) != 0) 
                      strMultivector+= "e(" + std::to_string(indice) + ") ^ ";
                
                indice++;
                masc = masc >> 1;
            }

			strMultivector = strMultivector.substr(0, strMultivector.size() - 2);
            strMultivector+= " ) ";
        }
    }

    strMultivector.size() == 0 ? strMultivector : strMultivector.substr(0, strMultivector.size()-1);
    out <<  strMultivector;
    return out;
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


template<typename T1>
bool is_blade(const Multivector<T1>& A){

    if (A.m.size() == 0)
        return false;

    auto it = A.m.begin();
    int grade = take_grade(it->first);

    for(it++; it != A.m.end(); it++){
        if(grade != take_grade(it->first))
            return false;
    }

    return true;
}

template<typename T1>
Multivector<T1> grade_extraction(const Multivector<T1>& mul, type _grade){
    Multivector<T1> R;

	for (auto it = mul.m.begin(); it != mul.m.end(); it++) 
		if (take_grade(it->first) == _grade)
			R.m[it->first] = it->second;

    return R;
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
Multivector<typename std::common_type<T1, T2, T3>::type> GP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){

     Multivector<typename std::common_type<T1, T2, T3>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
		for (auto itB = B.m.begin(); itB != B.m.end(); itB++){
            Multivector<typename std::common_type<T1, T2, T3>::type> D;
			D.m[itA->first^itB->first] = canonical_order(itA->first, itB->first)*metric_factor(itA->first & itB->first, ort)*itA->second*itB->second;
			C = C + D;
		}

	return C;
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> LConst(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
 
    Multivector<typename std::common_type<T1, T2, T3>::type> C = GP(A,B,ort);
    return grade_extraction(C, take_grade(B.m.begin()->first) - take_grade(A.m.begin()->first));
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> RConst(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
 
    Multivector<typename std::common_type<T1, T2, T3>::type> C = GP(A,B,ort);
    return grade_extraction(C, take_grade(A.m.begin()->first) - take_grade(B.m.begin()->first));
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> SCP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
 
    Multivector<typename std::common_type<T1, T2, T3>::type> C = GP(A,B,ort);
    return grade_extraction(C, 0);
}


template<typename T1>
Multivector<T1> Reverse(const Multivector<T1>& A){
  
    int r = take_grade(A.m.begin()->first);
    return ((-1)^( (r*(r-1)) >> 2))*A;
}



template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2>::type> SQR_Norm_Reverse(const Multivector<T1>& A, const Orthogonal<T2>& ort){
    return SCP(A, Reverse(A), ort);
}


template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> INV(const Multivector<T1>& A, const Orthogonal<T2>& ort){
  
    Multivector<typename std::common_type<T1, T2, T3>::type> S = SCP(A, Reverse(A), ort);
    Multivector<typename std::common_type<T1, T2, T3>::type> R = e(0)*(1.0/S.m.begin()->second);
    return GP(SQR_Norm_Reverse(R,ort),Reverse(A), ort);
}


template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> IGP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
  
    return GP(A,INV(B, ort), ort);
}

template <typename T>
Multivector<T> uminus(const Multivector<T>& A){
    return -1*A;
}
