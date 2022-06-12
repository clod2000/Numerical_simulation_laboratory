
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "../generatore/random.h"
#include "../blockfunctions.h"

using namespace std;

int main(int argc, char *argv[])
{

   Random rnd;
   rnd.StartGen();

   // es. 0.1.1.1

   int M = 100000, N = 100;

   vector<double> r; 
   for (auto i = 0; i < M; ++i)   
      r.push_back(rnd.Rannyu(0., 1.));

   PrintStatistic(r,N,M,"ex_1.1.1.dat");

   // es. 0.1.1.2

   r.resize(0);
   for (auto i = 0; i < M; ++i)
      r.push_back( pow(rnd.Rannyu() -0.5 , 2) ); 

   PrintStatistic(r,N,M,"ex_1.1.2.dat");


   // es.0.1.1.3

   M=100;
   int n= 10000;

   vector<double> rndm;
   for(auto i{0}; i < n*M; ++i) 
      rndm.push_back(rnd.Rannyu());

   auto chi = [] (double ni , int n, int M ) { //funzione per il calcolo del chi^2
      return (n == 0 ? 0 : ( pow(ni - n/M , 2 ) / ( n / M )) );
   };
  
   ofstream outChi("ex_1.1.3.dat");

   double chiVal{};

   vector<double> chi_progres(M);

   for(int k{}; k<100; ++k){
      chiVal = 0;
      for(int i{}; i < M; ++i) {
         int counter{};
         for(int j{}; j < n; ++j) {
            double rndmNum=rnd.Rannyu();
            if( (i)/double(M)<= rndmNum && rndmNum<=(i+1)/double(M) ) //verifico se rndmNum appartiene all'intervallo i-esimo
               counter++;
         } 
         chi_progres[i]=chi(counter,n,M);
         chiVal+= chi_progres[i];
      }
      outChi<<chiVal<<endl;
   }
   return 0;
}
