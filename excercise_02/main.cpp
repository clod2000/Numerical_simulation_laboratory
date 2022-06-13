
#include <iostream>
#include <fstream>
#include <string>
#include "../generatore/random.h"


#include <vector>
#include <cmath>


using namespace std;
 
int main (int argc, char *argv[]){

    Random rnd;
    rnd.StartGen();

/////exercise 2.1.



    auto func = [] (double x) {
        return M_PI/2.*cos(M_PI*x/2.) ;
    };

    auto p = [](double x) {
      
         return 2.*(1-x);
    }; 

    auto g = [&] (double x) {
        return func(x)/p(x);
    };

    auto inv = [&] (double y) {
        return 1-sqrt(1-(rnd.Rannyu()));
    };

    auto error = [] (vector<double> AV, vector<double> AV2, int n) { 
        return (n == 0 ? 0 : sqrt((AV2[n] - (AV[n]*AV[n]))/double(n)));
    };


    int M = 1000000;
    int N = 100;

    vector<double> ave(N), av2(N), sum_prog(N), su2_prog(N), err_prog(N);
    vector<double> aveT(N), av2T(N), sumT_prog(N), su2T_prog(N), errT_prog(N);
    for(auto i{0}; i < N; ++i) {

        double sum{};
        double sumT{};

        for(auto j{0};j<M/N; ++j) {
            sum += func(rnd.Rannyu());
            sumT += g(inv(rnd.Rannyu())); //distribution realized with inversion of the cumulative function
           // cerr<<sum<<endl;
        }

        ave[i] = sum/double(M/N);
        av2[i] = ave[i]*ave[i];

        aveT[i] = sumT/double(M/N);
        av2T[i] = aveT[i]*aveT[i];

    }

    ofstream out1("intUnif.dat");
    ofstream out2("intTayl.dat");

    for( int i{}; i < N; ++i) {
         for(auto j=0; j < i+1 ; ++j ) {
         sum_prog[i] += ave[j];
         su2_prog[i] += av2[j];

         sumT_prog[i] += aveT[j];
         su2T_prog[i] += av2T[j];
      }
      sum_prog[i]/=(i+1);
      out1 << sum_prog[i] << " ";

      su2_prog[i]/=(i+1);

      err_prog[i] = error(sum_prog,su2_prog,i);
      out1 << err_prog[i]<<endl;

      sumT_prog[i]/=(i+1);
      out2 << sumT_prog[i] << " ";

      su2T_prog[i]/=(i+1);

      errT_prog[i] = error(sumT_prog,su2T_prog,i);
      out2 << errT_prog[i]<<endl;
    }

  
return 0;
}