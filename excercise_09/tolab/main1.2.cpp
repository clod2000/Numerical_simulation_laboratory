#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "cmath"

#include "../generatore/random.h"


using namespace std;

int main () {

   Random rnd;

   rnd.StartGen();

    vector<int> N {1,2,10,100};

    ofstream out_unif ("unif_1.2.dat");
    ofstream out_exp ("exp_1.2.dat");
    ofstream out_lor ("lor_1.2.dat");

    vector<double>sum_unif(N.size());

    for(auto i{0}; i< 10000; ++i){

        for(size_t j{}; j < N.size(); ++j){

            double sum_unif{};
            double sum_exp{};
            double sum_lor{};

            for(auto k{0}; k<N[j]; ++k){

               sum_unif+=rnd.Rannyu(0.,1.)/double(N[j]);
               sum_exp+=rnd.Exp(1.)/double(N[j]);
               sum_lor+=rnd.Lorentz(0.,1.)/double(N[j]);

            }
            out_unif<<sum_unif<<" ";
            out_exp<<sum_exp<<" ";
            out_lor<<sum_lor<<" ";    
        }
        out_unif<<endl;
        out_exp<<endl;
        out_lor<<endl;
    }

return 0;
}