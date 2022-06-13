// simulation for the lattice


#include <iostream>
#include <fstream>
#include <string>

#include "walk.h"
#include "../funzioni.h"
#include "../generatore/random.h"

using namespace std;

int main() {

   Random rnd;
   rnd.StartGen();

    auto error = [] (vector<double> AV, vector<double> AV2, int n) { 
        return (n == 0 ? 0 : sqrt((AV2[n] - (AV[n]*AV[n]))/double(n)));
    };

    vector<double> rndm(300);
    int N{100}; // num of blocks
    int M{100}; // num of walks for block
    int steps{100};

    walk value(steps);

    vector<double> sum(M) , sum2(M), err_prog(N), r2(M), r_blk(N), r_blk2(N);

    for(int l{}; l < N ; ++l) {

        r2.assign(M,0.);

        for(int k{}; k < M; ++k){ //create a block

            for(size_t i{}; i < rndm.size(); ++i)
                rndm[i] = rnd.Rannyu();
        
            value.DoRandom(rndm);

            for(size_t j{}; j < r2.size(); ++j){
                    r2[j] += value.GetRadius2(j)/double(M);

        
            }
        }

        for(int i{}; i < N; ++i) {
            r_blk[i] += sqrt(r2[i])/double(N);
            r_blk2[i] += r2[i]/double(N);
        }

        for(int i{}; i < N; ++i) 
            err_prog[i] = error(r_blk,r_blk2,i);

    }

    ofstream out("Discrete_RW.dat");

    for(int i{}; i < N; ++i) 
        out<<r_blk[i] <<" "<<err_prog[i]<<endl;


    return 0;
}