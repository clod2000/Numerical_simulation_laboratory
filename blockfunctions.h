#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

double error(vector<double> AV, vector<double> AV2, int n) { 
    return (n == 0 ? 0 : sqrt((AV2[n] - (AV[n]*AV[n]))/double(n)));
}

void PrintStatistic( vector<double> r , int numOfBlocks, int numOfEl, string outFile) {

    int N = numOfBlocks;
    int M = numOfEl;
    int L = M/N;

    vector<double> ave(N), av2(N), sum_prog(N), su2_prog(N),err_prog(N);

    for(auto i=0; i < N ; ++i ) {
      
        double sum{};

        for(auto j{0}; j < L ; ++j ) {
            int k = j + i*L;
            sum += r[k]; 
        }
        ave[i] = sum/double(L);
        av2[i] = ave[i]*ave[i];
    }

   ofstream out(outFile);
   
    for(auto i=0; i < N ; ++i ) {
        for(auto j=0; j < i+1 ; ++j ) {
            sum_prog[i] += ave[j];
            su2_prog[i] += av2[j];
        }
        sum_prog[i]/=(i+1);
        su2_prog[i]/=(i+1);
        err_prog[i] = error(sum_prog,su2_prog,i);

        out << sum_prog[i] <<" "<< err_prog[i]<<endl;
   }

   


}