
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "../generatore/random.h"
#include "../blockfunctions.h"
#include "../funzioni.h"

using namespace std;

int main() {

    double S0 = 100;
    double T = 1;
    double K = 100;
    double r = 0.1;
    double sigma = 0.25;

    Random rnd;
    rnd.StartGen();


    double M{10000}; //steps
    double N{100}; //blocks
    double rndG;

    auto S = [] (double r, double z, double sigma, double t, double s0 ) { //Asset price function
        return s0*exp((r-sigma*sigma/2.)*(t) + sigma*z*sqrt(t));
    };

    vector<double> C_direct;
    vector<double> P_direct;
    double S_direct{};

//direct

    for(int i{}; i < M; ++i) {
        rndG = rnd.Gauss(0.,1.);
        S_direct = S(r,rndG,sigma,T,S0);
       // cerr<<S_direct<<endl;
        C_direct.push_back(exp(-r*T)*max(0.,S_direct -K));
        P_direct.push_back(exp(-r*T)*max(0.,K - S_direct));
    }

    PrintStatistic(C_direct,N,M,"Call_direct.dat");
    PrintStatistic(P_direct,N,M,"Put_direct.dat");

    double Call_mean = CalcolaMedia(C_direct);
    double Put_mean = CalcolaMedia(P_direct);


    cout<<"Call "<<Call_mean<<" "<<"Put_mean "<<Put_mean<<endl;

//discrete

    int steps{100};
    double deltat= T/double(steps);

    vector<double> C_discr;
    vector<double> P_discr;


    for(int i{}; i < M; ++i) {

    double oldS{S0}; 
    double newS;
    
    for(int j{0}; j<steps; ++j) {
        rndG = rnd.Gauss(0.,1.);
        newS = S(r,rndG,sigma,deltat,oldS);
        //cerr<<newS<<" "<<oldS<<endl;
        oldS = newS;
    }
    C_discr.push_back(exp(-r*T)*max(0.,newS -K));
    P_discr.push_back(exp(-r*T)*max(0.,K - newS));
    
    }

    PrintStatistic(C_discr,N,M,"Call_discr.dat");
    PrintStatistic(P_discr,N,M,"Put_discr.dat");




    return 0;
}