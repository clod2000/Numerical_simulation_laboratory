#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "../funzioni.h"
#include "../generatore/random.h"

#include "Population.h"
#include <armadillo>

using namespace std;
using namespace arma;

void Update_file(Population P);


int main() {

    Random rnd;
    rnd.StartGen();

    double r{1.};
    mat pos (2,ncity);
    double theta{};

    for(int i{}; i < ncity; ++i) {

        if(type_sim){
            pos(0,i) =rnd.Rannyu();
            pos(1,i) =rnd.Rannyu();
        }else{
            theta = rnd.RanAngle();
            pos(0,i) =r*cos(theta);
            pos(1,i) =r*sin(theta);
        }
    }
    
    pos.save("cities.dat", raw_ascii);


//starting_population

    Population Pop(pos,nind);
    Pop.Start_Pop();
    for(int i{}; i < nind ; ++i)
        Pop.Check_Bonds(r,i);
    Pop.SaveMatrix("Starting_Population.dat");



//first pop
   
// new population
    Population Pop2(Pop);
    Pop2.Try_Mutation();
    Pop2.Try_Crossover();

   // cout<< "Pop2 survived the crossover"<<endl;
    
    for(int i{}; i < nind; ++i){
        Pop2.Check_Bonds(r,i);
    }
   // Pop2.SaveMatrix("Pop2.dat");
    Pop2.Calculate_L();
    Pop2.Sort_for_L();
    //Pop2.SaveMatrix("Pop2Sorted.dat");
    Pop2.Calculate_L();
    //Pop2.Print_L();


    Update_file(Pop2);


    Population Pop3(Pop2,nind,rnd);   //create a population with nind individual selected by the selection algoritm
    //Pop3.SaveMatrix("Pop3.dat");
    Pop3.Sort_for_L();
  //  Pop3.SaveMatrix("Pop3Sorted.dat");
   // Pop3.Print_L();

    Update_file(Pop3);


    
    ///////////////////////////////////////// simulation

       for(int j{}; j< numgen; ++j){   
        cout<<"generation "<<j<<endl;
       
        Population P(Pop3,nind,rnd);
        //cout<<"new generation"<<endl;
        P.Try_Mutation();
        
        P.Try_Crossover();  // crossover is implemented in order to prefer the firsts elements of the matrix
        P.Sort_for_L();
        for(int i{}; i< nind; ++i)
           // P.Check_Bonds(r,i);

        P.Save_L_mean(j);
        Update_file(P);
        Pop3 = P;

        if( j == numgen -1) P.SaveMatrix("lastpop.dat");
    }

  



return 0;

}


void Update_file( Population P) {
    ofstream out_x, out_y;
    out_x.open("output_x.dat", ios::app );
    out_y.open("output_y.dat", ios::app );

    
    for ( int j = 0; j < P.Get_ncols(); ++j){
        out_x << setw(15) << P.Get_el(0, j)(0); 
        out_y << setw(15) << P.Get_el(0, j)(1);
    }
    out_x <<setw(15)<< P.Get_el(0, 0)(0) << endl;
    out_y <<setw(15)<< P.Get_el(0, 0)(1) << endl;
}


