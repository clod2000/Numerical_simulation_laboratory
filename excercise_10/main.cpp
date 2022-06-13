#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "../funzioni.h"
#include "../generatore/random.h"
#include "Population.h"
#include <mpi.h>
#include <armadillo>

using namespace std;
using namespace arma;

void Update_file_MPI(Population P,int myrnk);

int main(int argc, char* argv[]) {

    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    Random rnd;
/*
    rnd.StartGen();
        
    mat pos (2,ncity);

    double r{1.};
    double theta{};
    int myrank =rank;
 
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

*/

    for(int i{}; i < size; ++i){
        if(rank == i)                     // Starting each core from a different seed
            rnd.StartGenMPI(rank +1);
    }
 

    ifstream inCities;
    inCities.open("American_capitals.dat");
    mat AmPos(2,50);
    for(int i{}; i < 50; ++i)
       inCities>>AmPos(0,i) >> AmPos(1,i);
    inCities.close();
    AmPos.save("Amcities.dat", raw_ascii);

   // if(rank == 0 ) pos.save("cities.dat", raw_ascii);
   

    Population Pop(AmPos,nind, rnd);
    Pop.Start_Pop();
    Pop.Sort_for_L();

  

    Population Pop2(Pop, rnd);
    Pop2.Try_Mutation();
    Pop2.Sort_for_L();
    Pop2.Try_Crossover();
    
    Pop2.Calculate_L();
    Pop2.Sort_for_L();
    //Pop2.Calculate_L();
    //Pop2.Print_L();


     



    vector<double> x_trans(ncity);
    vector<double> y_trans(ncity);


    for(int j{}; j< numgen; ++j){   
        //cerr<<"generation "<<j<<endl;

        Population P(Pop2,nind,rnd);
       // cerr<<"new generation"<<endl;
          
        P.Try_Mutation();
        
        P.Sort_for_L();
       // cerr<<"checkpoint L "<<endl;
        P.Try_Crossover();  
        //cerr<<"checkpoint Crossover "<<endl;
        P.Sort_for_L();
       
        P.Save_L_mean(j,rank);
        Update_file_MPI(P,rank);


   
    
        if( j % 30 == 0) {
            for(int nodo{}; nodo < size ; ++nodo){


                for(int i{}; i < ncity; ++i){   //creo un vettore da inviare
                    x_trans[i] = P.Get_el(0,i)(0);
                    y_trans[i] = P.Get_el(0,i)(1);
                    
                }
            /*
                cerr<<"node "<<nodo<<" prima:" <<endl;
                for(int i{0}; i < 10; ++i)
                    cerr<< x_trans[i] << '\t';
                cerr<<endl;
            */
                MPI_Bcast(&x_trans.front(), x_trans.size(), MPI_DOUBLE, nodo, MPI_COMM_WORLD);
                MPI_Bcast(&y_trans.front(), y_trans.size(), MPI_DOUBLE, nodo, MPI_COMM_WORLD);
            /*
                cerr<<"node "<< nodo <<"dopo:"<<endl;
                for(int i{0}; i < 10; ++i)
                    cerr<< x_trans[i] << '\t';
                cerr<<endl;
            */ 
               for(int i{}; i < ncity; ++i){
                    P.Set_el(P.Get_nrows() -1 -nodo, i, 0 ,x_trans[i] );
                    P.Set_el(P.Get_nrows() -1 -nodo, i, 1 ,y_trans[i] );
                }
                
            P.Sort_for_L();
            }
   
        }
    
        Pop2 = P;

        //if( j == numgen -1) P.SaveMatrix("lastpop.dat");
    }

    MPI_Finalize();

    return 0;





}








void Update_file_MPI( Population P, int myrnk) {
    ofstream out_x, out_y;
    out_x.open("output_x." + to_string(myrnk), ios::app );
    out_y.open("output_y." + to_string(myrnk), ios::app );

    
    for ( int j = 0; j < P.Get_ncols(); ++j){
        out_x << setw(15) << P.Get_el(0, j)(0); 
        out_y << setw(15) << P.Get_el(0, j)(1);
    }
    out_x <<setw(15)<< P.Get_el(0, 0)(0) << endl;
    out_y <<setw(15)<< P.Get_el(0, 0)(1) << endl;
}


