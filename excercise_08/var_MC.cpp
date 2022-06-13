#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "var_MC.h"
#include "../funzioni.h"


using namespace std;

int main() {


    Input();

    if(SA != 0) {
        SimulatedAnnealing();
        ReadOptimizedParameters();
        mu = mu_min;
        sigma = sigma_min;
        cout<<" new mu ="<<mu<<"\t new sigma = "<<sigma<<endl;
        }
    

    ofstream outx;
    outx.open("psi_distribution.dat");

  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; istep++)
    {
      Move(mu,sigma);
      outx<<xold<<endl;
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
 
  return 0;
}

double Gaussiana(double x, double x0, double sigma) {
    return exp(-pow( (x-x0)/sigma, 2.) /2.);
}

double Psi_trial(double x,double sigma, double mu) {
 return Gaussiana(x,mu,sigma) + Gaussiana(x,-mu,sigma);
}
double Psi2(double x,double sigma,double mu){
    return Psi_trial(x,sigma,mu)*Psi_trial(x,sigma,mu) ;
}
double Kinetic(double x, double sigma, double mu) {
    
    double a = (x-mu)*(x-mu)/(sigma*sigma);
    double b = (x+mu)*(x+mu)/(sigma*sigma);
    return 0.5/(sigma*sigma)*( 1. -( a*Gaussiana(x,mu,sigma) + b*Gaussiana(x,-mu,sigma) )/Psi_trial(x,sigma,mu) );
}

double Potenzial(double x) {
    return pow(x,4) -5./2.*pow(x,2);
}
double Hamiltonian(double x, double sigma,double mu) {
    return Potenzial(x) + Kinetic(x,sigma, mu);
}



void Input() {

    ifstream ReadInput;

    cout << "Quantum and variational MC simulation" << endl;
    cout << "The program uses hbar=1 and m=1 units " << endl;

    //Read input informations
    ReadInput.open("input.in");
    ReadInput >> SA; // if=0 standard simulation else Simulated Annealing 

    ReadInput >> temp;
    beta = 1./temp;
    if(SA) cout<<"performing Simulated Annealing with starting temperature = "<<temp<<endl;
    


    ReadInput >> mu;
    cout << "starting mu value = " << mu << endl;

    ReadInput >> sigma;
    cout << "starting sigma value = " << sigma << endl;

    ReadInput >> xold;
    cout << "Starting x position = " << xold << endl;

    ReadInput >> delta;
    ReadInput >> nblk;
    ReadInput >> nstep;
 
    //preparing for simulation
    rnd.StartGen();
    ie = 0;
    n_props = 1;

    ReadInput >> delta_opt;
    ReadInput >> delta_temp;
    ReadInput >> temp_rep;

   


}

void SimulatedAnnealing(void){

    double n_opt_steps = 1000;

    H_old = 0;
    double H_new = 0;
    double mu_new, sigma_new;
    double x = xold;
    int wd = 20;


    ofstream opt_out;
    opt_out.open("Results.dat");

   

    while(true){
        
        for(int k{}; k < temp_rep; ++k){

            H_old = 0;
            H_new = 0;

            xold = x;  //Reset of the starting position

            for(int i{}; i< n_opt_steps; ++i){
                Move(mu,sigma);
            H_old += Hamiltonian(xold,sigma,mu);
            }
            H_old /= n_opt_steps;

            mu_new =fabs( mu + delta_opt*(rnd.Rannyu() -0.5) );
            sigma_new =fabs( sigma + delta_opt*(rnd.Rannyu() - 0.5) );
            //in our situation the sign of mu and sigma is equivalent

            xold = x; //Reset of the starting position

            for(int i{}; i< n_opt_steps; ++i){
                Move(mu_new,sigma_new);
                H_new += Hamiltonian(xold,sigma_new,mu_new);
            }
            H_new /= n_opt_steps;
            
            //calculate H with mu and sigma, then evolve mu and sigma and recalculate H with the new parameters
            //then use Metropolis with Boltzmann weight, if accept, use mu and sigma as new parameters
            //start each time from x0
            //every few cycles lower the temperature 

            double p = exp( 1./temp*(H_old-H_new)) ;
            double a = rnd.Rannyu();
            //cerr<<"p , a   =   "<<p<<" , "<<a<<endl; 
            if(p >= a ){
                mu = mu_new;
                sigma = sigma_new;
                SA_accepted++;
                //opt_out<< temp<<setw(wd) << mu << setw(wd)<< sigma << setw(wd) << H_new <<endl;
            }
            SA_attempted++;

        }

        
        opt_out<< temp<<setw(wd) << mu << setw(wd)<< sigma << setw(wd) << H_new <<endl;
        
        temp -= delta_temp;
        //cout<<"actual temperature = Numerical"<<temp<<endl;

        if(temp <= delta_temp) break;
    }

}


void ReadOptimizedParameters(void) {
    ifstream in;
    double t, m , s, h, h_min;
    double k = righeNelFile("Results.dat");
    in.open("Results.dat");

    in>>t>>mu_min>>sigma_min>>h_min;

    for(int i{1}; i < k; ++i){
        in>>t>>m>>s>>h;
        if(h < h_min) {
            mu_min = m;
            sigma_min = s;
            h_min = h;
        }
    } 
    in.close();  
}




void Move(double mu, double sigma) {

    accepted = 0;
    attempted = 0;

    for(int i{}; i < nstep ; ++i) {
       
        double xnew;
        double p, psi2_old, psi2_new;

        attempted++;
        
        xnew =(xold + delta*(rnd.Rannyu() -0.5));
      
        psi2_old = Psi2(xold,sigma, mu);
        psi2_new = Psi2(xnew,sigma,mu);

        //Metropolis test
        p = min(1.,psi2_new/psi2_old);

        if( p >= rnd.Rannyu() ) {
            xold = xnew;
            accepted += 1;
        }
        
       // cerr<<"acceptance ratio: " << (double)accepted/(double)attempted<<endl;

       // cerr<<"---------------------------------------------------------------"<<endl;
    
    }
}

void Measure() {

    /*  //to check distribution is correct
    ofstream outx;
    outx.open("xfill.dat",ios::app);
    outx<<xold<<endl;
    outx.close();
    */
    double E_loc = Hamiltonian(xold,sigma,mu);
   /*
    cout<<"E_loc = "<<E_loc<<endl;
    cout<<"potenzial = "<<Potenzial(xold)<<endl;
    cout<<"kinetic = "<<Kinetic(xold,sigma,mu)<<endl;
    */
    walker[ie] = E_loc;

    return;
}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i){
     blk_av[i] += walker[i];
   }
   blk_norm ++;
}


void Averages(int iblk) //Print results for current block
{
    
    ofstream H_output;
    const int wd=12, wde=20, pde=10;
    
    //if(iblk%10000 == 0)
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
   // cout<< "accepted: "<<accepted<<'\t'<<"attemped: "<<attempted<<endl;
    
    H_output.open("output_H.dat",ios::app);
    
    stima_H = blk_av[ie]/blk_norm; //Total energy
    glob_av[ie] += stima_H;
    glob_av2[ie] += stima_H*stima_H;
    err_H=Error(glob_av[ie],glob_av2[ie],iblk);

    H_output << setw(wd) << iblk <<  setw(wde)<< setprecision(pde) << stima_H << setw(wde)<<setprecision(pde) << glob_av[ie]/(double)iblk << setw(wde) << setprecision(pde) << err_H << endl;

    cout << "----------------------------" << endl << endl;

    H_output.close();
   
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

