
#include <iostream>
#include <fstream>
#include <string>
#include "../generatore/random.h"


#include <vector>
#include <cmath>


using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
  
////////exercise 1.3


   auto error = [] (vector<double> AV, vector<double> AV2, int n) { 
      return (n == 0 ? 0 : sqrt((AV2[n] - (AV[n]*AV[n]))/double(n)));
   }; 


   vector<double> needle(2);

   double L = 7; // stick lenght
   double d = 15;  // distance between lines

   int M = 1000000; //number of throws
   int Nu = 100;  //number of blocks
   int K = M/Nu ; //elements in each block


   vector<double> pi(Nu), pi2(Nu), error_prog(Nu), pi_prog(Nu), pi2_prog(Nu); //save the progressive sum of pi to get the mean

   ofstream piOut("pi.dat");

   double theta{};

   

   for(int j{}; j < Nu; ++j){
      double hit{}, tot{};

      for(int i=0; i < K; ++i) {
         needle[0] = rnd.Rannyu(0.,1000*d);   //create initial point of the stick
         // theta = rnd.Rannyu(0., 4 * atan(1.)); 
         //cerr<<theta<<endl;
         theta = rnd.RanAngle();
         needle[1] = needle[0] + cos(theta)*L; //create final point of the stick  considering its 
                                             // projection on the axis perpendicular to the lines 
         if(int(needle[0]/d) != int(needle[1]/d) )
            hit++;
         tot++;
      }

      double P= hit/double(tot);

      pi[j] = 2*L / (P*d);
      pi2[j] = pi[j]* pi[j];
   }

   for(auto i=0; i < Nu ; ++i ) {
      for(auto j=0; j < i+1 ; ++j ) {

      pi_prog[i] += pi[j];
      pi2_prog[i] += pi2[j];
      }

      pi_prog[i]/=(i+1);
      pi2_prog[i]/=(i+1);
      error_prog[i] = error(pi_prog,pi2_prog,i);

      piOut<<pi[i]<<" ";
      piOut<<pi_prog[i]<<" ";
      piOut<<error_prog[i]<<endl;
   }

   return 0;
}