

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "MD_MC.h"

using namespace std;

int main()
{ 
       
    counter=0;

    Input();
    double temptest = temp;

    cerr<<"temptest "<<temptest<<endl;

   
    double globaltemp=0;

    while(true){
        if(counter > 0) {
            if(globaltemp < temptest)
                temp+= temptest/( (double) counter +1.);
            else temp-= temptest/( (double) counter + 1.);
            cerr <<"glob av it: "<<globaltemp<<endl;
        }
        counter++;

    Input(); //Inizialization

        cerr<<endl<<"inizio ciclo numero "<<counter<<endl;
        
        for(int iblk=1; iblk <= nblk; iblk++) //Simulation
        {
            Reset(iblk);   //Reset block averages
            for(int istep=1; istep <= nstep; istep++)
            {
            Move();
            Measure();
            Accumulate(); //Update block averages
            
            }
            Averages(iblk);   //Print results for current block
            globaltemp = glob_av[it]/(double)iblk;
            
        } 
            if(fabs( globaltemp - temptest) < exp10(-4.))
                break;
    } 

    cerr<<endl<<"valore della temperatura iniziale da impostare:" <<temp<<endl;

  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;
/*
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;
*/
//Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

//Read input informations
  ReadInput.open("input.in");

  ReadInput >> iNVET;
  ReadInput >> restart;

  if(restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  if(counter == 0){
    ReadInput >> temp;
    beta = 1.0/temp;
    cout << "Temperature = " << temp << endl;

    ReadInput >> npart;
    //cout << "Number of particles = " << npart << endl;

    ReadInput >> rho;
    // cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    box = pow(vol,1.0/3.0);
    // cout << "Volume of the simulation box = " << vol << endl;
    // cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    // cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
        
    ReadInput >> delta;

    ReadInput >> nblk;

    ReadInput >> nstep;

    cout << "The program perform Metropolis moves with uniform translations" << endl;
    cout << "Moves parameter = " << delta << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();
  }
//Prepare arrays for measurements
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  ip = 4;
  n_props = 5; //Number of observables

//Read initial configuration
  cout << "Read initial configuration" << endl << endl;
 
    ReadConf.open("config.in");
    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i)
    {
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    cout << "velocity scale factor: " << fs << endl << endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  

  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();

  for (int i=0; i<npart; ++i)
  {
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    
  }
  
//Evaluate properties of the initial configuration
  Measure();

//Print initial values for measured properties

  cout << "Initial temperature      = " << walker[it] << endl;
  
  return;
}


void Move()
{

  double xnew, ynew, znew;

 
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){ //Verlet integration scheme

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;
      attempted = attempted + 1.0;
    }
  
  return;
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Measure() //Properties measurement
{
  double  kin=0.0;
      
  for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
 
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

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Epot, Ekin, Etot, Temp, Press;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
   // Temp.open("output_temp.dat",ios::app);
   
    stima_temp = blk_av[it]/blk_norm; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);

    

//Temperature
    cout << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
    cout << "----------------------------" << endl << endl;
  //  Temp.close();
  
}

/*
void ConfFinal(void)
{
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("config.out");
  WriteVelocity.open("velocity.out");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}
*/

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}
