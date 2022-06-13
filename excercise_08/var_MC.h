
#pragma once
#include "../generatore/random.h"

#include <vector>

// Random generator
Random rnd;

// parameters, observables
const int m_props = 1000;
int n_props, ie;
double walker[m_props];
int counter;
double delta_opt, delta_temp;
int temp_rep;
int SA_accepted, SA_attempted;

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_H;
double err_H;

// variables

double beta, temp, H_old, H_new;
double mu, sigma;
double xold;
double mu_min, sigma_min;

// simulation*/
int SA, nstep, nblk; // SA= simulated annealing
double delta;

// pigreco
const double pi = 3.1415927;

// functions
void Input(void);
void Accumulate(void);
void Reset(int);
void Averages(int);
void Move(double, double);
void Measure(void);
double Error(double, double, int);

double Potenzial(double);
double Gaussiana(double, double, double);
double Psi_trial(double, double, double);

void SimulatedAnnealing(void);
void ReadOptimizedParameters(void);