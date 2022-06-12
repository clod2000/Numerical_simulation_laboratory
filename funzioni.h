#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>

using namespace std;

template <typename T> double CalcolaMedia ( const vector<T> & v) {

     T accumulo = 0;

  if(v.size() == 0)
    return accumulo;

  for ( size_t k = 0 ; k < v.size() ; k++ ) { 
    accumulo = double(k) / double (k + 1)*accumulo + 1. / double (k + 1) * v[k];
  }
  
  return accumulo;
}

template <typename T> double CalcolaVarianza (const vector<T> & v) {

  T varianza = 0;
  
  if(v.size() == 0)
    return varianza;

  double average = 0;

  for ( size_t k = 0 ; k < v.size() ; k++ ) {

   double old_average = average;

    average = CalcolaMedia(v);

    varianza = 1. / (double)(k + 1) * (k * varianza + v[k] * v[k] + k * old_average * old_average) - average*average;
  }
  return varianza ;
  
}

template <typename T> double CalcolaMediana ( vector<T> v) {
  
  vector<T> w = v;
  
  T mediana = 0;

  sort(v.begin(),v.end());

  if( v.size() % 2 == 0 ) {
  mediana = (v[v.size() /2 -1] + v[v.size() /2 ] )/ (double) 2 ;
  } else {
    mediana = v[v.size() /2 ];
  }

  return mediana;
} 

template <typename T> void print(const vector<T>& v) {
    for(size_t i = 0; i < v.size() ; i++)
      if(v[i] != 0)  
        cout <<v[i] << endl;
}

template <typename T> void print(const char * Filename , const vector<T>& v) {
     ofstream out(Filename);

    out <<"Dati riordinati: "<<endl;

    for(size_t i = 0 ; i < v.size() ; i++)
        out << v[i]<< endl;

} 

template <typename T> vector<T> Read(int ndata, const char* filename) {
    
    vector<T> v;
  
  ifstream in (filename);

  if(!in) {
    cerr<< "errore apertura del file"<<endl;
    exit(2);
  }
  
  for (int i = 0; i < ndata ; i++ ) {
    
    T value{};
    in >> value; 
    v.push_back(value);
    
  }
  return v;
}


bool are_close(double calculated, double expected, double epsilon = 1e-1) {
  return fabs(calculated - expected) < epsilon;
}

void test_functions() {  //test assert

    vector<double> data{1,2,3,4};
    vector<double> odd_data{1,2,3};
    
    assert(are_close(CalcolaMedia(data) , 2.5));
    assert(are_close(CalcolaVarianza(data) , 1.25));
    assert(are_close(CalcolaMediana(data) ,2.5));
    assert(are_close(CalcolaMediana(odd_data) , 2));
}



int righeNelFile(string file) {
    ifstream in;
    in.open( file );
    if ( in.fail() )
      return -1;
    string s;
    int k = 0;
    getline( in, s );
    while ( !in.eof() ) {
      k++;
      getline( in, s );
    }
    in.close();
    return k;
}