#include <iostream> 
#include <fstream>

using namespace std;

int main() {

    int numblk = 50;

    ifstream inEn("output.ene.0");
    ifstream inHeat("output.heat.0");
    ifstream inMag("output.mag.0");
    ifstream inChi("output.chi.0");

    

    string h,c,m,x;
    for(int i{}; i < numblk ; ++i) {
        getline(inEn, h);
        getline(inHeat, c);
        getline(inMag, m);
        getline(inChi, x);
    }

    ofstream ofEn("ResEn.dat", ios::app);
    ofstream ofHeat("ResHeat.dat", ios::app);
    ofstream ofMag("ResMag.dat", ios::app);
    ofstream ofChi("ResChi.dat", ios::app);

    ofEn<<h<<endl;
    ofHeat<<c<<endl;
    ofMag<<m<<endl;
    ofChi<<x<<endl;


    return 0;


}