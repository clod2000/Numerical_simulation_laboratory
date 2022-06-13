#pragma once

#include <armadillo>
#include <cassert>
#include <vector>
#include <list>
#include <string>
#include <iomanip>

#include <algorithm>

#include "../generatore/random.h"

// global variables

const int ncity = 50;
const int nind = 2000;  // numero di individui per generazione
const int numgen = 200; // numero di generazioni
const int type_sim = 0; // if 0 circle, else square

const double pm = 0.1;  // mutation probability
const double pc = 0.8;  // crossover probability
const double p_sel = 3 ; // esponente operatore di selezione

using namespace arma;
using namespace std;

// population is a matrix with  1 individuals for each row

class Population : public Random
{

public:
    // constructors
    Population(arma::mat pos, int nindividuals, Random rnd) : m_L(nindividuals) , m_rnd(rnd)
    {
      //  m_rnd.StartGen();
        m_pop.set_size(nindividuals, ncity);
        arma::rowvec slave(2);
        for (int i{}; i < ncity; ++i)
        {
            m_pop(0, i) = pos.col(i);
        }
    }
    Population(Population &P, Random rnd) : m_L(P.Get_nrows()) , m_rnd(rnd)
    {
       // m_rnd.StartGen();
        m_pop.set_size(P.Get_nrows(), P.Get_ncols());
        m_pop = P.Get_field_pop();
    }

    Population(Population &P, int num_rows, Random rnd) : m_L(num_rows), m_rnd(rnd)
    {
        // m_rnd.StartGen();
        // cout<<"creating new generation"<<endl;
        m_pop.set_size(num_rows, P.Get_ncols());

        Row_copy_from(0, 0, P); // copio sempre la prima riga

        //cerr<<"checkpoint "<<endl;
        for (int i{1}; i < m_pop.n_rows; ++i)
            Row_copy_from(i, P.selection_operator(), P);
    }

    // start function
    void StartGenerator() { m_rnd.StartGen(); }
    void SaveMatrix(std::string nmfile);
    void Start_Pop();
    void Check_Bonds(double r, int row);

    // axcess metods
    int Get_ncols() { return m_pop.n_cols; }
    int Get_nrows() { return m_pop.n_rows; }
    arma::field<mat> Get_field_pop() { return m_pop; }
    mat Get_el(int i, int j) { return m_pop(i, j); }
    void Set_el(int i, int j, int matel, double x) { m_pop(i,j)(matel) = x; }

    void Row_copy_from(int r1, int r2, Population &P)
    {
        assert(r1 < m_pop.n_rows);
        assert( r2 < m_pop.n_rows);
        for (int i{}; i < m_pop.n_cols; ++i)
            m_pop(r1, i) = P.Get_el(r2, i);
    }

    void Row_copy_to_vec(int row, vector<mat> &v)
    {
        assert(m_pop.n_cols == v.size());
        for (unsigned int i{}; i < m_pop.n_cols; ++i)
            v[i] = m_pop(row, i);
    }
    // for L
    double Get_Distance2(int row, int a, int b);
    void Calculate_L();
    double Get_L(int i) { return m_L[i]; }
    void Print_L()
    {
        for (size_t i{}; i < m_L.size(); ++i)
            cout << m_L[i] << endl;
    }
    double Get_L_mean()
    {
        double Lm{};
        for (size_t i{}; i < m_L.size() / 2.; ++i)
            Lm += m_L[i];
        return 2 * Lm / (double)m_L.size();
    }
    void Save_L_mean(int j, int myrnk);

    // mutations
    void Try_Mutation();
    void Mutation1(int nrow); // permutation of pairs
    void Mutation2(int nrow); // shift m elements for n steps starting from position a
    void Mutation3(int nrow); // permutation of m pairs
    void Mutation4(int nrow); // invert m elements

    // selection functions
    void Sort_for_L();
    int selection_operator();
    int selection_best_L();

    // crossover
    void Crossover(int row1, int row2);
    void Try_Crossover(void);

    // variables
private:
    arma::field<mat> m_pop;
    Random m_rnd;
    vector<double> m_L;
};

//////////////////////////////////// mutations ///////////////////////////

void Population::Mutation1(int nrow)
{
    // cout<<"Mutation 1"<<endl;
    unsigned int a = m_rnd.Rannyu(1., ncity);
    unsigned int b = m_rnd.Rannyu(1., ncity);
    //cout<< "a , b = "<<a<<","<<b<<endl;

    m_pop(nrow, a).swap(m_pop(nrow, b));
}

void Population::Mutation2(int nrow)

{
    // cout<<"Mutation 2"<<endl;
    int n = m_rnd.Rannyu(0., ncity - 1);
    int m = m_rnd.Rannyu(0., ncity - 2);
    int a = m_rnd.Rannyu(0., ncity - 2 - m);
    // cout<<"n , m , a "<<n<<" "<<m<<" "<<" "<<a<<endl;
    // cout<<"riga n " << nrow <<endl;

    vector<mat> vcopy(ncity);

    for (size_t i = 0; i < vcopy.size(); ++i)
        vcopy[i] = m_pop(nrow, i);

    vector<int> lcopy(ncity);

    for (size_t i = 0; i < lcopy.size(); ++i)
        lcopy[i] = i;

    lcopy.erase(lcopy.begin()); // tolgo il primo perchè deve restare fisso

    vector<int> index(m); // creo un vettore di indici che andranno ricollocati
    for (size_t i{}; i < index.size(); ++i)
        index[i] = (a + i + n) % (ncity - 1);

    vector<int> copiator(m); // vettore su cui copiare temporaneamente i dati
    for (size_t i{}; i < copiator.size(); ++i)
        copiator[i] = lcopy[a + i];

    for (int i{}; i < m; ++i) // elimino gli elementi che vanno da a ad a+m
        lcopy.erase(lcopy.begin() + a);

    while (lcopy.size() < ncity - 1)
    { // inserisco gli elementi eliminati nella posizione corretta
        for (int i{}; i < m; ++i)
        {
            if ((size_t)index[i] <= lcopy.size() && lcopy.size() < ncity - 1)
                lcopy.emplace(lcopy.begin() + index[i], copiator[i]);
        }
    }

    // cout<<"checkpoint"<<endl;

    for (int i{1}; i < ncity; ++i)
    {
        m_pop(nrow, i) = vcopy[lcopy[i - 1]];
        // cout<<i<<endl;
    }
}

void Population::Mutation3(int nrow)
{
    // cout<<"Mutation 3"<<endl;
    double med = (double)ncity / 2.;
    int m = m_rnd.Rannyu(1., med - 1);
    vector<mat> vcopy(ncity);
    for (int i = 0; i < ncity; ++i)
    {
        vcopy[i] = m_pop(nrow, i);
    }
    int a = m_rnd.Rannyu(1., med - m);      // genero la posizione di partenza degli m elementi
    int b = m_rnd.Rannyu(med, 2 * med - m); // genero la posizione di partenza degli altri m elementi
                                            //  cout<<"a = "<<a<<'\t'<<"b ="<<b<<endl;

    for (int i = a, j = b; i < a + m; ++i, ++j)
    {
        m_pop(nrow, j) = vcopy[i];
        m_pop(nrow, i) = vcopy[j];
    }
}

void Population::Mutation4(int nrow)
{
    // cout<<"Mutation 4"<<endl;
    int m = m_rnd.Rannyu(1., ncity - 2);
    vector<mat> vcopy(ncity);
    for (int i = 0; i < ncity; ++i)
    {
        vcopy[i] = m_pop(nrow, i);
    }
    int a = m_rnd.Rannyu(1., ncity - m);
    // cout<< a<<endl;
    for (int i{a}, j{a + m - 1}; i < a + m; ++i, --j)
    {
        m_pop(nrow, i) = vcopy[j];
    }
}

void Population::Try_Mutation()
{
    for (int nrow{1}; nrow < Get_nrows(); ++nrow)
    {
        if (m_rnd.Rannyu() <= pm)
            Mutation1(nrow);
        //cerr<<"mutation 1 "<<endl;
        if (m_rnd.Rannyu() <= pm)
            Mutation2(nrow);
        //cerr<<"mutation 2 "<<endl;
        if (m_rnd.Rannyu() <= pm)
            Mutation3(nrow);
        //cerr<<"mutation 3 "<<endl;
        if (m_rnd.Rannyu() <= pm)
            Mutation4(nrow);
        //cerr<<"mutation 4 "<<endl;
    }
}

/////////////////////////////// starts///////////////////////////////

void Population::Start_Pop()
{

    //cout << m_pop.n_rows << endl;
    for (size_t i{1}; i < m_pop.n_rows; ++i)
    {
        for (int j{0}; j < ncity; ++j)
            m_pop(i, j) = m_pop(i - 1, j);
        //cout<<"ok copy row "<< i<<endl;
        Mutation1(i);
        //cout<<"ok mutation row "<<i<<endl;
    }
}

void Population::Check_Bonds(double r, int row)
{
    double epsilon = 1E-10;
    double x, y = 0.;
    for (int i{}; i < ncity; ++i)
    {
        x = m_pop(row, i)(0);
        y = m_pop(row, i)(1);
    }
    if (type_sim)
        assert((x - r < epsilon) && (y - r < epsilon) && "checkbonds not respected!");
    else
        assert((x * x + y * y - r * r < epsilon) && "checkbonds not respected!");
}

void Population::SaveMatrix(std::string nmfile)
{

    ofstream out;
    out.open(nmfile);

    for (unsigned int i = 0; i < m_pop.n_rows; ++i)
    {
        for (unsigned int j = 0; j < m_pop.n_cols; ++j)
            out << setw(25) << m_pop(i, j)(0) << "," << m_pop(i, j)(1);
        out << endl;
    }

    out.close();
}

///////////////////// fitness function //////////////////

double Population::Get_Distance2(int row, int a, int b)
{
    assert(row < Get_nrows() );
    double x1 = m_pop(row, a)(0);
    double y1 = m_pop(row, a)(1);
    double x2 = m_pop(row, b)(0);
    double y2 = m_pop(row, b)(1);
    return pow(x1 - x2, 2) + pow(y1 - y2, 2);
}

void Population::Calculate_L()
{
    //cerr<<"into calculate_L..."<<endl;
    double L{};
    for (size_t j{}; j < m_L.size(); ++j)
    {
        L = 0.;
        for (int i{1}; i < ncity; ++i)
            L += Get_Distance2(j, i, i - 1);
        L += Get_Distance2(j, 0, ncity - 1);
        m_L[j] = L;
    }
}

void Population::Save_L_mean(int j, int myrnk)
{
    ofstream outL;
    outL.open("output.L_means." + to_string(myrnk), ios::app);
    outL << setw(15) << j + 1 << setw(15) << Get_L(0) << setw(15) << Get_L_mean() << endl;
}

/////////////// selection functions ////////////////////

int Population::selection_operator()
{
    int j= int(ncity * pow(m_rnd.Rannyu(), p_sel)) + 1;
    //cerr<<j<<endl;
    return j;
}
int Population::selection_best_L()
{
    double Lmin = m_L.back();
    int k = 0;
    for (size_t i{}; i < m_L.size(); ++i)
    {
        assert(m_L[i] != 0 && "L values not calculated yet!!");
        if (m_L[i] <= Lmin)
        {
            Lmin = m_L[i];
            k = i;
        }
    }
    return k;
}

void Population::Sort_for_L()
{
    // cout<<"start"<<endl;
    Calculate_L();
    //cerr<<"Calculate_L works"<<endl;
    arma::field<mat> fcopy(m_pop.n_rows, m_pop.n_cols); // creo field di copia
    // se il programma è lento fare una copia list
    vector<double> vcopy;
  
    
    fcopy = m_pop;
    
    //cerr<<"fcopy done" <<endl;

    int index;
    //cerr<<"starting for cycle"<<endl;
    for (unsigned int i{}; i < fcopy.n_rows; ++i)
    {  //cerco il valore minore di L
        index = selection_best_L();
        vcopy.push_back(m_L[index]);
         //cout<<index<<endl;

        for (unsigned int j{}; j < fcopy.n_cols; ++j) // metto l'individuo migliore in posizione i
            m_pop(i, j) = fcopy(index, j);
        m_L[index] += 1E5; // questo è per far si che al  prossimo giro del ciclo l'elemento migliore non venga più preso

    }
    m_L = vcopy; // riposiziono m_L nell'ordine corretto
    //cerr<<"exiting from Sort"<<endl;
}

///////////////////////// crossover ///////////////////////

void Population::Crossover(int rw1, int rw2)
{

    // cout<<"into crossover"<<endl;
    int row1 = rw1;
    int row2 = rw2;
    if(m_rnd.Rannyu()<0.5 || row1== 0)               //probabilità del 50% di salvare il padre
        row1 += Get_nrows() / 2;
    if(m_rnd.Rannyu()<0.5)     
        row2 += (rw2 < Get_nrows() / 2 ? Get_nrows() / 2 : 0);

    vector<mat> v1(ncity), v2(ncity);
    Row_copy_to_vec(row1, v1);
    Row_copy_to_vec(row2, v2);
    int cutter = m_rnd.Rannyu(1, ncity - 2);
    // cout<<"cutter "<<cutter<<endl;
    int cut_counter;

    // cout<<"copy vec ok" <<endl;

    for (size_t i{(size_t)cutter}, j{1}; i < v1.size() || j < v1.size(); ++j)
    { // aggiorno riga 2
        cut_counter = 1;
        // cout<<"cut_counter = "<<cut_coun ter<<endl;
        for (int k{1}; k < cutter; ++k)
            if (v1[j].front() != v2[k].front() && v1[j].back() != v2[k].back())
                cut_counter++;
        //  cout<<"cut counter 2 = "<<cut_counter<<endl;
        if (cut_counter == cutter)
        {

            m_pop(row2, i) = v1[j];
            ++i;
        }
        // cout<<"i , j  = "<< i<<","<<j<<endl;
    }

    // cout<<"crossover of v1 done"<<endl;

    for (size_t i{(size_t)cutter}, j{1}; i < v1.size() || j < v1.size(); ++j)
    { // aggiorno riga 2
        cut_counter = 1;
        // cout<<"cut_counter = "<<cut_coun ter<<endl;
        for (int k{1}; k < cutter; ++k)
            if (v2[j].front() != v1[k].front() && v2[j].back() != v1[k].back())
                cut_counter++;
        //  cout<<"cut counter 2 = "<<cut_counter<<endl;
        if (cut_counter == cutter)
        {

            m_pop(row1, i) = v2[j];
            ++i;
        }
    }
    // cout<<"crossover of v2 done"<<endl;
}

void Population ::Try_Crossover(void)
{

    for (int i{}; i < Get_nrows() / 2; ++i)
    {
        int j = selection_operator(); //  m_rnd.Rannyu(0,Get_nrows() -1);
        if (i != j && m_rnd.Rannyu() <= pc)
        {
            // cout<< "crossover tra "<< i <<" e "<<j <<endl;
            Crossover(i, j);
        }
    }
}
