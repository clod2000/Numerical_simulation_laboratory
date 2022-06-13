#pragma once
#include <iostream>
#include <cassert>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;


class position{

    public: 
        position() : m_pos(3), m_cont(3){
            
        }

        position(double x, double y, double z) : m_pos(3) , m_cont(3){
                m_pos[0] = x;
                m_pos[1] = y; 
                m_pos[2] = z;
        }

      
        void SetStepLenght(double a) { m_a = a;}
        double GetPos(int i) const { 
            assert(i>=0 && i<3);
            return m_pos[i];}

        void SetPos(int i, double x) {
            assert(i>=0 && i<3) ;
            m_pos[i] = x;}

        void AddPos(int i, double x) {
            assert(i>=0 && i<3) ;
            m_pos[i] += x;}

        void CopyPos( position p) { 
            for(int i{}; i<3; ++i) 
                m_pos[i] = p.GetPos(i) ;
            }
    
        double GetStepLenght() {return m_a;}
        double GetRadius2(){ return (m_pos[0]*m_pos[0] +m_pos[1]*m_pos[1] + m_pos[2]*m_pos[2]);}
        void SetRadius2(){m_r2 =(m_pos[0]*m_pos[0] + m_pos[1]*m_pos[1] + m_pos[2]*m_pos[2]); }
        
    protected:
        vector<double> m_pos;   // position for lattice x, y, z
        vector<double> m_cont; //continuum position, r,theta,phi
        double m_a = 1.;
        double m_r2; 

};


class walk : public position {

    public:

        walk() : m_step {100} , m_stride(m_step) {}
      
        walk(int step) : m_step {step} , m_stride(m_step) {}

        void DoRandom(std::vector<double> rnd) ; //generate a discrete random walk filling the vector stride
        void DoRandomCont(vector<double> rnd); //generate random walk in the continuum

        void SetStride(int i, position p) { 
            assert(i>=0 && i< int( m_stride.size() ) );
            m_stride[i].CopyPos(p);
           // std::cout<<"copy of "<<i<<" done"<<std::endl;
        }

        double GetRadius2(int i) {
            return m_stride[i].GetRadius2();
        }

        void SetRadius2(int i) {
            m_stride[i].SetRadius2();
        }

        double GetStep() {return m_step;}

        void Print() {
            for(position value : m_stride) {
                for(int i{}; i< 3; ++i)
                    std::cout<<value.GetPos(i)<<"\t";
                std::cout<<std::endl;
            }
        }
       
        walk& operator= ( const walk& w) {
            m_step = w.m_step;
            m_stride.resize(0); 

            for(position value : w.m_stride) {
                m_stride.push_back(value);
            }
            return *this;
        }
        
      
    private:

        int m_step;
        std::vector<position> m_stride;

};


void walk::DoRandom(std::vector<double> rnd)  { //vector with random numbers uniformly distributed in [0,1)

    position slave1, slave2;

    for(size_t i{}; i< m_stride.size(); ++i) {
        

        for(int i{}; i < 3; ++i)
            slave1.SetPos(i,0.);

        double rnd_step = 3*rnd[i]; 
        double rnd_sign = 1 + rnd[m_stride.size() + i]; //I take the number after  m_stride.size() in order to avoid correlations
        double sign = (rnd_sign > 1.5 ? 1 : -1);

        slave1.SetPos(int(rnd_step), sign*slave1.GetStepLenght() );
    
        for(size_t i{}; i < m_pos.size(); ++i )
            slave2.AddPos(i,slave1.GetPos(i));
        
        SetStride(i,slave2);
        GetRadius2(i); //charge the value of  r2 iin each walk position
    }
}


void walk::DoRandomCont(std::vector<double> rnd)  { 
    position slave1, slave2;

    for(size_t i{}; i< m_stride.size(); ++i) {

      
        for(int i{}; i < 3; ++i)
            slave1.SetPos(i,0.);

        double rnd_theta = 2*M_PI*rnd[i]; 
        double rnd_phi = acos( 1 -2*rnd[m_stride.size() + i] ); //Sampling of a solid angle
       
        slave1.SetPos(0 , m_a*sin(rnd_phi)*cos(rnd_theta)  ); //x
        slave1.SetPos(1 , m_a*sin(rnd_phi)*sin(rnd_theta)  ); //y
        slave1.SetPos(2 , m_a*cos(rnd_phi)  );                //z

        for(size_t i{}; i < m_pos.size(); ++i )
             slave2.AddPos(i,slave1.GetPos(i));


        SetStride(i,slave2);
        SetRadius2(i); //charge the value of  r2 in each walk position
    }
}