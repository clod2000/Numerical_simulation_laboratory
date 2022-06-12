# pragma once

#include <iostream>
#include<cmath>
#include<cassert>

class Position {

    public:
   // Position(double r, double theta) : m_x{r*cos(theta)} , m_y{r*sin(theta)}{}
    Position() { m_x = 0.; m_y = 0.;}
    Position(double x, double y) : m_x{x}, m_y{y} {}

    double Get_x() { return m_x; }
    double Get_y() { return m_y;} 
    void Set_x(double x){ m_x = x;}
    void Set_y(double y){m_y = y;}
    void Print_xy();
    double GetDistance2( Position P) { return pow(m_x - P.Get_x() , 2 ) + pow(m_y - P.Get_y() , 2 ); }

    void checkbonds(double r);

    Position &operator=( Position &p){
        m_x = p.Get_x();
        return *this;
    }

    private:
    double m_x,m_y;

};

void Position:: Print_xy() { std::cout<<"("<<Get_x()<<","<<Get_y()<<")"<<std::endl;}

void Position::checkbonds(double r) {
    double epsilon = 1E-5;
    assert( (Get_x()*Get_x() + Get_y()*Get_y() - 1 < epsilon ) && "checkbonds not respected!");
}

void Swap(Position& p1, Position& p2) {
    Position t = p1;
    p1 = p2;
    p2 = t;
}

