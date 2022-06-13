#pragma once

#include "Position.h"
#include <vector>

class Individual: public Position {

    public:

    Individual() : m_pos(34) {}

    void Mutation1(int a , int b); //permutation of pairs

    private:
    std::vector<Position> m_pos;

};


void Individual::Mutation1(int a, int b) {
    if( a < m_pos.size() && b < m_pos.size() ){
        Swap(m_pos);
    }
}
