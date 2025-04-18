#ifndef BARRA_HPP
#define BARRA_HPP

#include <armadillo>
#include "Nodo.hpp"

using namespace arma;

class Barra {
public: 
    double xf, yf, xn, yn;
    double dem, largo;
    double cosx, cosy;
    int nodos[2];
    mat T, k, k_ex;

    Barra(Nodo nodo1 = Nodo(), Nodo nodo2 = Nodo());
};

#endif

