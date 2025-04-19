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
	double E, A, EA;
    int nodos[2];
    mat T, k_loc, k_glob;

    Barra(Nodo nodo1 = Nodo(), Nodo nodo2 = Nodo(), double E_in = 0, double A_in=0);
	
	double fuerza_interna(arma::vec desp);
};

#endif

