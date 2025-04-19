#include "../include/Barra.hpp"

Barra::Barra(Nodo nodo1, Nodo nodo2, double E_in, double A_in) {
    xf = nodo2.x;
    yf = nodo2.y;
    xn = nodo1.x;
    yn = nodo1.y;
	E = E_in;
	A = A_in;
	EA = E*A;

    nodos[0] = nodo1.num_nodo;
    nodos[1] = nodo2.num_nodo;

    largo = sqrt(pow((xf - xn), 2) + pow((yf - yn), 2));
    cosx = (xf - xn) / largo;
    cosy = (yf - yn) / largo;

    T = {
        {cosx, cosy, 0, 0},
        {0, 0, cosx, cosy}
    };

    k_loc = {
        {1, -1},
        {-1, 1}
    };

    k_loc = E*A * k_loc / largo;
    k_glob = T.t() * k_loc * T;
}

double Barra::fuerza_interna(arma::vec desp){
	vec fuerza = k_loc * T * desp;
	return fuerza(1);
}

