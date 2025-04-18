#include "../include/Barra.hpp"

Barra::Barra(Nodo nodo1, Nodo nodo2) {
    xf = nodo2.x;
    yf = nodo2.y;
    xn = nodo1.x;
    yn = nodo1.y;

    nodos[0] = nodo1.num_nodo;
    nodos[1] = nodo2.num_nodo;

    largo = sqrt(pow((xf - xn), 2) + pow((yf - yn), 2));
    cosx = (xf - xn) / largo;
    cosy = (yf - yn) / largo;

    T = {
        {cosx, cosy, 0, 0},
        {0, 0, cosx, cosy}
    };

    k = {
        {1, -1},
        {-1, 1}
    };

    k = k / largo;
    k = T.t() * k * T;
}

