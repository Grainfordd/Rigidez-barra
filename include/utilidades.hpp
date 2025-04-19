#ifndef UTILIDADES_H
#define UTILIDADES_H

#include <armadillo>
#include "Barra.hpp"

void mostrar_barras(Barra barras[], int num_elementos);
void escribir_resultado(vec datos, std::string dato);
arma::mat eliminar_fila(arma::mat A, int fila);
arma::mat eliminar_columna(arma::mat A, int columna);
arma::vec eliminar_elemento(arma::vec A, int indice);

#endif

