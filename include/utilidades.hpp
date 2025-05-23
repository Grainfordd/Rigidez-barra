#ifndef UTILIDADES_H
#define UTILIDADES_H

#include <armadillo>
#include "Barra.hpp"
/* #include <map> */

using PhysicalMap = std::map<int, std::string>;
void mostrar_barras(Barra barras[], int num_elementos);
void escribir_resultado(vec datos, std::string dato);
void escribir_vtk(Nodo nodos[], Barra elementos[], vec disp, vec esf, std::string nombre_archivo, int num_nodos, int num_elem);
void escribir_esf(vec esf);
arma::mat eliminar_fila(arma::mat A, int fila);
arma::mat eliminar_columna(arma::mat A, int columna);
arma::vec eliminar_elemento(arma::vec A, int indice);
/* std::tuple<arma::mat, arma::imat, arma::imat> leerMallaMSH(const std::string& archivo); */
std::tuple<arma::mat, arma::imat, arma::imat, PhysicalMap> leerMallaMSH(const std::string& archivo) ;

#endif

