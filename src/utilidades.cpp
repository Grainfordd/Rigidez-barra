#include "../include/utilidades.hpp"
#include <iostream>

using namespace arma;
using namespace std;

void mostrar_barras(Barra barras[], int num_elementos) {
    cout << "********************************************************************************" << endl;
    for (int i = 0 ; i < num_elementos ; i++) {
        cout << "------------------------ Miembro " << i + 1 << " ----------------------------------------" << endl;
        barras[i].k.print();
        cout << "-------------------------------------------------------------------------" << endl;
    }
    cout << "********************************************************************************" << endl;
}

void escribir_resultado(vec datos, string dato){
	string nombre_archivo = "../res_" + dato + ".txt";
	ofstream archivo(nombre_archivo);

	if (archivo.is_open()){
		archivo << "Nodo \t x \t y\n";
		int j = 1;
		for (int i = 0; i < datos.n_elem; i++){
			if (i % 2 == 0) {
				archivo << (j) << "\t\t"<<datos(i) << "\t";
				j++;
			}
			else{
				archivo << datos(i) << "\n";
			}
		}

	}

}

mat eliminar_fila(mat A, int fila) {
    int largo = A.n_rows;
    if (fila == 0) return A.rows(fila + 1, largo - 1);
    if (fila + 1 == largo) return A.rows(0, fila - 1);
    return join_vert(A.rows(0, fila - 1), A.rows(fila + 1, largo - 1));
}

mat eliminar_columna(mat A, int columna) {
    int largo = A.n_cols;
    if (columna == 0) return A.cols(columna + 1, largo - 1);
    if (columna + 1 == largo) return A.cols(0, columna - 1);
    return join_horiz(A.cols(0, columna - 1), A.cols(columna + 1, largo - 1));
}

vec eliminar_elemento(vec A, int indice) {
    if (indice == 0) return A(span(1, A.n_elem - 1));
    if (indice == A.n_elem - 1) return A(span(0, A.n_elem - 2));
    return join_cols(A(span(0, indice - 1)), A(span(indice + 1, A.n_elem - 1)));
}

