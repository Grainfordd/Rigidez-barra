#include "../include/utilidades.hpp"
#include <iostream>
#include <string>

using namespace arma;
using namespace std;

void mostrar_barras(Barra barras[], int num_elementos) {
    cout << "********************************************************************************" << endl;
    for (int i = 0 ; i < num_elementos ; i++) {
        cout << "------------------------ Miembro " << i + 1 << " ----------------------------------------" << endl;
        barras[i].k_glob.print();
        cout << "-------------------------------------------------------------------------" << endl;
    }
    cout << "********************************************************************************" << endl;
}

void escribir_vtk(Nodo nodos[], Barra elementos[], vec disp, vec esf, string nombre_archivo, int num_nodos, int num_elem){
	 
	ofstream archivo("../" + nombre_archivo);

	if (archivo.is_open()){
		archivo << "# vtk DataFile Version 2.0\n"; // Línea que siempre tiene que tener hvk
		archivo << "Armadura\n";
		archivo << "ASCII\n";
		archivo << "DATASET UNSTRUCTURED_GRID\n\n";

		// -------------- Nodos ----------------------
		archivo << "POINTS " << num_nodos << " float\n";
		for (int i = 0; i < num_nodos ; i++){
			archivo << nodos[i].x << " " << nodos[i].y << " 0\n";
		}

		archivo << "\n";
		archivo << "CELLS " << num_elem << " " << num_elem*3 << "\n";
		for (int i = 0; i < num_elem ; i++) {
			archivo << "2 " << elementos[i].nodos[0]-1 << " " << elementos[i].nodos[1]-1 << "\n";
		}

		archivo << "\n";
		archivo << "CELL_TYPES " << num_elem << "\n";
		for (int i = 0 ;i < num_elem ; i++){
			archivo << "3\n";
		}

		archivo << "\n";
		archivo << "POINT_DATA " << num_nodos << "\n";
		archivo << "VECTORS Desplazamiento float\n";
		for (int i = 0; i < disp.n_elem ; i++){
			if (i%2 == 0){
				archivo << disp[i] << " ";
			}
			else {
				archivo << disp[i] << " 0\n";
			}
		}

		archivo << "\n\n";
		archivo << "CELL_DATA " << num_elem << "\n";
		archivo << "SCALARS Esfuerzo float 1\n";
		archivo << "LOOKUP_TABLE default\n";

		for (int i = 0; i < num_elem; i++) {
			archivo << esf(i) << "\n";
		}

	}

}
void escribir_esf(vec esf){
	ofstream archivo("../esfuerzos.txt");

	if (archivo.is_open()){
		archivo << "Elemento \t x\n";

		for (int i = 0 ; i < esf.n_elem ; i++){
			archivo << i+1 << "\t\t" << esf(i) << "\n";
		}
		archivo.close();

	}
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
		archivo.close();
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
