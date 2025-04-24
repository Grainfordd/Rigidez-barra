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


tuple<arma::mat, arma::imat, arma::imat, PhysicalMap> leerMallaMSH(const string& archivo) {
    ifstream file(archivo);
    arma::mat nodos;
    arma::imat elementos_punto, elementos_linea;
    PhysicalMap physical_names;

    if (!file.is_open()) {
        cerr << "Error al abrir: " << archivo << endl;
        return make_tuple(nodos, elementos_punto, elementos_linea, physical_names);
    }

    string linea;
    bool en_nodos = false, en_elementos = false;

    while (getline(file, linea)) {
        // Leer PhysicalNames
        if (linea == "$PhysicalNames") {
            int n;
            file >> n;
            file.ignore(); // Ignorar salto de línea
            for (int i = 0; i < n; ++i) {
                int dim, tag;
                string nombre_completo;
                file >> dim >> tag;
                getline(file, nombre_completo);
                size_t first_quote = nombre_completo.find("\"");
                size_t last_quote = nombre_completo.rfind("\"");
                string nombre = nombre_completo.substr(first_quote + 1, last_quote - first_quote - 1);
                physical_names[tag] = nombre;
            }
            continue;
        }

        if (linea == "$Nodes") {
            en_nodos = true;
            getline(file, linea);  // Leer número de nodos
            int num_nodos = stoi(linea);
            nodos.set_size(num_nodos, 3);
            continue;
        }
        if (linea == "$EndNodes") en_nodos = false;

        if (linea == "$Elements") {
            en_elementos = true;
            getline(file, linea);  // Leer número de elementos
            continue;
        }
        if (linea == "$EndElements") en_elementos = false;

        if (en_nodos) {
            istringstream iss(linea);
            int id;
            double x, y, z;
            iss >> id >> x >> y >> z;
            nodos(id - 1, 0) = x;
            nodos(id - 1, 1) = y;
            nodos(id - 1, 2) = z;
        }

        if (en_elementos) {
            istringstream iss(linea);
            vector<int> datos;
            int v;
            while (iss >> v) datos.push_back(v);

            int id = datos[0];
            int tipo = datos[1];
            int num_tags = datos[2];

            if (tipo == 15 && datos.size() >= 4 + num_tags - 1) {
                int physical_tag = (num_tags >= 1) ? datos[3] : -1;
                int nodo = datos[3 + num_tags - 1];
                elementos_punto.insert_rows(elementos_punto.n_rows, imat{{id, physical_tag, nodo}});
            }
            else if (tipo == 1 && datos.size() >= 4 + num_tags) {
                int nodo1 = datos[3 + num_tags];
                int nodo2 = datos[3 + num_tags + 1];
                elementos_linea.insert_rows(elementos_linea.n_rows, imat{{id, nodo1, nodo2}});
            }
        }
    }

    file.close();
    return make_tuple(nodos, elementos_punto, elementos_linea, physical_names);
}
