#include <iostream>
#include <fstream>
#include <armadillo>
#include <sstream>
#include <map>

using namespace std;
using namespace arma;

// tag → nombre
using PhysicalMap = map<int, string>;

/**
 * Lee un archivo .msh y extrae nodos, elementos y nombres físicos.
 * @param archivo Ruta del archivo .msh
 * @return Tupla con:
 *   - arma::mat nodos (x,y,z)
 *   - arma::imat elementos_punto (id, physical_tag, nodo)
 *   - arma::imat elementos_linea (id, nodo_i, nodo_j)
 *   - map<tag, nombre> physical_names
 */
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
                int nodo1 = datos[3 + num_tags - 2];
                int nodo2 = datos[3 + num_tags - 1];
                elementos_linea.insert_rows(elementos_linea.n_rows, imat{{id, nodo1, nodo2}});
            }
        }
    }

    file.close();
    return make_tuple(nodos, elementos_punto, elementos_linea, physical_names);
}

int main() {
    auto [nodos, puntos, lineas, nombres] = leerMallaMSH("../data/malla.msh");

    cout << "=== Physical Names ===\n";
    for (auto& [tag, nombre] : nombres) {
        cout << "Tag " << tag << " → " << nombre << endl;
    }

    cout << "\n=== Elementos Punto (con Physical Tag) ===\n";
    puntos.print();

    cout << "\n=== Elementos Línea (sin Tag) ===\n";
    lineas.print();

	cout << "\n=== Nodos " <<endl;
	nodos.print();
    return 0;
}


