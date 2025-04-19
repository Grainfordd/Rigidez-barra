#include <iostream>
#include <armadillo>
#include "../include/Nodo.hpp"
#include "../include/Barra.hpp"
#include "../include/utilidades.hpp"

using namespace std;
using namespace arma;


int main(void){
	// ------------------------ Lectura de datos -------------------------------

	string nodos_file = "../data/nodos.txt";
	string elementos_file = "../data/elementos.txt";
	string fuerzas_file = "../data/fuerzas.txt";
	string disp_file = "../data/disp.txt";



	mat nodos_info;
	mat elementos_info;
	nodos_info.load(nodos_file);
	elementos_info.load(elementos_file);

	// Asignación de nodos 
	int num_nodos = nodos_info.n_rows;
	Nodo nodos[num_nodos];
	for (int i = 0; i < num_nodos ; i++){
		double x = nodos_info.row(i)[0];
		double y = nodos_info.row(i)[1];
		nodos[i] = Nodo(i+1, x, y);
	}

	// Asignación elementos
	int num_elementos = elementos_info.n_rows;
	Barra barras[num_elementos];
	for (int i = 0; i < num_elementos ; i++){
		Nodo nodo1 = nodos[((int)elementos_info.row(i)[0]) - 1];
		Nodo nodo2 = nodos[((int)elementos_info.row(i)[1]) - 1];

		double E = elementos_info.row(i)[2];
		double A = elementos_info.row(i)[3];
		barras[i] = Barra(nodo1, nodo2, E, A);
	}
	//--------------------------------------------------------------------------
	
	mostrar_barras(barras, 2);

	int n = num_nodos*2;


	mat K_global = zeros(n, n);


	// ---------------------------- Ensamblar matriz ----------------------------------
	for (int k = 0; k < num_elementos; k++){
		for (int i = 1; i<= n ; i++){
			bool cond1 = i != barras[k].nodos[0]*2 - 1 && i != barras[k].nodos[0]*2;
			bool cond2 = barras[k].nodos[1]*2 - 1 != i && i != barras[k].nodos[1]*2;

			bool cond = cond1 && cond2;
			if (cond){
				barras[k].k_glob.insert_rows(i-1, 1);
				barras[k].k_glob.insert_cols(i-1, 1);
			}
		}
	}

	for (int i = 0 ; i < num_elementos ; i++){
		K_global += barras[i].k_glob;
	}

	// ----------------------------------------------------------------------------------
	
	// ---------------------------- Reducción matriz global ----------------------------
	mat nodos_disp;
	nodos_disp.load(disp_file);

	vec disp = vectorise(nodos_disp.t());
	mat K_red = K_global;

	for (int i = num_nodos - 1 ; i >= 0 ; i--){
		int nodo = (i+1);

		if (nodos_disp(i, 1) == 0) {
			K_red = eliminar_fila(K_red, nodo*2 - 1);
			K_red = eliminar_columna(K_red, nodo*2 - 1);

		}
		if (nodos_disp(i, 0) == 0) {
			K_red = eliminar_fila(K_red, nodo*2 - 2);
			K_red = eliminar_columna(K_red, nodo*2 - 2);
		}
	}
	// --------------------------------------------------------------------------------


	// ---------------------------- Calculo fuerzas y desplazamientos ------------------
	mat f_read;
	f_read.load(fuerzas_file);

	vec f_red = vectorise(f_read.t());
	f_red.print();
	vec disp_red;
	disp_red = solve(K_red, f_red);

	int k = disp_red.n_elem - 1;
	for (int i = num_nodos - 1 ; i >= 0 ; i--){
		int nodo = (i+1);
		if (nodos_disp(i, 1) == 1) {
			disp[nodo*2 - 1] = disp_red[k];
			k--;
		}

		if (nodos_disp(i, 0) == 1) {
			disp[nodo*2 - 2] = disp_red[k];
			k--;
		}
	}

	mat F = K_global * disp;
	// --------------------------------------------------------------------------------
	// -------------------- Cálculo esfuerzos -----------------------------------------
	

	vec fuerza_int = zeros(num_elementos);
	vec esfuerzos = zeros(num_elementos);

	for (int i = 0; i < num_elementos ; i++){
		vec desp_loc = zeros(4);
		
		// Primer nodo
		int j = 2 * barras[i].nodos[0];

		desp_loc(0) = disp(j-2);
		desp_loc(1) = disp(j-1);

		// Segundo nodo

		j = 2 * barras[i].nodos[1];
		desp_loc(2) = disp(j-2);
		desp_loc(3) = disp(j-1);

		double fuerza_loc = barras[i].fuerza_interna(desp_loc);
		fuerza_int(i) = fuerza_loc;
	}

	for (int i = 0; i < num_elementos; i++){
		esfuerzos(i) = fuerza_int(i) / barras[i].A;
	}
	
	
	// --------------------------------------------------------------------------------


	cout << "-------- Fuerzas ----------" << endl;
	F.print();
	cout << "------------------" << endl;

	cout << "-------- Desplazamientos ----------" << endl;
	disp.print();
	cout << "------------------" << endl;

	cout << "-------- Fuerzas internas----------" << endl;
	fuerza_int.print();
	cout << "------------------" << endl;
	cout << "-------- esfuerzos  ---------------" << endl;
	esfuerzos.print();
	cout << "------------------" << endl;


	escribir_resultado(disp, "disp");
	escribir_resultado(F, "fuerzas");
	escribir_esf(esfuerzos);
	escribir_vtk(nodos, barras, disp, esfuerzos, "resultados.vtk", num_nodos, num_elementos);
}
