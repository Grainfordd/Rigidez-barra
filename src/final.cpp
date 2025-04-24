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


    auto [nodos_info, puntos, lineas, nombres] = leerMallaMSH("../data/malla.msh");

	// Asignación de nodos 
	int num_nodos = nodos_info.n_rows;
	Nodo nodos[num_nodos];
	for (int i = 0; i < num_nodos ; i++){
		double x = nodos_info.row(i)[0];
		double y = nodos_info.row(i)[1];
		nodos[i] = Nodo(i+1, x, y);
	}

	// Asignación elementos
	int num_elementos = lineas.n_rows;
	Barra barras[num_elementos];
	for (int i = 0; i < num_elementos ; i++){
		Nodo nodo1 = nodos[((int)lineas.row(i)[1]) - 1];
		Nodo nodo2 = nodos[((int)lineas.row(i)[2]) - 1];

		barras[i] = Barra(nodo1, nodo2, 1, 1);
	}
	//--------------------------------------------------------------------------
	
	mostrar_barras(barras, num_elementos);

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

	/* mostrar_barras(barras, num_elementos); */
	for (int i = 0 ; i < num_elementos ; i++){
		K_global += barras[i].k_glob;
	}
	K_global.print();

	// ----------------------------------------------------------------------------------
	
	// ---------------------------- Reducción matriz global ----------------------------
	mat K_red = K_global;
	cout << "-----------------------------------------"<< endl;
	K_global.print();
	cout << "-----------------------------------------"<< endl;

	vec disp = ones(n); 
	int num_cond = puntos.n_rows;
	// Leo Las condiciones de borde
	for (int i = 0; i < num_cond ; i++){
		int tag = puntos.row(i)[1];
		int nodo = puntos.row(i)[2]; 

		if (nombres[tag] == "Fijo"){
			disp[nodo*2 - 2] = 0;
			disp[nodo*2 - 1] = 0;
		}
		if (nombres[tag] == "Rod_x") {
			disp[nodo*2 - 2] = 0;
		} 
		if (nombres[tag] == "Rod_y") {
			disp[nodo*2 - 1] = 0;
		} 
	}

	for (int i = n-1 ; i >= 0 ; i--){
		if (disp[i] == 0){
			K_red = eliminar_fila(K_red, i);
			K_red = eliminar_columna(K_red, i);
		}
	}
	// --------------------------------------------------------------------------------


	// ---------------------------- Calculo fuerzas y desplazamientos ------------------
	mat f_read;
	f_read.load(fuerzas_file);

	vec f_red = vectorise(f_read.t());
	vec disp_red;
	disp_red = solve(K_red, f_red);

	int k = disp_red.n_elem - 1;
	for (int i = disp.n_rows-1; i >= 0; i--){
		if (disp[i] != 0){
			disp[i] = disp_red[k];
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

