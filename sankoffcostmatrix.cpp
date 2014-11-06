/*
 * sankoffcostmatrix.cpp
 *
 *  Created on: Nov 6, 2014
 *      Author: diep
 */

#include "sankoffcostmatrix.h"
#include <iostream>
#include <fstream>
using namespace std;

SankoffCostMatrix::SankoffCostMatrix() {
    columns = NULL;
    nstates = 0;

}

// to work as Fitch cost
SankoffCostMatrix::SankoffCostMatrix(int num_states) {
    nstates = num_states;

    // allocate memory for columns
    columns = new UINT*[nstates];
    for(int i = 0; i < nstates; i++){
        columns[i] = new UINT[nstates];
    }

    // set every cells except for the diagonal to be 1
    for(int i = 0; i < nstates; i++){
        for(int j = 0; j < nstates; j++){
        	if(j == i)
        		columns[j][i] = 0;
        	else
        		columns[j][i] = 1;
        }
    }
}

SankoffCostMatrix::SankoffCostMatrix(char * cost_matrix_file){
    columns = NULL;
    nstates = 0;
    ifstream fin(cost_matrix_file);
    if(!fin.is_open()){
        cout << "ERROR while reading cost matrix file. Please check your input file!" << endl;
        exit(1);
    }
    fin >> nstates;

    // allocate memory for columns
    columns = new UINT*[nstates];
    for(int i = 0; i < nstates; i++){
        columns[i] = new UINT[nstates];
    }

    // read numbers from file
    for(int i = 0; i < nstates; i++){
        for(int j = 0; j < nstates; j++)
            fin >> columns[j][i];
    }

    fin.close();
}

SankoffCostMatrix::SankoffCostMatrix(const SankoffCostMatrix& existing_scm){
    nstates = existing_scm.nstates;
    // allocate memory for columns
    columns = new UINT*[nstates];
    for(int i = 0; i < nstates; i++){
        columns[i] = new UINT[nstates];
    }

    // copy numbers from the other matrix
    for(int i = 0; i < nstates; i++){
        for(int j = 0; j < nstates; j++)
            columns[i][j] = existing_scm.columns[i][j];
    }
}

SankoffCostMatrix::~SankoffCostMatrix() {
    for(int i = 0; i < nstates; i++){
        delete [] columns[i];
        columns[i] = NULL;
    }
    delete [] columns;
    columns = NULL;
}

void SankoffCostMatrix::print(ostream& out){
    if(!columns){
        cout << "Cost matrix is empty!" << endl;
        exit(1);
    }
    for(int i = 0; i < nstates; i++){
        for(int j = 0; j < nstates; j++)
            out << columns[j][i] << "\t";
        out << endl;
    }
}




