/*
 * sankoffcostmatrix.h
 *
 *  Created on: Nov 6, 2014
 *      Author: diep
 */

#include "tools.h"
#include <iostream>
using namespace std;

#ifndef SANKOFFCOSTMATRIX_H_
#define SANKOFFCOSTMATRIX_H_

class SankoffCostMatrix {
public:
    /**
     * default constructor
     */
    SankoffCostMatrix();

    /**
     * construct Fitch cost matrix
     * @param num_states: number of states in the alignment data
     */
    SankoffCostMatrix(int num_states);

    /**
     * construct matrix from file
     * @param cost_matrix_file: name of the file containing cost matrix
     * - first line is the matrix size
     * - following lines contain matrix cell values
     */
    SankoffCostMatrix(char * cost_matrix_file);

    /**
     * copy constructor
     */
    SankoffCostMatrix(const SankoffCostMatrix& existing_scm);

    virtual ~SankoffCostMatrix();

    /**
     * print to output device the Sankoff cost matrix
     */
    void print(ostream& out);

    int nstates; // number of states read from file
    UINT** columns; // cost read from file
};

#endif /* SANKOFFCOSTMATRIX_H_ */
