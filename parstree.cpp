/*
 * parstree.cpp
 *
 *  Created on: Nov 6, 2014
 *      Author: diep
 */

#include <cstring>
#include "parstree.h"
#include "tools.h"

ParsTree::ParsTree(): IQTree(){
    cost_matrix = NULL;
}

ParsTree::ParsTree(Alignment *alignment): IQTree(alignment){
    cost_matrix = NULL;
}

ParsTree::~ParsTree() {
    if(cost_matrix){
        delete cost_matrix;
        cost_matrix = NULL;
    }

    if (central_partial_pars)
        delete[] central_partial_pars;
    central_partial_pars = NULL;
}

void ParsTree::loadCostMatrixFile(char * file_name){
    if(cost_matrix){
        delete cost_matrix;
        cost_matrix = NULL;
    }
    if(strcmp(file_name, "fitch") == 0)
//    if(file_name == NULL)
    	cost_matrix = new SankoffCostMatrix(aln->num_states);
    else
    	cost_matrix = new SankoffCostMatrix(file_name);
}

/**
compute the tree parsimony score
@return parsimony score of the tree
*/
int ParsTree::computeParsimony(){
//	assert(root && root->isLeaf());
//    cout << "nstates = " << aln->num_states << endl;
//	clearAllPartialLH();
    assert(root->isLeaf());
    PhyloNeighbor *nei = ((PhyloNeighbor*) root->neighbors[0]);
    current_it = nei;
    assert(current_it);
    current_it_back = (PhyloNeighbor*) nei->node->findNeighbor(root);
    assert(current_it_back);

    int nptn = aln->size();
    if(_pattern_pars == NULL) _pattern_pars = aligned_alloc<BootValTypePars>(nptn + VCSIZE_USHORT);

    return computeParsimonyBranch((PhyloNeighbor*) root->neighbors[0], (PhyloNode*) root);
}

UINT computeCostFromChild(UINT child_cost, UINT transition_cost){
    return (child_cost == UINT_MAX) ? UINT_MAX : (child_cost + transition_cost);
}

/**
compute partial parsimony score of the subtree rooted at dad
@param dad_branch the branch leading to the subtree
@param dad its dad, used to direct the traversal
*/
void ParsTree::computePartialParsimony(PhyloNeighbor *dad_branch, PhyloNode *dad){
	// don't recompute the parsimony
    if (dad_branch->partial_lh_computed & 2)
        return;

    Node *node = dad_branch->node;
    //assert(node->degree() <= 3);
    int ptn;
    if(aln->num_states != cost_matrix->nstates){
        cout << "Your cost matrix is not compatible with the alignment"
            << " in terms of number of states. Please check!" << endl;
        exit(1);
    }
    int nstates = aln->num_states;
    assert(dad_branch->partial_pars);

    int pars_block_size = getParsBlockSize();

    if (node->isLeaf() && dad) {
//        cout << "############# leaf!" << endl;
        // external node
        for(int i = 0; i < pars_block_size - 1; i++)
            dad_branch->partial_pars[i] = UINT_MAX;
        dad_branch->partial_pars[pars_block_size - 1] = 0; // reserved for corresponding subtree pars
        for (ptn = 0; ptn < aln->size(); ptn++){
            // ignore const ptn because it does not affect pars score
            if (!aln->at(ptn).is_const) {
                int ptn_start_index = ptn * nstates;
                char state;
                if (node->name == ROOT_NAME) {
                    state = aln->STATE_UNKNOWN;
                } else {
                    assert(node->id < aln->getNSeq());
                    state = (aln->at(ptn))[node->id];
                }

                if (state < nstates) {
                    dad_branch->partial_pars[ptn_start_index + state] = 0;
                } else {
                    // unknown, ambiguous character
//                    cout << "####### ambigous state = " << int(state) << endl;
                    initLeafSiteParsForAmbiguousState(state, dad_branch->partial_pars + ptn_start_index);
                }
            }
        }
    } else {
//        cout << "############# internal!" << endl;
        // internal node
        UINT i, j, ptn, min_ptn_pars, min_child_ptn_pars;
        int ptn_start_index;
        UINT * partial_pars = dad_branch->partial_pars;
        for(int i = 0; i < pars_block_size; i++)
            partial_pars[i] = 0;

        FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME) {
            computePartialParsimony((PhyloNeighbor*) (*it), (PhyloNode*) node);
            UINT *partial_pars_child = ((PhyloNeighbor*) (*it))->partial_pars;
            for (ptn = 0; ptn < aln->size(); ptn++){
                // ignore const ptn because it does not affect pars score
                if (aln->at(ptn).is_const) continue;
                ptn_start_index = ptn * nstates;

                for(i = 0; i < nstates; i++){
                    // min(j->i) from child_branch
                    min_child_ptn_pars = computeCostFromChild(partial_pars_child[ptn_start_index], cost_matrix->columns[i][0]);
                    for(j = 1; j < nstates; j++)
                        if(computeCostFromChild(partial_pars_child[ptn_start_index + j], cost_matrix->columns[i][j]) < min_child_ptn_pars)
                            min_child_ptn_pars = computeCostFromChild(partial_pars_child[ptn_start_index + j], cost_matrix->columns[i][j]);
                    partial_pars[ptn_start_index + i] = computeCostFromChild(partial_pars[ptn_start_index + i], min_child_ptn_pars);
                }
            }
        }
        // calc subtree pars
        for (ptn = 0; ptn < aln->size(); ptn++){
            // ignore const ptn because it does not affect pars score
            if (aln->at(ptn).is_const) continue;
            ptn_start_index = ptn * nstates;
            min_ptn_pars = partial_pars[ptn_start_index];
            for(i = 1; i < nstates; i++){
                if(partial_pars[ptn_start_index + i] < min_ptn_pars)
                    min_ptn_pars = partial_pars[ptn_start_index + i];
            }
            partial_pars[pars_block_size - 1] += min_ptn_pars * aln->at(ptn).frequency;
        }
    }

    dad_branch->partial_lh_computed |= 2;
}

void ParsTree::initLeafSiteParsForAmbiguousState(char state, UINT * site_partial_pars){
    int i, nstates = aln->num_states;
    if(state < nstates) return; // no need for manipulate normal state

    if (state == aln->STATE_UNKNOWN){
        for(i = 0; i < nstates; i++) site_partial_pars[i] = 0;
        return;
    }

    if (state == STATE_INVALID){
        outError("Alignment contains invalid state. Please check your data!");
    }

    for(i = 0; i < nstates; i++) site_partial_pars[i] = UINT_MAX;

    switch (nstates) {
        case 2:
        	outError("Alignment contains invalid state. Please check your data!");
        	break;
        case 4: // DNA
			switch (state) {
				case 1+4+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('G')] = 0;
					return; // A or G, Purine
				case 2+8+3:
					site_partial_pars[aln->convertState('C')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return; // C or T, Pyrimidine
				case 1+8+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return; // A or T, Weak
				case 2+4+3:
					site_partial_pars[aln->convertState('G')] = 0;
					site_partial_pars[aln->convertState('C')] = 0;
					return; // G or C, Strong
				case 1+2+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('C')] = 0;
					return; // A or C, Amino
				case 4+8+3:
					site_partial_pars[aln->convertState('G')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return; // G or T, Keto
				case 2+4+8+3:
					site_partial_pars[aln->convertState('C')] = 0;
					site_partial_pars[aln->convertState('G')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return;// C or G or T
				case 1+2+8+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('C')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return; // A or C or T
				case 1+4+8+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('G')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return; // A or G or T
				case 1+2+4+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('G')] = 0;
					site_partial_pars[aln->convertState('C')] = 0;
					return; // A or G or C
				default:
					outError("Alignment contains invalid state. Please check your data!");
					return;
			}
			break;
        case 20: // Protein
        	if (state == 4+8+19){
        		site_partial_pars[aln->convertState('D')] = 0;
        		site_partial_pars[aln->convertState('N')] = 0;
        		return; // Aspartic acid (D) or Asparagine (N)
        	}
        	else if (state == 32+64+19){
        		site_partial_pars[aln->convertState('Q')] = 0;
        		site_partial_pars[aln->convertState('E')] = 0;
        		return; // Glutamine (Q) or Glutamic acid (E)
        	}
        	else{
        		outError("Alignment contains invalid state. Please check your data!");
        		return;
        	}
        default:
        	// unknown
        	outError("Alignment contains invalid state. Please check your data!");
        	return;
    }

}

/**
compute tree parsimony score based on a particular branch
@param dad_branch the branch leading to the subtree
@param dad its dad, used to direct the traversal
@param branch_subst (OUT) if not NULL, the number of substitutions on this branch
@return parsimony score of the tree
*/
int ParsTree::computeParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);

    if (!central_partial_pars)
        initializeAllPartialPars();

    // DTH: I don't really understand what this is for. ###########
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
//        cout << "swapped\n";
    }

    if ((dad_branch->partial_lh_computed & 2) == 0)
        computePartialParsimony(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 2) == 0)
        computePartialParsimony(node_branch, node);

    // now combine likelihood at the branch
    tree_pars = 0;
    int nstates = aln->num_states;
    UINT i, j, ptn, min_ptn_pars, min_node_ptn_pars, min_dad_ptn_pars;

    UINT *ptn_partial_pars = new UINT[nstates];
    if(!ptn_partial_pars){
        outError("Could not allocate for ptn_partial_pars\n");
        exit(1);
    }

    int ptn_start_index;
    for (ptn = 0; ptn < aln->size(); ptn++){
        _pattern_pars[ptn] = 0;
        if (aln->at(ptn).is_const) continue;
        ptn_start_index = ptn * nstates;
        for(i = 0; i < nstates; i++){
            // min(j->i) from node_branch
            min_node_ptn_pars = computeCostFromChild(node_branch->partial_pars[ptn_start_index], cost_matrix->columns[i][0]);
            for(j = 1; j < nstates; j++)
                if(computeCostFromChild(node_branch->partial_pars[ptn_start_index + j], cost_matrix->columns[i][j]) < min_node_ptn_pars)
                    min_node_ptn_pars = computeCostFromChild(node_branch->partial_pars[ptn_start_index + j], cost_matrix->columns[i][j]);
            // min(j->i) from dad_branch
            min_dad_ptn_pars = computeCostFromChild(dad_branch->partial_pars[ptn_start_index], cost_matrix->columns[i][0]);
            for(j = 1; j < nstates; j++)
                if(computeCostFromChild(dad_branch->partial_pars[ptn_start_index + j], cost_matrix->columns[i][j]) < min_dad_ptn_pars)
                    min_dad_ptn_pars = computeCostFromChild(dad_branch->partial_pars[ptn_start_index + j], cost_matrix->columns[i][j]);
            ptn_partial_pars[i] = computeCostFromChild(min_node_ptn_pars, min_dad_ptn_pars);
        }
        min_ptn_pars = ptn_partial_pars[0];
        for(i = 1; i < nstates; i++)
            if(min_ptn_pars > ptn_partial_pars[i])
                min_ptn_pars = ptn_partial_pars[i];
        _pattern_pars[ptn] = min_ptn_pars;
        tree_pars += min_ptn_pars * aln->at(ptn).frequency;
    }

    if (branch_subst)
        *branch_subst = tree_pars;
    if(ptn_partial_pars){
        delete [] ptn_partial_pars;
        ptn_partial_pars = NULL;
    }
    return tree_pars;
}

void ParsTree::initializeAllPartialPars() {
	PhyloTree::initializeAllPartialPars();
//    if(params->maximum_parsimony && (!_pattern_pars))
//    	_pattern_pars = new UINT[aln->size()];
//    int index = 0;
//    initializeAllPartialPars(index);
}

size_t ParsTree::getParsBlockSize(){
    // the extra one is reserved for subtree pars score
    // TODO minus const ptn
    return aln->size() * aln->num_states + 1;
}

void ParsTree::initializeAllPartialPars(int &index, PhyloNode *node, PhyloNode *dad) {
    size_t pars_block_size = getParsBlockSize();
    if (!node) {
        node = (PhyloNode*) root;
        // allocate the big central partial pars memory
        if (!central_partial_pars) {
            int memsize = (aln->getNSeq() - 1) * 4 * pars_block_size; // Important for handling informal NEXUS format (e.g. output of TNT)
            if (verbose_mode >= VB_MED)
                cout << "Allocating " << memsize * sizeof(UINT) << " bytes for partial parsimony vectors" << endl;
            central_partial_pars = new UINT[memsize];
            if (!central_partial_pars)
                outError("Not enough memory for partial parsimony vectors");
        }
        index = 0;
    }
    if (dad) {
        // make memory alignment 16
        // assign a region in central_partial_lh to both Neihgbors (dad->node, and node->dad)
        PhyloNeighbor *nei = (PhyloNeighbor*) node->findNeighbor(dad);
        nei->partial_pars = central_partial_pars + (index * pars_block_size);
        nei = (PhyloNeighbor*) dad->findNeighbor(node);
        nei->partial_pars = central_partial_pars + ((index + 1) * pars_block_size);
        index += 2;
        //assert(index < nodeNum * 2 - 1);
    }
    FOR_NEIGHBOR_IT(node, dad, it) initializeAllPartialPars(index, (PhyloNode*) (*it)->node, node);
}

UINT* ParsTree::newPartialPars(){
	UINT *ret = new UINT[getParsBlockSize()];
    return ret;
}

void ParsTree::initParsData(Params* pars_params) {
	if(!pars_params) return;
    if(cost_matrix == NULL) loadCostMatrixFile(pars_params->sankoff_cost_file);
}
