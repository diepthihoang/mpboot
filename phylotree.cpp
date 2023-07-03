//
// C++ Implementation: phylotree
//
// Description:
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "phylotree.h"
#include "bionj.h"
//#include "rateheterogeneity.h"
#include "alignmentpairwise.h"
#include <algorithm>
#include <limits>
#include "timeutil.h"
#include "nnisearch.h"
#include "phylosupertree.h"
#include "parstree.h"
#include "sprparsimony.h"
//const static int BINARY_SCALE = floor(log2(1/SCALING_THRESHOLD));
//const static double LOG_BINARY_SCALE = -(log(2) * BINARY_SCALE);

/****************************************************************************
 SPRMoves class
 ****************************************************************************/

void SPRMoves::add(PhyloNode* prune_node, PhyloNode* prune_dad, PhyloNode* regraft_node, PhyloNode* regraft_dad,
    double score) {
    if (size() >= MAX_SPR_MOVES && score <= rbegin()->score)
        return;
    if (size() >= MAX_SPR_MOVES) {
        iterator it = end();
        it--;
        erase(it);
    }
    SPRMove spr;
    spr.prune_node = prune_node;
    spr.prune_dad = prune_dad;
    spr.regraft_node = regraft_node;
    spr.regraft_dad = regraft_dad;
    spr.score = score;
    insert(spr);
}

/****************************************************************************
 PhyloTree class
 ****************************************************************************/

PhyloTree::PhyloTree() : MTree() {
    init();
}

void PhyloTree::init() {
    aln = NULL;
    model = NULL;
    site_rate = NULL;
    optimize_by_newton = true;
    central_partial_lh = NULL;
    tip_partial_lh = NULL;
    tip_partial_lh_computed = false;
    central_scale_num = NULL;
    central_partial_pars = NULL;
    model_factory = NULL;
    tmp_partial_lh1 = NULL;
    tmp_partial_lh2 = NULL;
    tmp_anscentral_state_prob1 = NULL;
    tmp_anscentral_state_prob2 = NULL;
    //tmp_ptn_rates = NULL;
    //state_freqs = NULL;
    tmp_scale_num1 = NULL;
    tmp_scale_num2 = NULL;
    discard_saturated_site = true;
    _pattern_lh = NULL;
    _pattern_lh_cat = NULL;
    _pattern_pars = NULL;
    //root_state = STATE_UNKNOWN;
    root_state = 126;
    theta_all = NULL;
    ptn_freq = NULL;
    ptn_invar = NULL;
    subTreeDistComputed = false;
    dist_matrix = NULL;
    sse = LK_SSE;  // FOR TUNG: you forgot to initialize this variable!
    save_all_trees = 0;
    mlCheck = 0; // FOR: upper bounds
    nodeBranchDists = NULL;
    save_all_trees = NULL;
    add_row = false;
}

PhyloTree::PhyloTree(Alignment* aln) : MTree() {
    init();
    this->aln = aln;  // Diep: This requires a call to setAlignment(aln); after the tree is inputted !!!!!!!!!!!!!!!!!!
}

void PhyloTree::discardSaturatedSite(bool val) {
    discard_saturated_site = val;
}

PhyloTree::~PhyloTree() {
    if (central_partial_lh)
        delete[] central_partial_lh;
    central_partial_lh = NULL;
    if (central_scale_num)
        delete[] central_scale_num;
    central_scale_num = NULL;

    if (central_partial_pars)
        delete[] central_partial_pars;
    central_partial_pars = NULL;
    if (model_factory)
        delete model_factory;
    if (model)
        delete model;
    if (site_rate)
        delete site_rate;
    if (tmp_scale_num1)
        delete[] tmp_scale_num1;
    if (tmp_scale_num2)
        delete[] tmp_scale_num2;
    if (tmp_partial_lh1)
        delete[] tmp_partial_lh1;
    if (tmp_partial_lh2)
        delete[] tmp_partial_lh2;
    if (tmp_anscentral_state_prob1)
        delete[] tmp_anscentral_state_prob1;
    if (tmp_anscentral_state_prob2)
        delete[] tmp_anscentral_state_prob2;
    //if (tmp_ptn_rates)
    //	delete [] tmp_ptn_rates;
    if (_pattern_lh_cat)
        delete[] _pattern_lh_cat;
    if (_pattern_lh)
        aligned_free(_pattern_lh);
    if (_pattern_pars) {
        aligned_free(_pattern_pars);
        _pattern_pars = NULL;
    }
    //if (state_freqs)
    //	delete [] state_freqs;
    if (theta_all)
        aligned_free(theta_all);
    if (ptn_freq)
        aligned_free(ptn_freq);
    if (ptn_invar)
        aligned_free(ptn_invar);
    if (dist_matrix)
        delete[] dist_matrix;
    
    // delete[] save_states_dad;
}

void PhyloTree::assignLeafNames(Node* node, Node* dad) {
    if (!node)
        node = root;
    if (node->isLeaf()) {
        node->id = atoi(node->name.c_str());
        assert(node->id >= 0 && node->id < leafNum);
        node->name = aln->getSeqName(node->id).c_str();
    }
    FOR_NEIGHBOR_IT(node, dad, it)assignLeafNames((*it)->node, node);
}

void PhyloTree::copyTree(MTree* tree) {
    MTree::copyTree(tree);
    if (!aln)
        return;
    // reset the ID with alignment
    setAlignment(aln);
}

void PhyloTree::copyTree(MTree* tree, string& taxa_set) {
    MTree::copyTree(tree, taxa_set);
    if (!aln)
        return;
    // reset the ID with alignment
    setAlignment(aln);
}

void PhyloTree::copyPhyloTree(PhyloTree* tree) {
    MTree::copyTree(tree);
    if (!tree->aln)
        return;
    setAlignment(tree->aln);
}

void PhyloTree::setAlignment(Alignment* alignment) {
    aln = alignment;
    //alnSize = aln->size();
    //ptn_freqs.resize(alnSize);
    //numStates = aln->num_states;
    //tranSize = numStates * numStates;
    //ptn_freqs.resize(alnSize);
    //for (int ptn = 0; ptn < alnSize; ++ptn) {
    //    ptn_freqs[ptn] = (*aln)[ptn].frequency;
    //}
    //block = aln->num_states * numCat;
    //lh_size = aln->size() * block;

    int nseq = aln->getNSeq();
    for (int seq = 0; seq < nseq; seq++) {
        string seq_name = aln->getSeqName(seq);
        Node* node = findLeafName(seq_name);
        if (!node) {
            string str = "Alignment has a sequence name ";
            str += seq_name;
            str += " which is not in the tree";
            outError(str);
        }
        assert(node->isLeaf());
        node->id = seq;
    }
}

void PhyloTree::setRootNode(char* my_root) {
    string root_name;
    if (my_root)
        root_name = my_root;
    else
        root_name = aln->getSeqName(0);
    root = findNodeName(root_name);
    assert(root);
}

void PhyloTree::readTreeString(const string& tree_string) {
    stringstream str;
    str << tree_string;
    str.seekg(0, ios::beg);
    freeNode();
    readTree(str, rooted);
    setAlignment(aln);
    if (isSuperTree()) {
        ((PhyloSuperTree*)this)->mapTrees();
    }
    else {
        clearAllPartialLH();
    }
}

void PhyloTree::readTreeFile(const string& file_name) {
    ifstream str;
    str.open(file_name.c_str());
    //	str << tree_string;
    //	str.seekg(0, ios::beg);
    freeNode();
    readTree(str, rooted);
    setAlignment(aln);
    if (isSuperTree()) {
        ((PhyloSuperTree*)this)->mapTrees();
    }
    else {
        clearAllPartialLH();
    }
    str.close();
}

string PhyloTree::getTreeString() {
    stringstream tree_stream;
    printTree(tree_stream);
    return tree_stream.str();
}

string PhyloTree::getTopology() {
    stringstream tree_stream;
    // important: to make topology string unique
    setRootNode(params->root);
    printTree(tree_stream, WT_TAXON_ID + WT_SORT_TAXA);
    return tree_stream.str();
}

void PhyloTree::rollBack(istream& best_tree_string) {
    best_tree_string.seekg(0, ios::beg);
    freeNode();
    readTree(best_tree_string, rooted);
    assignLeafNames();
    initializeAllPartialLh();
    clearAllPartialLH();
}

void PhyloTree::setModel(ModelSubst* amodel) {
    model = amodel;
    //state_freqs = new double[numStates];
    //model->getStateFrequency(state_freqs);
}

void PhyloTree::setModelFactory(ModelFactory* model_fac) {
    model_factory = model_fac;
}

void PhyloTree::setRate(RateHeterogeneity* rate) {
    site_rate = rate;
    if (!rate)
        return;
    //numCat = site_rate->getNRate();
    //if (aln) {
    //    block = aln->num_states * numCat;
    //    lh_size = aln->size() * block;
    //}
}

RateHeterogeneity* PhyloTree::getRate() {
    return site_rate;
}

Node* PhyloTree::newNode(int node_id, const char* node_name) {
    return (Node*)(new PhyloNode(node_id, node_name));
}

Node* PhyloTree::newNode(int node_id, int node_name) {
    return (Node*)(new PhyloNode(node_id, node_name));
}

void PhyloTree::clearAllPartialLH() {
    if (!root)
        return;
    ((PhyloNode*)root->neighbors[0]->node)->clearAllPartialLh((PhyloNode*)root);
    tip_partial_lh_computed = false;
}

void PhyloTree::computeAllPartialLh(PhyloNode* node, PhyloNode* dad) {
    if (!node) node = (PhyloNode*)root;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((((PhyloNeighbor*)*it)->partial_lh_computed & 1) == 0)
            computePartialLikelihood((PhyloNeighbor*)*it, node);
        PhyloNeighbor* rev = (PhyloNeighbor*)(*it)->node->findNeighbor(node);
        if ((rev->partial_lh_computed & 1) == 0)
            computePartialLikelihood(rev, (PhyloNode*)(*it)->node);
        computeAllPartialLh((PhyloNode*)(*it)->node, node);
    }
}

string PhyloTree::getModelName() {
    string name = model->name;
    if (model_factory->unobserved_ptns.size() > 0)
        name += "+ASC";
    name += site_rate->name;
    if (model->getFreqType() == FREQ_EMPIRICAL)
        name += "+F";
    return name;
}

string PhyloTree::getModelNameParams() {
    string name = model->getNameParams();
    if (model_factory->unobserved_ptns.size() > 0)
        name += "+ASC";
    name += site_rate->getNameParams();
    if (model->getFreqType() == FREQ_EMPIRICAL)
        name += "+F";
    return name;
}

/****************************************************************************
 Parsimony function
 ****************************************************************************/

 /*
  double PhyloTree::computeCorrectedParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad) {
  //	double corrected_bran = 0;
  //	int parbran;
  //	int parscore = computeParsimonyBranch(node21_it, node2, &parbran);
  //	if (site_rate->getGammaShape() != 0) {
  //		corrected_bran = (aln->num_states - 1.0) / aln->num_states
  //				* site_rate->getGammaShape()
  //				* (pow( 1.0 - aln->num_states / (aln->num_states - 1.0) * ((double) parbran / aln->getNSite()),
  //						-1.0 / site_rate->getGammaShape()) - 1.0);
  //	} else {
  //		corrected_bran = -((aln->num_states - 1.0) / aln->num_states)
  //				* log(1.0 - (aln->num_states / (aln->num_states - 1.0)) * ((double) parbran / aln->getNSite()));
  //	}
  //	return corrected_bran;
  }
  */

void PhyloTree::initializeAllPartialPars() {
    int index = 0;
    initializeAllPartialPars(index);
    //assert(index == (nodeNum - 1)*2);
}

void PhyloTree::initializeAllPartialPars(int& index, PhyloNode* node, PhyloNode* dad) {
    size_t pars_block_size = getBitsBlockSize();
    if (!node) {
        node = (PhyloNode*)root;
        // allocate the big central partial pars memory
        if (!central_partial_pars) {
            int memsize = (aln->getNSeq() - 1) * 4 * pars_block_size;
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
        PhyloNeighbor* nei = (PhyloNeighbor*)node->findNeighbor(dad);
        nei->partial_pars = central_partial_pars + (index * pars_block_size);
        nei = (PhyloNeighbor*)dad->findNeighbor(node);
        nei->partial_pars = central_partial_pars + ((index + 1) * pars_block_size);
        index += 2;
        //assert(index < nodeNum * 2 - 1);
    }
    FOR_NEIGHBOR_IT(node, dad, it)initializeAllPartialPars(index, (PhyloNode*)(*it)->node, node);
}

size_t PhyloTree::getBitsBlockSize() {
    // reserve the last entry for parsimony score
    // size: _pars
    // size - nptn - 1 : _pattern_pars[0]
    // size - 2: _pattern_pars[nptn-1]
    return (aln->num_states * aln->size() + UINT_BITS - 1) / UINT_BITS + aln->size() + 1;
}

int PhyloTree::getBitsEntrySize() {
    // reserve the last entry for parsimony score
    return (aln->num_states + UINT_BITS - 1) / UINT_BITS;
}

UINT* PhyloTree::newBitsBlock() {
    return new UINT[getBitsBlockSize()];
}

void PhyloTree::getBitsBlock(UINT* bit_vec, int index, UINT*& bits_entry) {
    int nstates = aln->num_states;
    int myindex = (index * nstates);
    int bit_pos_begin = myindex >> BITS_DIV;
    int bit_off_begin = myindex & BITS_MODULO;
    int bit_pos_end = (myindex + nstates) >> BITS_DIV;
    int bit_off_end = (myindex + nstates) & BITS_MODULO;

    if (bit_pos_begin == bit_pos_end) {
        bits_entry[0] = (bit_vec[bit_pos_begin] >> bit_off_begin) & ((1 << nstates) - 1);
        return;
    }
    UINT part1 = (bit_vec[bit_pos_begin] >> bit_off_begin);
    int rest_bits = nstates;
    int id;
    for (id = 0; rest_bits >= UINT_BITS; id++, rest_bits -= UINT_BITS, bit_pos_begin++) {
        bits_entry[id] = part1;
        if (bit_off_begin > 0)
            bits_entry[id] |= (bit_vec[bit_pos_begin + 1] << (UINT_BITS - bit_off_begin));
        part1 = (bit_vec[bit_pos_begin + 1] >> bit_off_begin);
    }
    if (bit_pos_begin == bit_pos_end) {
        bits_entry[id] = (bit_vec[bit_pos_begin] >> bit_off_begin) & ((1 << rest_bits) - 1);
        return;
    }
    UINT part2 = bit_vec[bit_pos_end];
    if (bit_off_end < UINT_BITS)
        part2 &= ((1 << bit_off_end) - 1);
    bits_entry[id] = part1;
    if (bit_off_begin > 0)
        bits_entry[id] |= (part2 << (UINT_BITS - bit_off_begin));
}

void PhyloTree::setBitsBlock(UINT*& bit_vec, int index, UINT* bits_entry) {
    int nstates = aln->num_states;
    int myindex = (index * nstates);
    int bit_pos_begin = myindex >> BITS_DIV;
    int bit_off_begin = myindex & BITS_MODULO;
    int bit_pos_end = (myindex + nstates) >> BITS_DIV;
    int bit_off_end = (myindex + nstates) & BITS_MODULO;

    //assert(value <= allstates);

    if (bit_pos_begin == bit_pos_end) {
        // first clear the bit between bit_off_begin and bit_off_end
        int allstates = (1 << nstates) - 1;
        bit_vec[bit_pos_begin] &= ~(allstates << bit_off_begin);
        // now set the bit
        bit_vec[bit_pos_begin] |= bits_entry[0] << bit_off_begin;
        return;
    }
    int len1 = UINT_BITS - bit_off_begin;
    // clear bit from bit_off_begin to UINT_BITS
    bit_vec[bit_pos_begin] &= (1 << bit_off_begin) - 1;
    // set bit  from bit_off_begin to UINT_BITS
    bit_vec[bit_pos_begin] |= (bits_entry[0] << bit_off_begin);
    int rest_bits = nstates - len1;
    int id;
    for (id = 0; rest_bits >= UINT_BITS; bit_pos_begin++, id++, rest_bits -= UINT_BITS) {
        bit_vec[bit_pos_begin + 1] = (bits_entry[id + 1] << bit_off_begin);
        if (len1 < UINT_BITS)
            bit_vec[bit_pos_begin + 1] |= (bits_entry[id] >> len1);
    }

    assert(bit_pos_begin == bit_pos_end - 1);
    // clear bit from 0 to bit_off_end
    bit_vec[bit_pos_end] &= ~((1 << bit_off_end) - 1);
    // now set the bit the value
    if (len1 < UINT_BITS)
        bit_vec[bit_pos_end] |= (bits_entry[id] >> len1);
    rest_bits -= bit_off_begin;
    if (rest_bits > 0)
        bit_vec[bit_pos_end] |= (bits_entry[id + 1] << bit_off_begin);
}

bool PhyloTree::isEmptyBitsEntry(UINT* bits_entry) {
    int rest_bits = aln->num_states;
    int i;
    for (i = 0; rest_bits >= UINT_BITS; rest_bits -= UINT_BITS, i++)
        if (bits_entry[i])
            return false;
    if (bits_entry[i] & ((1 << rest_bits) - 1))
        return false;
    return true;
}

void PhyloTree::unionBitsEntry(UINT* bits_entry1, UINT* bits_entry2, UINT*& bits_union) {
    int rest_bits = aln->num_states;
    int i;
    for (i = 0; rest_bits > 0; rest_bits -= UINT_BITS, i++)
        bits_union[i] = bits_entry1[i] | bits_entry2[i];
}

void PhyloTree::setBitsEntry(UINT*& bits_entry, int id) {
    int bit_pos = id >> BITS_DIV;
    int bit_off = id & BITS_MODULO;
    bits_entry[bit_pos] |= (1 << bit_off);
}

bool PhyloTree::getBitsEntry(UINT*& bits_entry, int id) {
    int bit_pos = id >> BITS_DIV;
    int bit_off = id & BITS_MODULO;
    if (bits_entry[bit_pos] & (1 << bit_off))
        return true;
    return false;
}

void setBitsAll(UINT*& bit_vec, int num) {
    //int id;
    int size = num / UINT_BITS;
    memset(bit_vec, 255, size * sizeof(UINT));
    num &= BITS_MODULO;
    if (num)
        bit_vec[size] = (1 << num) - 1;
}


/*
 * this will map DNA states to 0-15 (4-bit), i.e.,
 * 0 (A) -> 1
 * 1 (C) -> 2
 * 2 (G) -> 4
 * 3 (T) -> 8
 * A|C -> 3, etc.
 * A|C|G|T -> 15
 */
UINT dna_state_map[128];
UINT prot_state_map[128];
/*
 * this will recompute for 2 DNA states as input and map to result of Fitch algorithm:
 * For example, in bits:
 * 0100, 1110 -> 0100 (take intersection)
 * 1001, 0100 -> 1101 (take union because intersection = empty)
 * End effect: an array of size 256, with 4 left bit of index for state 1 and 4 right bit for state 2
 */
UINT dna_fitch_result[256];
/*
 * either 0 or 1
 */
UINT dna_fitch_step[256];

void precomputeFitchInfo() {
    dna_state_map[0] = 1; // A
    dna_state_map[1] = 2; // C
    dna_state_map[2] = 4; // G
    dna_state_map[3] = 8; // T
    dna_state_map[1 + 4 + 3] = 1 + 4; // A or G, Purine
    dna_state_map[2 + 8 + 3] = 2 + 8; // C or T, Pyrimidine
    dna_state_map[18] = 15; // STATE_UNKNOWN FOR DNA
    dna_state_map[1 + 8 + 3] = 1 + 8; // A or T, Weak
    dna_state_map[2 + 4 + 3] = 2 + 4; // G or C, Strong
    dna_state_map[1 + 2 + 3] = 1 + 2; // A or C, Amino
    dna_state_map[4 + 8 + 3] = 4 + 8; // G or T, Keto
    dna_state_map[2 + 4 + 8 + 3] = 2 + 4 + 8; // C or G or T
    dna_state_map[1 + 2 + 8 + 3] = 1 + 2 + 8; // A or C or T
    dna_state_map[1 + 4 + 8 + 3] = 1 + 4 + 8; // A or G or T
    dna_state_map[1 + 2 + 4 + 3] = 1 + 2 + 4; // A or G or C

    int state;
    for (state = 0; state < 256; state++) {
        UINT state1 = state & 15;
        UINT state2 = state >> 4;
        UINT intersection = state1 & state2;
        if (intersection == 0) {
            dna_fitch_result[state] = state1 | state2;
            dna_fitch_step[state] = 1;
        }
        else {
            dna_fitch_result[state] = intersection;
            dna_fitch_step[state] = 0;
        }
    }

    for (state = 0; state < 20; state++)
        prot_state_map[state] = (1 << state);
    prot_state_map[22] = (1 << 20) - 1; // STATE_UNKNOWN FOR PROTEIN
    prot_state_map[20] = 4 + 8; // N or D
    prot_state_map[21] = 32 + 64; // Q or E
}

/*
 * encode max. 8 consecutive parsimony state sites into max. 5 consecutive integers (32-bit)
 * @param encode (OUT) destination integer array (size <= 5)
 * @param sites source parsimony state site array (size <= 8)
 * @param num_sites number of sites
 */
void encodeProtState(UINT* encode, UINT* sites, int num_sites) {
    switch (num_sites) {
    case 8:
        encode[0] = sites[0] | (sites[1] << 20);
        encode[1] = (sites[1] >> 12) | (sites[2] << 8) | (sites[3] << 28);
        encode[2] = (sites[3] >> 4) | (sites[4] << 16);
        encode[3] = (sites[4] >> 16) | (sites[5] << 4) | (sites[6] << 24);
        encode[4] = (sites[6] >> 8) | (sites[7] << 12);
        break;
    case 1:
        encode[0] = sites[0];
        break;
    case 2:
        encode[0] = sites[0] | (sites[1] << 20);
        encode[1] = (sites[1] >> 12);
        break;
    case 3:
        encode[0] = sites[0] | (sites[1] << 20);
        encode[1] = (sites[1] >> 12) | (sites[2] << 8);
        break;
    case 4:
        encode[0] = sites[0] | (sites[1] << 20);
        encode[1] = (sites[1] >> 12) | (sites[2] << 8) | (sites[3] << 28);
        encode[2] = (sites[3] >> 4);
        break;
    case 5:
        encode[0] = sites[0] | (sites[1] << 20);
        encode[1] = (sites[1] >> 12) | (sites[2] << 8) | (sites[3] << 28);
        encode[2] = (sites[3] >> 4) | (sites[4] << 16);
        encode[3] = (sites[4] >> 16);
        break;
    case 6:
        encode[0] = sites[0] | (sites[1] << 20);
        encode[1] = (sites[1] >> 12) | (sites[2] << 8) | (sites[3] << 28);
        encode[2] = (sites[3] >> 4) | (sites[4] << 16);
        encode[3] = (sites[4] >> 16) | (sites[5] << 4);
        break;
    case 7:
        encode[0] = sites[0] | (sites[1] << 20);
        encode[1] = (sites[1] >> 12) | (sites[2] << 8) | (sites[3] << 28);
        encode[2] = (sites[3] >> 4) | (sites[4] << 16);
        encode[3] = (sites[4] >> 16) | (sites[5] << 4) | (sites[6] << 24);
        encode[4] = (sites[6] >> 8);
        break;
    default:
        outError(__func__);
        break;
    }
}

/*
 * decode max. 8 consecutive parsimony state sites from max. 5 consecutive integers (32-bit)
 * @param encode source integer array (size <= 5)
 * @param sites (OUT) parsimony state site array (size <= 8)
 * @param num_sites number of sites
 */
void decodeProtState(UINT* encode, UINT* sites, int num_sites) {
    sites[0] = encode[0] & ((1 << 20) - 1);
    if (num_sites == 8) {
        sites[1] = (encode[0] >> 20) | ((encode[1] & ((1 << 8) - 1)) << 12);
        sites[2] = (encode[1] >> 8) & ((1 << 20) - 1);
        sites[3] = (encode[1] >> 28) | ((encode[2] & ((1 << 16) - 1)) << 4);
        sites[4] = (encode[2] >> 16) | ((encode[3] & ((1 << 4) - 1)) << 16);
        sites[5] = (encode[3] >> 4) & ((1 << 20) - 1);
        sites[6] = (encode[3] >> 24) | ((encode[4] & ((1 << 12) - 1)) << 8);
        sites[7] = encode[4] >> 12;
        return;
    }
    if (num_sites >= 2)
        sites[1] = (encode[0] >> 20) | ((encode[1] & ((1 << 8) - 1)) << 12);
    if (num_sites >= 3)
        sites[2] = (encode[1] >> 8) & ((1 << 20) - 1);
    if (num_sites >= 4)
        sites[3] = (encode[1] >> 28) | ((encode[2] & ((1 << 16) - 1)) << 4);
    if (num_sites >= 5)
        sites[4] = (encode[2] >> 16) | ((encode[3] & ((1 << 4) - 1)) << 16);
    if (num_sites >= 6)
        sites[5] = (encode[3] >> 4) & ((1 << 20) - 1);
    if (num_sites >= 7)
        sites[6] = (encode[3] >> 24) | ((encode[4] & ((1 << 12) - 1)) << 8);
    if (num_sites >= 8)
        sites[7] = encode[4] >> 12;
}

void PhyloTree::computePartialParsimony(PhyloNeighbor* dad_branch, PhyloNode* dad) {
    // don't recompute the parsimony
    if (dad_branch->partial_lh_computed & 2)
        return;
    Node* node = dad_branch->node;
    //assert(node->degree() <= 3);
    int ptn;
    int nstates = aln->num_states;
    int pars_size = getBitsBlockSize();
    int entry_size = getBitsEntrySize();
    int nptn = aln->size();
    int ptn_pars_start_id = pars_size - nptn - 1;

    assert(dad_branch->partial_pars);
    dad_branch->partial_lh_computed |= 2;

    if (nstates == 4 && aln->seq_type == SEQ_DNA && (node->isLeaf() || node->degree() == 3)) {
        // ULTRAFAST VERSION FOR DNA, assuming that UINT is 32-bit integer
        if (node->isLeaf() && dad) {
            // external node
            for (ptn = 0; ptn < aln->size(); ptn += 8) {
                UINT states = 0;
                int maxi = aln->size() - ptn;
                if (maxi > 8) maxi = 8;
                for (int i = 0; i < maxi; i++) {
                    UINT bit_state = dna_state_map[(aln->at(ptn + i))[node->id]];
                    states |= (bit_state << (i * 4));
                    dad_branch->partial_pars[ptn_pars_start_id + ptn + i] = 0;
                }
                dad_branch->partial_pars[ptn / 8] = states;
            }
            //			// the remaining bits
            //			UINT states = 0;
            //			int maxi = aln->size() - ptn;
            //			for (int i = 0; i< maxi; i++) {
            //				UINT bit_state = dna_state_map[(aln->at(ptn+i))[node->id]];
            //				states |= (bit_state << (i*4));
            //			}
            //			dad_branch->partial_pars[(ptn/8)] = states;
            dad_branch->partial_pars[pars_size - 1] = 0; // set subtree score = 0
        }
        else {
            // internal node
            memset(dad_branch->partial_pars + ptn_pars_start_id, 0, nptn * sizeof(int));
            UINT* left = NULL, * right = NULL;
            int pars_steps = 0;
            FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME) {
                computePartialParsimony((PhyloNeighbor*)(*it), (PhyloNode*)node);
                if (!left)
                    left = ((PhyloNeighbor*)(*it))->partial_pars;
                else
                    right = ((PhyloNeighbor*)(*it))->partial_pars;
                pars_steps += ((PhyloNeighbor*)(*it))->partial_pars[pars_size - 1];
                for (int p = 0; p < nptn; p++)
                    dad_branch->partial_pars[ptn_pars_start_id + p] += ((PhyloNeighbor*)(*it))->partial_pars[ptn_pars_start_id + p];
            }
            for (ptn = 0; ptn < aln->size(); ptn += 8) {
                UINT states_left = left[ptn / 8];
                UINT states_right = right[ptn / 8];
                UINT states_dad = 0;
                int maxi = aln->size() - ptn;
                if (maxi > 8) maxi = 8;
                for (int i = 0; i < maxi; i++) {
                    UINT state_left = (states_left >> (i * 4)) & 15;
                    UINT state_right = (states_right >> (i * 4)) & 15;
                    UINT state_both = state_left | (state_right << 4);
                    states_dad |= dna_fitch_result[state_both] << (i * 4);
                    pars_steps += dna_fitch_step[state_both] * aln->at(ptn + i).frequency;
                    dad_branch->partial_pars[ptn_pars_start_id + ptn + i] += dna_fitch_step[state_both];
                }
                dad_branch->partial_pars[ptn / 8] = states_dad;
            }
            //            // remaining bits
            //			UINT states_left = left[ptn/8];
            //			UINT states_right = right[ptn/8];
            //			UINT states_dad = 0;
            //			int maxi = aln->size() - ptn;
            //			for (int i = 0; i< maxi; i++) {
            //				UINT state_left = (states_left >> (i*4)) & 15;
            //				UINT state_right = (states_right >> (i*4)) & 15;
            //				UINT state_both = state_left | (state_right << 4);
            //				states_dad |= dna_fitch_result[state_both] << (i*4);
            //				pars_steps += dna_fitch_step[state_both] * aln->at(ptn+i).frequency;
            //			}
            //			dad_branch->partial_pars[ptn/8] = states_dad;
            dad_branch->partial_pars[pars_size - 1] = pars_steps;
        }
        return;
    } // END OF DNA VERSION

    if (nstates == 20 && aln->seq_type == SEQ_PROTEIN && (node->isLeaf() || node->degree() == 3)) {
        // ULTRAFAST VERSION FOR protein, assuming that UINT is 32-bit integer
        if (node->isLeaf() && dad) {
            // external node
            UINT bit_ptn[8];
            int id = 0;
            for (ptn = 0; ptn < aln->size(); ptn += 8, id += 5) {
                int maxi = aln->size() - ptn;
                if (maxi > 8) maxi = 8;
                for (int i = 0; i < maxi; i++) {
                    bit_ptn[i] = prot_state_map[(aln->at(ptn + i))[node->id]];
                    dad_branch->partial_pars[ptn_pars_start_id + ptn + i] = 0;
                }
                encodeProtState(dad_branch->partial_pars + id, bit_ptn, maxi);
            }
            dad_branch->partial_pars[pars_size - 1] = 0; // set subtree score = 0
        }
        else {
            // internal node
            memset(dad_branch->partial_pars + ptn_pars_start_id, 0, nptn * sizeof(int));
            UINT* left = NULL, * right = NULL;
            int pars_steps = 0;
            FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME) {
                computePartialParsimony((PhyloNeighbor*)(*it), (PhyloNode*)node);
                if (!left)
                    left = ((PhyloNeighbor*)(*it))->partial_pars;
                else
                    right = ((PhyloNeighbor*)(*it))->partial_pars;
                pars_steps += ((PhyloNeighbor*)(*it))->partial_pars[pars_size - 1];
                for (int p = 0; p < nptn; p++)
                    dad_branch->partial_pars[ptn_pars_start_id + p] += ((PhyloNeighbor*)(*it))->partial_pars[ptn_pars_start_id + p];
            }
            UINT state_left[8], state_right[8], bit_ptn[8];
            int id = 0;
            for (ptn = 0; ptn < aln->size(); ptn += 8, id += 5) {
                int maxi = aln->size() - ptn;
                if (maxi > 8) maxi = 8;
                decodeProtState(left + id, state_left, maxi);
                decodeProtState(right + id, state_right, maxi);
                for (int i = 0; i < maxi; i++) {
                    bit_ptn[i] = state_left[i] & state_right[i];
                    if (!bit_ptn[i]) {
                        bit_ptn[i] = state_left[i] | state_right[i];
                        pars_steps += aln->at(ptn + i).frequency;
                        dad_branch->partial_pars[ptn_pars_start_id + ptn + i] += 1;
                    }
                }
                encodeProtState(dad_branch->partial_pars + id, bit_ptn, maxi);
            }

            dad_branch->partial_pars[pars_size - 1] = pars_steps;
        }
        return;
    } // END OF PROTEIN VERSION

    UINT* bits_entry = new UINT[entry_size];
    UINT* bits_entry_child = new UINT[entry_size];
    //UINT *bits_entry1 = new UINT[entry_size];
    //UINT *bits_entry2 = new UINT[entry_size];

    if (node->isLeaf() && dad) {
        // external node
        setBitsAll(dad_branch->partial_pars, nstates * aln->size());
        dad_branch->partial_pars[pars_size - 1] = 0;
        for (ptn = 0; ptn < aln->size(); ptn++)
            if (!aln->at(ptn).is_const) {
                char state;
                if (node->name == ROOT_NAME) {
                    state = aln->STATE_UNKNOWN;
                }
                else {
                    assert(node->id < aln->getNSeq());
                    state = (aln->at(ptn))[node->id];
                }
                if (state == aln->STATE_UNKNOWN) {
                    // fill all entries with bit 1
                    //setBitsBlock(dad_branch->partial_pars, ptn, (1 << nstates) - 1);
                }
                else if (state < nstates) {
                    memset(bits_entry, 0, sizeof(UINT) * entry_size);
                    setBitsEntry(bits_entry, state);
                    setBitsBlock(dad_branch->partial_pars, ptn, bits_entry);
                }
                else {
                    // ambiguous character, for DNA, RNA
                    state = state - (nstates - 1);
                    memset(bits_entry, 0, sizeof(UINT) * entry_size);
                    bits_entry[0] = state;
                    setBitsBlock(dad_branch->partial_pars, ptn, bits_entry);
                }
                dad_branch->partial_pars[ptn_pars_start_id + ptn] = 0;
            }
    }
    else {
        // internal node
        memset(dad_branch->partial_pars, 255, (pars_size - nptn - 1) * sizeof(int));
        memset(dad_branch->partial_pars + ptn_pars_start_id, 0, nptn * sizeof(int));
        UINT* partial_pars_dad = dad_branch->partial_pars;
        int partial_pars = 0;
        //UINT *partial_pars_child1 = NULL, *partial_pars_child2 = NULL;
        // take the intersection of two child states (with &= bit operation)
        FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME) {
            computePartialParsimony((PhyloNeighbor*)(*it), (PhyloNode*)node);
            /*
             if (!partial_pars_child1)
             partial_pars_child1 = ((PhyloNeighbor*) (*it))->partial_pars;
             else
             partial_pars_child2 = ((PhyloNeighbor*) (*it))->partial_pars;
             */
            UINT* partial_pars_child = ((PhyloNeighbor*)(*it))->partial_pars;
            for (int i = 0; i < pars_size - 1; i++)
                partial_pars_dad[i] &= partial_pars_child[i];
            partial_pars += partial_pars_child[pars_size - 1];
            for (int p = 0; p < nptn; p++)
                dad_branch->partial_pars[ptn_pars_start_id + p] += ((PhyloNeighbor*)(*it))->partial_pars[ptn_pars_start_id + p];
        }
        //assert(partial_pars_child1 && partial_pars_child2);
        // take the intersection of two bits block
        //for (int i = 0; i < pars_size - 1; i++)
        //    partial_pars_dad[i] = partial_pars_child1[i] & partial_pars_child2[i];
        //int partial_pars = partial_pars_child1[pars_size - 1] + partial_pars_child2[pars_size - 1];
        // now check if some intersection is empty, change to union (Fitch algorithm) and increase the parsimony score
        memset(bits_entry, 0, entry_size * sizeof(UINT));
        for (ptn = 0; ptn < aln->size(); ptn++)
            if (!aln->at(ptn).is_const) {
                getBitsBlock(partial_pars_dad, ptn, bits_entry);
                if (isEmptyBitsEntry(bits_entry)) {
                    FOR_NEIGHBOR_IT(node, dad, it2)if ((*it2)->node->name != ROOT_NAME) {
                        UINT* partial_pars_child = ((PhyloNeighbor*)(*it2))->partial_pars;
                        getBitsBlock(partial_pars_child, ptn, bits_entry_child);
                        unionBitsEntry(bits_entry, bits_entry_child, bits_entry);
                    }
                    //getBitsBlock(partial_pars_child2, ptn, bits_entry2);
                    //unionBitsEntry(bits_entry1, bits_entry2, bits_entry);
                    //cout << bits_entry[0] << " " << bits_entry[1] << endl;
                    setBitsBlock(partial_pars_dad, ptn, bits_entry);
                    partial_pars += aln->at(ptn).frequency;
                    dad_branch->partial_pars[ptn_pars_start_id + ptn] += 1;
                }
            }

        /*
         for (ptn = 0; ptn < aln->size(); ptn++)
         if (!aln->at(ptn).is_const) {
         getBitsBlock(partial_pars_dad, ptn, bits_entry);
         if (isEmptyBitsEntry(bits_entry)) {
         getBitsBlock(partial_pars_child1, ptn, bits_entry1);
         getBitsBlock(partial_pars_child2, ptn, bits_entry2);
         unionBitsEntry(bits_entry1, bits_entry2, bits_entry);
         //cout << bits_entry[0] << " " << bits_entry[1] << endl;
         setBitsBlock(partial_pars_dad, ptn, bits_entry);
         partial_pars += aln->at(ptn).frequency;
         }
         }*/
        partial_pars_dad[pars_size - 1] = partial_pars;
    }
    //delete[] bits_entry2;
    //delete[] bits_entry1;
    delete[] bits_entry_child;
    delete[] bits_entry;
}

int PhyloTree::computeParsimonyBranch(PhyloNeighbor* dad_branch, PhyloNode* dad, int* branch_subst) {
    PhyloNode* node = (PhyloNode*)dad_branch->node;
    PhyloNeighbor* node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
    assert(node_branch);
    if (!central_partial_pars)
        initializeAllPartialPars();
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode* tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor* tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }

    int nptn = aln->size();
    if (!_pattern_pars) _pattern_pars = aligned_alloc<BootValTypePars>(nptn + VCSIZE_USHORT);
    memset(_pattern_pars, 0, sizeof(BootValTypePars) * (nptn + VCSIZE_USHORT));

    if ((dad_branch->partial_lh_computed & 2) == 0)
        computePartialParsimony(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 2) == 0)
        computePartialParsimony(node_branch, node);
    // now combine likelihood at the branch

    int pars_size = getBitsBlockSize();
    int entry_size = getBitsEntrySize();
    int ptn_pars_start_id = pars_size - nptn - 1;
    //int nstates = aln->num_states;
    int i, ptn;
    int tree_pars = 0;

    if (aln->num_states == 4 && aln->seq_type == SEQ_DNA) {
        // ULTRAFAST VERSION FOR DNA
        for (ptn = 0; ptn < aln->size(); ptn += 8) {
            UINT states_left = node_branch->partial_pars[ptn / 8];
            UINT states_right = dad_branch->partial_pars[ptn / 8];
            UINT states_dad = 0;
            int maxi = aln->size() - ptn;
            if (maxi > 8) maxi = 8;
            for (i = 0; i < maxi; i++) {
                UINT state_left = (states_left >> (i * 4)) & 15;
                UINT state_right = (states_right >> (i * 4)) & 15;
                UINT state_both = state_left | (state_right << 4);
                states_dad |= dna_fitch_result[state_both] << (i * 4);
                tree_pars += dna_fitch_step[state_both] * aln->at(ptn + i).frequency;
                _pattern_pars[ptn + i] = node_branch->partial_pars[ptn_pars_start_id + ptn + i] +
                    dad_branch->partial_pars[ptn_pars_start_id + ptn + i] + dna_fitch_step[state_both];
                if(add_row)
                    save_branch_states_dad[ptn/8] = states_dad, save_branch_fitch_result[ptn+i] = dna_fitch_step[state_both];
            }
        }
        //		// the remaining bits
        //		UINT states_left = node_branch->partial_pars[ptn/8];
        //		UINT states_right = dad_branch->partial_pars[ptn/8];
        //		int maxi = aln->size() - ptn;
        //		for (i = 0; i< maxi; i++) {
        //			UINT state_left = (states_left >> (i*4)) & 15;
        //			UINT state_right = (states_right >> (i*4)) & 15;
        //			UINT state_both = state_left | (state_right << 4);
        //			_pattern_pars[ptn + i] += dna_fitch_step[state_both];
        //			tree_pars += dna_fitch_step[state_both] * aln->at(ptn+i).frequency;
        //		}
    }
    else if (aln->num_states == 20 && aln->seq_type == SEQ_PROTEIN) {
        // ULTRAFAST VERSION FOR PROTEIN
        UINT state_left[8], state_right[8];
        int id = 0;
        for (ptn = 0; ptn < aln->size(); ptn += 8, id += 5) {
            int maxi = aln->size() - ptn;
            if (maxi > 8) maxi = 8;
            decodeProtState(node_branch->partial_pars + id, state_left, maxi);
            decodeProtState(dad_branch->partial_pars + id, state_right, maxi);
            for (i = 0; i < maxi; i++) {
                if (!(state_left[i] & state_right[i])) {
                    tree_pars += aln->at(ptn + i).frequency;
                    _pattern_pars[ptn + i] = node_branch->partial_pars[ptn_pars_start_id + ptn + i] +
                        dad_branch->partial_pars[ptn_pars_start_id + ptn + i] + 1;
                }
                else
                    _pattern_pars[ptn + i] = node_branch->partial_pars[ptn_pars_start_id + ptn + i] +
                    dad_branch->partial_pars[ptn_pars_start_id + ptn + i];
            }
        }
    }
    else {
        // NORMAL VERSION FOR ALL #STATES
        UINT* partial_pars = newBitsBlock();
        UINT* bits_entry = new UINT[entry_size];
        for (i = 0; i < pars_size - 1; i++)
            partial_pars[i] = (node_branch->partial_pars[i] & dad_branch->partial_pars[i]);

        for (ptn = 0; ptn < aln->size(); ptn++)
            if (!aln->at(ptn).is_const) {
                getBitsBlock(partial_pars, ptn, bits_entry);
                if (isEmptyBitsEntry(bits_entry)) {
                    tree_pars += aln->at(ptn).frequency;
                    _pattern_pars[ptn] = node_branch->partial_pars[ptn_pars_start_id + ptn] +
                        dad_branch->partial_pars[ptn_pars_start_id + ptn] + 1;
                }
                else
                    _pattern_pars[ptn + i] = node_branch->partial_pars[ptn_pars_start_id + ptn] +
                    dad_branch->partial_pars[ptn_pars_start_id + ptn];
            }
        delete[] bits_entry;
        delete[] partial_pars;
    }
    if (branch_subst)
        *branch_subst = tree_pars;
    
    // cout << "fitch: " << tree_pars << '\n';

    tree_pars += node_branch->partial_pars[pars_size - 1] + dad_branch->partial_pars[pars_size - 1];

    return tree_pars;
}

int PhyloTree::computeParsimony() {
    assert(root->isLeaf());
    PhyloNeighbor* nei = ((PhyloNeighbor*)root->neighbors[0]);
    current_it = nei;
    assert(current_it);
    current_it_back = (PhyloNeighbor*)nei->node->findNeighbor(root);
    assert(current_it_back);

    int nptn = aln->size();
    if (_pattern_pars == NULL) _pattern_pars = aligned_alloc<BootValTypePars>(nptn + VCSIZE_USHORT);

    return computeParsimonyBranch((PhyloNeighbor*)root->neighbors[0], (PhyloNode*)root);
}

void PhyloTree::printParsimonyStates(PhyloNeighbor* dad_branch, PhyloNode* dad) {
    if (!dad) {
        dad = (PhyloNode*)root;
        dad_branch = (PhyloNeighbor*)root->neighbors[0];
        cout << "Parsimonious states for every node and site: " << endl;
    }
    int site;
    cout << "States for node ";
    int max_len = aln->getMaxSeqNameLength();
    if (max_len < 3)
        max_len = 3;
    cout.width(max_len);
    if (!dad_branch->node->name.empty())
        cout << left << dad_branch->node->name;
    else
        cout << left << dad_branch->node->id;
    cout << " are ";
    UINT* bits_entry = new UINT[getBitsEntrySize()];
    for (site = 0; site < aln->getNSite(); site++) {
        int ptn = aln->getPatternID(site);
        getBitsBlock(dad_branch->partial_pars, ptn, bits_entry);
        if (aln->at(ptn).is_const) {
            int state = aln->at(ptn)[0];
            if (state < aln->num_states)
                setBitsEntry(bits_entry, state);
            else {
                memset(bits_entry, 0, sizeof(UINT) * getBitsEntrySize());
                bits_entry[0] = state - (aln->num_states - 1);
                ;
            }
        }
        cout << "{";
        bool first = true;
        for (int i = 0; i < aln->num_states; i++)
            if (getBitsEntry(bits_entry, i)) {
                cout << ((!first) ? "," : "") << i;
                first = false;
            }
        cout << "}\t";
    }
    cout << endl;
    delete[] bits_entry;
    FOR_NEIGHBOR_IT(dad_branch->node, dad, it)printParsimonyStates((PhyloNeighbor*)(*it), (PhyloNode*)(dad_branch->node));
}

int PhyloTree::computeParsimonyScore(int ptn, int& states, PhyloNode* node, PhyloNode* dad) {
    int score = 0;
    states = 0;
    if (!node)
        node = (PhyloNode*)root;
    if (node->degree() > 3)
        outError("Does not work with multifurcating tree");
    if (verbose_mode == VB_DEBUG)
        cout << ptn << " " << node->id << "  " << node->name << endl;

    if (node->isLeaf()) {
        char state;
        if (node->name == ROOT_NAME) {
            state = aln->STATE_UNKNOWN;
        }
        else {
            assert(node->id < aln->getNSeq());
            state = (*aln)[ptn][node->id];
        }
        if (state == aln->STATE_UNKNOWN) {
            states = (1 << aln->num_states) - 1;
        }
        else if (state < aln->num_states)
            states = (1 << state);
        else {
            // ambiguous character, for DNA, RNA
            states = state - 3;
        }
    }
    if (!node->isLeaf() || node == root) {
        int union_states = 0;
        int intersect_states = (1 << aln->num_states) - 1;
        if (states != 0) {
            union_states = states;
            intersect_states = states;
        }

        FOR_NEIGHBOR_IT(node, dad, it) {
            int states_child;
            int score_child = computeParsimonyScore(ptn, states_child, (PhyloNode*)((*it)->node), node);
            union_states |= states_child;
            intersect_states &= states_child;
            score += score_child;
        }
        if (intersect_states)
            states = intersect_states;
        else {
            states = union_states;
            score++;
        }
    }
    return score;
}

int PhyloTree::computeParsimonyScore() {
    assert(root && root->isLeaf());

    int score = 0;
    for (int ptn = 0; ptn < aln->size(); ptn++)
        if (!aln->at(ptn).is_const) {
            int states;
            score += computeParsimonyScore(ptn, states) * (*aln)[ptn].frequency;
        }
    return score;
}

/****************************************************************************
 Nearest Neighbor Interchange with parsimony
 ****************************************************************************/

double PhyloTree::swapNNI(double cur_score, PhyloNode* node1, PhyloNode* node2) {
    assert(node1->degree() == 3 && node2->degree() == 3);
    FOR_NEIGHBOR_DECLARE(node1, node2, it1)
        break;
    Node* node1_nei = (*it1)->node;

    FOR_NEIGHBOR_IT(node2, node1, it2) {
        // do the NNI swap
        Node* node2_nei = (*it2)->node;
        node1->updateNeighbor(node1_nei, node2_nei);
        node1_nei->updateNeighbor(node1, node2);
        node2->updateNeighbor(node2_nei, node1_nei);
        node2_nei->updateNeighbor(node2, node1);

        // compute the score of the swapped topology
        double score = computeParsimonyScore();
        // if better: return
        if (score < cur_score) return score;
        // else, swap back
        node1->updateNeighbor(node2_nei, node1_nei);
        node1_nei->updateNeighbor(node2, node1);
        node2->updateNeighbor(node1_nei, node2_nei);
        node2_nei->updateNeighbor(node1, node2);
    }
    return cur_score;
}

double PhyloTree::searchNNI(double cur_score, PhyloNode* node, PhyloNode* dad) {
    if (!node)
        node = (PhyloNode*)root;
    if (!node->isLeaf() && dad && !dad->isLeaf()) {
        double score = swapNNI(cur_score, node, dad);
        if (score < cur_score)
            return score;
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        double score = searchNNI(cur_score, (PhyloNode*)(*it)->node, node);
        if (score < cur_score) return score;
    }
    return cur_score;
}

void PhyloTree::searchNNI() {
    cout << "Search with Nearest Neighbor Interchange..." << endl;
    double cur_score = computeParsimonyScore();
    do {
        double score = searchNNI(cur_score);
        if (score >= cur_score)
            break;
        cout << "Better score found: " << score << endl;
        cur_score = score;
    } while (true);
}

/****************************************************************************
 Stepwise addition (greedy) by maximum parsimony
 ****************************************************************************/

 // random generator function:
 //ptrdiff_t myrandom(ptrdiff_t i) {
 //    return random_int(i);
 //}

 // pointer object to it:
 //ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;

void PhyloTree::computeParsimonyTree(const char* out_prefix, Alignment* alignment) {
    //    cout << "Computing parsimony tree by random stepwise addition..." << endl;
    double start_time = getCPUTime();
    aln = alignment;
    int size = aln->getNSeq();
    if (size < 3)
        outError(ERR_FEW_TAXA);

    root = newNode();
    Node* new_taxon;

    IntVector taxon_order;
    taxon_order.resize(size);
    for (int i = 0; i < size; i++)
        taxon_order[i] = i;
    // randomize the addition order
//    random_shuffle(taxon_order.begin(), taxon_order.end(), p_myrandom);
    my_random_shuffle(taxon_order.begin(), taxon_order.end());

    // create initial tree with 3 taxa
    for (leafNum = 0; leafNum < 3; leafNum++) {
        if (verbose_mode >= VB_MAX)
            cout << "Add " << aln->getSeqName(taxon_order[leafNum]) << " to the tree" << endl;
        new_taxon = newNode(taxon_order[leafNum], aln->getSeqName(taxon_order[leafNum]).c_str());
        root->addNeighbor(new_taxon, -1.0);
        new_taxon->addNeighbor(root, -1.0);
    }
    root = findNodeID(taxon_order[0]);

    // stepwise adding the next taxon
    for (leafNum = 3; leafNum < size; leafNum++) {
        if (verbose_mode >= VB_MAX)
            cout << "Add " << aln->getSeqName(taxon_order[leafNum]) << " to the tree";
        initializeAllPartialPars();
        clearAllPartialLH();
        // allocate a new taxon and a new adjacent internal node
        new_taxon = newNode(taxon_order[leafNum], aln->getSeqName(taxon_order[leafNum]).c_str());
        Node* added_node = newNode();
        added_node->addNeighbor(new_taxon, -1.0);
        new_taxon->addNeighbor(added_node, -1.0);
        ((PhyloNeighbor*)added_node->findNeighbor(new_taxon))->partial_pars = newBitsBlock();
        ((PhyloNeighbor*)new_taxon->findNeighbor(added_node))->partial_pars = newBitsBlock();
        // preserve two neighbors
        added_node->addNeighbor((Node*)1, -1.0);
        added_node->addNeighbor((Node*)2, -1.0);

        Node* target_node = NULL;
        Node* target_dad = NULL;
        int score = addTaxonMPFast(added_node, target_node, target_dad, root->neighbors[0]->node, root);
        delete[]((PhyloNeighbor*)new_taxon->findNeighbor(added_node))->partial_pars;
        delete[]((PhyloNeighbor*)added_node->findNeighbor(new_taxon))->partial_pars;
        if (verbose_mode >= VB_MAX)
            cout << ", score = " << score << endl;
        // now insert the new node in the middle of the branch node-dad
        //double len = target_dad->findNeighbor(target_node)->length;
        target_node->updateNeighbor(target_dad, added_node, -1.0);
        target_dad->updateNeighbor(target_node, added_node, -1.0);
        added_node->updateNeighbor((Node*)1, target_node, -1.0);
        added_node->updateNeighbor((Node*)2, target_dad, -1.0);
        // compute the likelihood
        //clearAllPartialLh();
        //optimizeAllBranches();
        //optimizeNNI();
    }

    nodeNum = 2 * leafNum - 2;
    setAlignment(alignment);
    initializeAllPartialPars();
    clearAllPartialLH();
    fixNegativeBranch(true);
    //    cout << "Time taken: " << getCPUTime() - start_time << " sec" << endl;
    if (out_prefix) {
        string file_name = out_prefix;
        file_name += ".parstree";
        printTree(file_name.c_str(), WT_NEWLINE);
    }
}

int PhyloTree::addTaxonMPFast(Node* added_node, Node*& target_node, Node*& target_dad, Node* node, Node* dad) {
    Neighbor* dad_nei = dad->findNeighbor(node);
    //Node *added_taxon = added_node->neighbors[0]->node;
    Node* added_taxon = NULL;
    for (int i = 0; i < 3; i++) {
        if (added_node->neighbors[i]->node != (Node*)1 && added_node->neighbors[i]->node != (Node*)2)
            added_taxon = added_node->neighbors[i]->node;
    }

    //    Node *added_taxon;
    //    for (NeighborVec::iterator it = (added_node)->neighbors.begin(); it != (added_node)->neighbors.end(); it++) {
    //    	if ( (*it)->node->isLeaf() ) {
    //    		added_taxon = (*it)->node;
    //    	}
    //    }

    // now insert the new node in the middle of the branch node-dad
    double len = dad_nei->length;
    node->updateNeighbor(dad, added_node, len / 2.0);
    dad->updateNeighbor(node, added_node, len / 2.0);
    added_node->updateNeighbor((Node*)1, node, len / 2.0);
    added_node->updateNeighbor((Node*)2, dad, len / 2.0);
    ((PhyloNeighbor*)added_node->findNeighbor(node))->partial_pars =
        ((PhyloNeighbor*)dad->findNeighbor(added_node))->partial_pars;
    ((PhyloNeighbor*)added_node->findNeighbor(dad))->partial_pars =
        ((PhyloNeighbor*)node->findNeighbor(added_node))->partial_pars;
    ((PhyloNeighbor*)added_node->findNeighbor(node))->partial_lh_computed = ((PhyloNeighbor*)dad->findNeighbor(
        added_node))->partial_lh_computed;
    ((PhyloNeighbor*)added_node->findNeighbor(dad))->partial_lh_computed = ((PhyloNeighbor*)node->findNeighbor(
        added_node))->partial_lh_computed;
    // compute the likelihood
    //clearAllPartialLh();
    ((PhyloNeighbor*)added_taxon->findNeighbor(added_node))->clearPartialLh();
    int best_score = computeParsimonyBranch((PhyloNeighbor*)added_node->neighbors[0], (PhyloNode*)added_node);
    target_node = node;
    target_dad = dad;
    // remove the added node
    node->updateNeighbor(added_node, dad, len);
    dad->updateNeighbor(added_node, node, len);
    added_node->updateNeighbor(node, (Node*)1, len);
    added_node->updateNeighbor(dad, (Node*)2, len);

    // now tranverse the tree downwards

    FOR_NEIGHBOR_IT(node, dad, it) {
        Node* target_node2;
        Node* target_dad2;
        double score = addTaxonMPFast(added_node, target_node2, target_dad2, (*it)->node, node);
        if (score < best_score) {
            best_score = score;
            target_node = target_node2;
            target_dad = target_dad2;
        }
    }
    return best_score;

}

int PhyloTree::addTaxonMP(Node* added_node, Node*& target_node, Node*& target_dad, Node* node, Node* dad) {
    Neighbor* dad_nei = dad->findNeighbor(node);

    // now insert the new node in the middle of the branch node-dad
    double len = dad_nei->length;
    node->updateNeighbor(dad, added_node, len / 2.0);
    dad->updateNeighbor(node, added_node, len / 2.0);
    added_node->updateNeighbor((Node*)1, node, len / 2.0);
    added_node->updateNeighbor((Node*)2, dad, len / 2.0);
    // compute the likelihood
    //clearAllPartialLh();
    int best_score = computeParsimonyScore();
    target_node = node;
    target_dad = dad;
    // remove the added node
    node->updateNeighbor(added_node, dad, len);
    dad->updateNeighbor(added_node, node, len);
    added_node->updateNeighbor(node, (Node*)1, len);
    added_node->updateNeighbor(dad, (Node*)2, len);

    // now tranverse the tree downwards

    FOR_NEIGHBOR_IT(node, dad, it) {
        Node* target_node2;
        Node* target_dad2;
        double score = addTaxonMP(added_node, target_node2, target_dad2, (*it)->node, node);
        if (score < best_score) {
            best_score = score;
            target_node = target_node2;
            target_dad = target_dad2;
        }
    }
    return best_score;
}

void PhyloTree::growTreeMP(Alignment* alignment) {

    cout << "Stepwise addition using maximum parsimony..." << endl;
    aln = alignment;
    int size = aln->getNSeq();
    if (size < 3)
        outError(ERR_FEW_TAXA);

    root = newNode();
    Node* new_taxon;

    // create initial tree with 3 taxa
    for (leafNum = 0; leafNum < 3; leafNum++) {
        if (verbose_mode >= VB_MAX)
            cout << "Add " << aln->getSeqName(leafNum) << " to the tree" << endl;
        new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
        root->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(root, 1.0);
    }
    root = findNodeID(0);
    //optimizeAllBranches();

    // stepwise adding the next taxon
    for (leafNum = 3; leafNum < size; leafNum++) {
        if (verbose_mode >= VB_MAX)
            cout << "Add " << aln->getSeqName(leafNum) << " to the tree";
        // allocate a new taxon and a new ajedcent internal node
        new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
        Node* added_node = newNode();
        added_node->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(added_node, 1.0);

        // preserve two neighbors
        added_node->addNeighbor((Node*)1, 1.0);
        added_node->addNeighbor((Node*)2, 1.0);

        Node* target_node = NULL;
        Node* target_dad = NULL;
        int score = addTaxonMP(added_node, target_node, target_dad, root->neighbors[0]->node, root);
        if (verbose_mode >= VB_MAX)
            cout << ", score = " << score << endl;
        // now insert the new node in the middle of the branch node-dad
        double len = target_dad->findNeighbor(target_node)->length;
        target_node->updateNeighbor(target_dad, added_node, len / 2.0);
        target_dad->updateNeighbor(target_node, added_node, len / 2.0);
        added_node->updateNeighbor((Node*)1, target_node, len / 2.0);
        added_node->updateNeighbor((Node*)2, target_dad, len / 2.0);
        // compute the likelihood
        //clearAllPartialLh();
        //optimizeAllBranches();
        //optimizeNNI();
    }

    nodeNum = 2 * leafNum - 2;
}

/****************************************************************************
 likelihood function
 ****************************************************************************/

void PhyloTree::initializeAllPartialLh() {
    // TODO: this will cause bug for partition model
    if (params->maximum_parsimony) {
        initializeAllPartialPars();
        return;
    }
    int index, indexlh;
    int numStates = model->num_states;
    // Minh's question: why getAlnNSite() but not getAlnNPattern() ?
    //size_t mem_size = ((getAlnNSite() % 2) == 0) ? getAlnNSite() : (getAlnNSite() + 1);
    size_t nptn = getAlnNPattern() + numStates; // extra #numStates for ascertainment bias correction
#ifdef __AVX
    size_t mem_size = ((nptn + 3) / 4) * 4;
#else
    size_t mem_size = ((nptn % 2) == 0) ? nptn : (nptn + 1);
#endif
    size_t block_size = mem_size * numStates * site_rate->getNRate();
    if (!tmp_partial_lh1) {
        tmp_partial_lh1 = newPartialLh();
        // FOR TUNG: below is wrong because you lost the actual pointer to be deleted afterwards
        //if (((intptr_t) tmp_partial_lh1) % 16 != 0)
        //    tmp_partial_lh1 = tmp_partial_lh1 + 1;
    }
    if (!tmp_partial_lh2) {
        tmp_partial_lh2 = newPartialLh();
        //if (((intptr_t) tmp_partial_lh2) % 16 != 0)
        //    tmp_partial_lh2 = tmp_partial_lh2 + 1;
    }

    if (!tmp_anscentral_state_prob1)
        tmp_anscentral_state_prob1 = new double[numStates];
    if (!tmp_anscentral_state_prob2)
        tmp_anscentral_state_prob2 = new double[numStates];
    //if (!tmp_ptn_rates)
    //	tmp_ptn_rates = new double[alnSize]
    if (!tmp_scale_num1)
        tmp_scale_num1 = newScaleNum();
    if (!tmp_scale_num2)
        tmp_scale_num2 = newScaleNum();
    // make sure _pattern_lh size is divisible by 4 (e.g., 9->12, 14->16)
    if (!_pattern_lh)
        _pattern_lh = aligned_alloc_double(mem_size);
    if (!_pattern_lh_cat)
        _pattern_lh_cat = new double[mem_size * site_rate->getNDiscreteRate()];
    if (!theta_all)
        theta_all = aligned_alloc_double(block_size);
    if (!ptn_freq)
        ptn_freq = aligned_alloc_double(mem_size);
    if (!ptn_invar)
        ptn_invar = aligned_alloc_double(mem_size);
    initializeAllPartialLh(index, indexlh);
    assert(index == (nodeNum - 1) * 2);
    if (sse == LK_EIGEN || sse == LK_EIGEN_SSE)
        assert(indexlh == (nodeNum - 1) * 2 - leafNum);
    else
        assert(indexlh == (nodeNum - 1) * 2);
    clearAllPartialLH();

}

void PhyloTree::deleteAllPartialLh() {
    if (central_partial_lh) {
        delete[] central_partial_lh;
    }
    if (central_scale_num) {
        delete[] central_scale_num;
    }
    if (central_partial_pars)
        delete[] central_partial_pars;
    central_partial_lh = NULL;
    central_scale_num = NULL;
    central_partial_pars = NULL;
    clearAllPartialLH();
}

uint64_t PhyloTree::getMemoryRequired() {
    size_t nptn = aln->getNPattern() + aln->num_states; // +num_states for ascertainment bias correction
#ifdef __AVX
    // block size must be divisible by 4
    uint64_t block_size = ((nptn + 3) / 4) * 4;
#else
    // block size must be divisible by 2
    uint64_t block_size = ((nptn % 2) == 0) ? nptn : (nptn + 1);
#endif
    block_size = block_size * aln->num_states;
    if (site_rate)
        block_size *= site_rate->getNRate();
    uint64_t mem_size = ((uint64_t)leafNum * 4 - 6) * block_size + 2 + (leafNum - 1) * 4 * nptn * sizeof(UBYTE);
    if (sse == LK_EIGEN || sse == LK_EIGEN_SSE)
        mem_size -= ((uint64_t)leafNum) * ((uint64_t)block_size - nptn * sizeof(UBYTE));

    return mem_size;
}

void PhyloTree::initializeAllPartialLh(int& index, int& indexlh, PhyloNode* node, PhyloNode* dad) {
    size_t pars_block_size = getBitsBlockSize();
    size_t nptn = aln->size() + aln->num_states; // +num_states for ascertainment bias correction
#ifdef __AVX
    // block size must be divisible by 4
    size_t block_size = ((nptn + 3) / 4) * 4;
#else
    // block size must be divisible by 2
    size_t block_size = ((nptn % 2) == 0) ? nptn : (nptn + 1);
#endif
    size_t scale_block_size = nptn;

    block_size = block_size * model->num_states * site_rate->getNRate();
    if (!node) {
        node = (PhyloNode*)root;
        // allocate the big central partial likelihoods memory
        if (!central_partial_lh) {
            uint64_t tip_partial_lh_size = aln->num_states * (aln->STATE_UNKNOWN + 1);
            /*
            switch (aln->seq_type) {
            case SEQ_DNA: tip_partial_lh_size *=16; break; // including ambiguous nt and gap
            case SEQ_PROTEIN: tip_partial_lh_size *=23; break; // including 2 ambiguous aa and gap
            default: tip_partial_lh_size *=(aln->num_states+1); break; // including gap
            }*/
            uint64_t mem_size = ((uint64_t)leafNum * 4 - 6) * (uint64_t)block_size + 2 + tip_partial_lh_size;
            if (sse == LK_EIGEN || sse == LK_EIGEN_SSE)
                mem_size -= (uint64_t)leafNum * (uint64_t)block_size;
            try {
                central_partial_lh = new double[mem_size];
            }
            catch (std::bad_alloc& ba) {
                outError("Not enough memory for partial likelihood vectors (bad_alloc)");
            }
            //central_partial_lh = (double*) Eigen::internal::conditional_aligned_malloc<true>((leafNum-1)*4*block_size);
            if (!central_partial_lh)
                outError("Not enough memory for partial likelihood vectors");
            size_t mem_shift = 0;
            if (((intptr_t)central_partial_lh) % MEM_ALIGNMENT != 0)
                mem_shift = (MEM_ALIGNMENT - (((intptr_t)central_partial_lh) % MEM_ALIGNMENT)) / sizeof(double);
            if (sse == LK_EIGEN || sse == LK_EIGEN_SSE)
                tip_partial_lh = central_partial_lh + (((nodeNum - 1) * 2 - leafNum) * block_size + mem_shift);
            else
                tip_partial_lh = central_partial_lh + (((nodeNum - 1) * 2) * block_size + mem_shift);
        }

        if (!central_scale_num) {
            uint64_t mem_size = (leafNum - 1) * 4 * scale_block_size;
            if (sse == LK_EIGEN || sse == LK_EIGEN_SSE)
                mem_size -= (uint64_t)leafNum * (uint64_t)scale_block_size;
            if (verbose_mode >= VB_MED)
                cout << "Allocating " << mem_size * sizeof(UBYTE) << " bytes for scale num vectors" << endl;
            try {
                central_scale_num = new UBYTE[mem_size];
            }
            catch (std::bad_alloc& ba) {
                outError("Not enough memory for scale num vectors (bad_alloc)");
            }
            if (!central_scale_num)
                outError("Not enough memory for scale num vectors");
        }

        if (!central_partial_pars) {
            if (verbose_mode >= VB_MED)
                cout << "Allocating " << (leafNum - 1) * 4 * pars_block_size * sizeof(UINT)
                << " bytes for partial parsimony vectors" << endl;
            try {
                central_partial_pars = new UINT[(leafNum - 1) * 4 * pars_block_size];
            }
            catch (std::bad_alloc& ba) {
                outError("Not enough memory for partial parsimony vectors (bad_alloc)");
            }
            if (!central_partial_pars)
                outError("Not enough memory for partial parsimony vectors");
        }
        index = 0;
        indexlh = 0;
    }
    if (dad) {
        // make memory alignment
        size_t mem_shift = 0;
        if (((intptr_t)central_partial_lh) % MEM_ALIGNMENT != 0)
            mem_shift = (MEM_ALIGNMENT - (((intptr_t)central_partial_lh) % MEM_ALIGNMENT)) / sizeof(double);
        // assign a region in central_partial_lh to both Neihgbors (dad->node, and node->dad)
        PhyloNeighbor* nei = (PhyloNeighbor*)node->findNeighbor(dad);
        //assert(!nei->partial_lh);
        if (nei->node->isLeaf() && (sse == LK_EIGEN || sse == LK_EIGEN_SSE)) {
            nei->partial_lh = NULL; // do not allocate memory for tip, use tip_partial_lh instead
            nei->scale_num = NULL;
        }
        else {
            nei->scale_num = central_scale_num + (indexlh * scale_block_size);
            nei->partial_lh = central_partial_lh + (indexlh * block_size + mem_shift);
            indexlh++;
        }
        nei->partial_pars = central_partial_pars + (index * pars_block_size);
        index++;
        nei = (PhyloNeighbor*)dad->findNeighbor(node);
        //assert(!nei->partial_lh);
        if (nei->node->isLeaf() && (sse == LK_EIGEN || sse == LK_EIGEN_SSE)) {
            nei->partial_lh = NULL; // do not allocate memory for tip, use tip_partial_lh instead
            nei->scale_num = NULL;
        }
        else {
            nei->scale_num = central_scale_num + ((indexlh)*scale_block_size);
            nei->partial_lh = central_partial_lh + (indexlh * block_size + mem_shift);
            indexlh++;
        }
        nei->partial_pars = central_partial_pars + (index * pars_block_size);
        index++;
        assert(index < nodeNum * 2 - 1);
    }
    FOR_NEIGHBOR_IT(node, dad, it) initializeAllPartialLh(index, indexlh, (PhyloNode*)(*it)->node, node);
}

double* PhyloTree::newPartialLh() {
    double* ret = new double[(aln->size() + aln->num_states + 3) * aln->num_states * site_rate->getNRate()];
    return ret;
}

UINT* PhyloTree::newPartialPars() {
    size_t pars_block_size = getBitsBlockSize();
    UINT* ret = new UINT[pars_block_size];
    return ret;
}

int PhyloTree::getPartialLhBytes() {
    return ((aln->size() + aln->num_states + 3) * aln->num_states * site_rate->getNRate()) * sizeof(double);
}

int PhyloTree::getScaleNumBytes() {
    return (aln->size() + aln->num_states) * sizeof(UBYTE);
}

UBYTE* PhyloTree::newScaleNum() {
    return new UBYTE[aln->size() + aln->num_states];
}

double PhyloTree::computeLikelihood(double* pattern_lh) {
    if (params->maximum_parsimony) {
        clearAllPartialLH();
        return -computeParsimony();
    }
    assert(model);
    assert(site_rate);
    assert(root->isLeaf());
    PhyloNeighbor* nei = ((PhyloNeighbor*)root->neighbors[0]);
    current_it = nei;
    assert(current_it);
    current_it_back = (PhyloNeighbor*)nei->node->findNeighbor(root);
    assert(current_it_back);

    double score;
    string root_name = ROOT_NAME;
    Node* vroot = findLeafName(root_name);
    if (root_state != aln->STATE_UNKNOWN && vroot) {
        if (verbose_mode >= VB_DEBUG)
            cout << __func__ << " HIT ROOT STATE " << endl;
        score = computeLikelihoodRooted((PhyloNeighbor*)vroot->neighbors[0], (PhyloNode*)vroot);
    }
    else {
        score = computeLikelihoodBranch(nei, (PhyloNode*)root, pattern_lh);
    }
    if (pattern_lh && nei->lh_scale_factor < 0.0) {
        int nptn = aln->getNPattern();
        //double check_score = 0.0;
        for (int i = 0; i < nptn; i++) {
            pattern_lh[i] += max(nei->scale_num[i], UBYTE(0)) * LOG_SCALING_THRESHOLD;
            //check_score += (pattern_lh[i] * (aln->at(i).frequency));
        }
        /*       if (fabs(score - check_score) > 1e-6) {
         cout << "score = " << score << " check_score = " << check_score << endl;
         outError("Scaling error ", __func__);
         }*/
    }
    return score;
}

double PhyloTree::computeLikelihoodRooted(PhyloNeighbor* dad_branch, PhyloNode* dad) {
    double score = computeLikelihoodBranchNaive(dad_branch, dad);
    if (verbose_mode >= VB_DEBUG) {
        printTransMatrices(dad_branch->node, dad);
        /*
         FOR_NEIGHBOR_IT(dad_branch->node, dad, it) {
         PhyloNeighbor *pit = (PhyloNeighbor*)(*it);
         cout << pit->node->name << "\t" << pit->partial_lh[0] << endl;

         }*/
    }
    double* state_freq = new double[aln->num_states];
    model->getStateFrequency(state_freq);
    score -= log(state_freq[(int)root_state]);
    delete[] state_freq;
    return score;
}

void PhyloTree::computePatternLikelihood(double* ptn_lh, double* cur_logl, double* ptn_lh_cat) {
    if (params->maximum_parsimony) {
        computePatternParsimony(ptn_lh, cur_logl);
        return;
    }
    /*	if (!dad_branch) {
     dad_branch = (PhyloNeighbor*) root->neighbors[0];
     dad = (PhyloNode*) root;
     }*/
    int nptn = aln->getNPattern();
    int i;
    int ncat = site_rate->getNDiscreteRate();
    if (ptn_lh_cat) {
        // Right now only Naive version store _pattern_lh_cat!
        if (sse == LK_NORMAL || sse == LK_SSE)
            computeLikelihoodBranchNaive(current_it, (PhyloNode*)current_it_back->node);
        else {
            switch (aln->num_states) {
            case 4: computeLikelihoodBranchEigen<4>(current_it, (PhyloNode*)current_it_back->node); break;
            case 20: computeLikelihoodBranchEigen<20>(current_it, (PhyloNode*)current_it_back->node); break;
            case 2: computeLikelihoodBranchEigen<2>(current_it, (PhyloNode*)current_it_back->node); break;
            default: outError("Option unsupported yet for this sequence type. Contact author if you really need it."); break;
            }
        }

    }
    double sum_scaling = current_it->lh_scale_factor + current_it_back->lh_scale_factor;
    //double sum_scaling = 0.0;
    if (sum_scaling < 0.0) {
        if (current_it->lh_scale_factor == 0.0) {
            for (i = 0; i < nptn; i++) {
                ptn_lh[i] = _pattern_lh[i] + (max(UBYTE(0), current_it_back->scale_num[i])) * LOG_SCALING_THRESHOLD;
            }
        }
        else if (current_it_back->lh_scale_factor == 0.0) {
            for (i = 0; i < nptn; i++) {
                ptn_lh[i] = _pattern_lh[i] + (max(UBYTE(0), current_it->scale_num[i])) * LOG_SCALING_THRESHOLD;
            }
        }
        else {
            for (i = 0; i < nptn; i++) {
                ptn_lh[i] = _pattern_lh[i] + (max(UBYTE(0), current_it->scale_num[i]) +
                    max(UBYTE(0), current_it_back->scale_num[i])) * LOG_SCALING_THRESHOLD;
            }
        }
    }
    else {
        memmove(ptn_lh, _pattern_lh, nptn * sizeof(double));
    }
    if (ptn_lh_cat) {
        int offset = 0;
        if (sum_scaling == 0.0) {
            int nptncat = nptn * ncat;
            for (i = 0; i < nptncat; i++) {
                ptn_lh_cat[i] = log(_pattern_lh_cat[i]);
            }
        }
        else if (current_it->lh_scale_factor == 0.0) {
            for (i = 0; i < nptn; i++) {
                double scale = (max(UBYTE(0), current_it_back->scale_num[i])) * LOG_SCALING_THRESHOLD;
                for (int j = 0; j < ncat; j++, offset++)
                    ptn_lh_cat[offset] = log(_pattern_lh_cat[offset]) + scale;
            }
        }
        else if (current_it_back->lh_scale_factor == 0.0) {
            for (i = 0; i < nptn; i++) {
                double scale = (max(UBYTE(0), current_it->scale_num[i])) * LOG_SCALING_THRESHOLD;
                for (int j = 0; j < ncat; j++, offset++)
                    ptn_lh_cat[offset] = log(_pattern_lh_cat[offset]) + scale;
            }
        }
        else {
            for (i = 0; i < nptn; i++) {
                double scale = (max(UBYTE(0), current_it->scale_num[i]) +
                    max(UBYTE(0), current_it_back->scale_num[i])) * LOG_SCALING_THRESHOLD;
                for (int j = 0; j < ncat; j++, offset++)
                    ptn_lh_cat[offset] = log(_pattern_lh_cat[offset]) + scale;
            }
        }
    }
    //    if (cur_logl) {
    //        double check_score = 0.0;
    //        for (int i = 0; i < nptn; i++) {
    //            check_score += (ptn_lh[i] * (aln->at(i).frequency));
    //        }
    //        if (fabs(check_score - *cur_logl) > 0.01) {
    //            cout << *cur_logl << " " << check_score << endl;
    //            outError("Wrong PhyloTree::", __func__);
    //        }
    //    }
        //double score = computeLikelihoodBranch(dad_branch, dad, pattern_lh);
        //return score;
}

void PhyloTree::computePatternParsimony(double* ptn_npars, double* cur_npars) {
    if (!ptn_npars)
        outError("ERROR: No space allocated for the pattern parsimony vector.\n");

    if (!params)
        outError("No params detected!");
    //	if(params->snni && params->spr_parsimony){
    //		// Need to call pllComputePatternParsimony of sprparsimony.cpp
    //		pllComputePatternParsimony(dynamic_cast<IQTree*>(this)->pllInst, dynamic_cast<IQTree*>(this)->pllPartitions , ptn_npars, cur_npars);
    //	}
        // As IQTree::computeParsimonBranch already computed _pattern_pars
        // just copy from that
    int nptn = aln->getNPattern();
    int i;
    for (i = 0; i < nptn; i++) {
        // TODO: this is a bit inefficient, should change everthing to int operations
        ptn_npars[i] = -double(_pattern_pars[i]);
    }

    //	// If cur_npars is not NULL then check the correctness of pattern parsimony
    //	// IMHO, this is kind of a waste of time if not for debugging
    //
    //	if(cur_npars){
    //		int my_tree_npars = 0;
    //		for(i = 0; i < nptn; i++) my_tree_npars += _pattern_pars[i] * aln->at(i).frequency;
    //		if( *cur_npars != -double(my_tree_npars))
    //			outError("Wrong PhyloTree::computePatternParsimony()!");
    //		ofstream pars_log;
    //		pars_log.open("pars.log", ios::out|ios::app);
    //		pars_log << endl << *cur_npars << "\t" << my_tree_npars;
    //		if( *cur_npars != my_tree_npars) pars_log << " >> ERROR";
    //		pars_log.close();
    //	}
}

double PhyloTree::computeLogLVariance(double* ptn_lh, double tree_lh) {
    int i;
    int nptn = getAlnNPattern();
    int nsite = getAlnNSite();
    double* pattern_lh = ptn_lh;
    if (!ptn_lh) {
        pattern_lh = new double[nptn];
        computePatternLikelihood(pattern_lh);
    }
    IntVector pattern_freq;
    aln->getPatternFreq(pattern_freq);
    if (tree_lh == 0.0) {
        for (i = 0; i < nptn; i++)
            tree_lh += pattern_lh[i] * pattern_freq[i];
    }
    double avg_site_lh = tree_lh / nsite;
    double variance = 0.0;
    for (i = 0; i < nptn; i++) {
        double diff = (pattern_lh[i] - avg_site_lh);
        variance += diff * diff * pattern_freq[i];
    }
    if (!ptn_lh)
        delete[] pattern_lh;
    if (nsite <= 1)
        return 0.0;
    return variance * ((double)nsite / (nsite - 1.0));
}

double PhyloTree::computeLogLDiffVariance(double* pattern_lh_other, double* ptn_lh) {
    int i;
    int nptn = getAlnNPattern();
    int nsite = getAlnNSite();
    double* pattern_lh = ptn_lh;
    if (!ptn_lh) {
        pattern_lh = new double[nptn];
        computePatternLikelihood(pattern_lh);
    }
    IntVector pattern_freq;
    aln->getPatternFreq(pattern_freq);

    double avg_site_lh_diff = 0.0;
    for (i = 0; i < nptn; i++)
        avg_site_lh_diff += (pattern_lh[i] - pattern_lh_other[i]) * pattern_freq[i];
    avg_site_lh_diff /= nsite;
    double variance = 0.0;
    for (i = 0; i < nptn; i++) {
        double diff = (pattern_lh[i] - pattern_lh_other[i] - avg_site_lh_diff);
        variance += diff * diff * pattern_freq[i];
    }
    if (!ptn_lh)
        delete[] pattern_lh;
    if (nsite <= 1)
        return 0.0;
    return variance * ((double)nsite / (nsite - 1.0));
}

double PhyloTree::computeLogLDiffVariance(PhyloTree* other_tree, double* pattern_lh) {
    double* pattern_lh_other = new double[getAlnNPattern()];
    other_tree->computePatternLikelihood(pattern_lh_other);
    delete[] pattern_lh_other;
    double res = computeLogLDiffVariance(pattern_lh_other, pattern_lh);
    return res;
}

void PhyloTree::getUnmarkedNodes(PhyloNodeVector& unmarkedNodes, PhyloNode* node, PhyloNode* dad) {
    if (!node) {
        node = (PhyloNode*)root;
    }

    if (markedNodeList.find(node->id) == markedNodeList.end()) {
        int numUnmarkedNei = 0;
        for (NeighborVec::iterator it = (node)->neighbors.begin(); it != (node)->neighbors.end(); it++) {
            if (markedNodeList.find((*it)->node->id) == markedNodeList.end())
                numUnmarkedNei++;
        }
        if (numUnmarkedNei == 1)
            unmarkedNodes.push_back(node);
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        getUnmarkedNodes(unmarkedNodes, (PhyloNode*)(*it)->node, node);
    }
}

double PhyloTree::optimizeOneBranchLS(PhyloNode* node1, PhyloNode* node2) {
    if (!subTreeDistComputed) {
        if (params->ls_var_type == WLS_PAUPLIN) {
            computeNodeBranchDists();
            for (int i = 0; i < leafNum; i++)
                for (int j = 0; j < leafNum; j++)
                    var_matrix[i * leafNum + j] = pow(2.0, nodeBranchDists[i * nodeNum + j]);
        }
        computeSubtreeDists();
    }
    double A, B, C, D;
    A = B = C = D = 0;
    PhyloNode* nodeA = NULL, * nodeB = NULL, * nodeC = NULL, * nodeD = NULL;
    double lsBranch;

    // One of the node is a leaf
    if (node1->isLeaf() || node2->isLeaf()) {
        if (node1->isLeaf()) {
            // nodeA and nodeB are children of node2
            FOR_NEIGHBOR_IT(node2, node1, it) {
                if (A == 0) {
                    A = getNumTaxa((*it)->node, node2);
                    nodeA = (PhyloNode*)(*it)->node;
                }
                else {
                    B = getNumTaxa((*it)->node, node2);
                    nodeB = (PhyloNode*)(*it)->node;
                }
            }
            // nodeC is now node1
            nodeC = node1;
        }
        else {
            // nodeA and nodeB are children of node1
            FOR_NEIGHBOR_IT(node1, node2, it) {
                if (A == 0) {
                    A = getNumTaxa((*it)->node, node1);
                    nodeA = (PhyloNode*)(*it)->node;
                }
                else {
                    B = getNumTaxa((*it)->node, node1);
                    nodeB = (PhyloNode*)(*it)->node;
                }
            }
            // nodeC is now node1
            nodeC = node2;
        }
        assert(A != 0);
        assert(B != 0);
        string keyAC = nodePair2String(nodeA, nodeC);
        assert(subTreeDists.count(keyAC));
        double distAC = subTreeDists[keyAC];
        double weightAC = subTreeWeights[keyAC];
        string keyBC = nodePair2String(nodeB, nodeC);
        assert(subTreeDists.count(keyBC));
        double distBC = subTreeDists[keyBC];
        double weightBC = subTreeWeights[keyBC];
        string keyAB = nodePair2String(nodeA, nodeB);
        assert(subTreeDists.count(keyAB));
        double distAB = subTreeDists[keyAB];
        double weightAB = subTreeWeights[keyAB];
        if (params->ls_var_type == OLS/* || params->ls_var_type == FIRST_TAYLOR || params->ls_var_type == FITCH_MARGOLIASH
                || params->ls_var_type == SECOND_TAYLOR*/) {
            lsBranch = 0.5 * (distAC / A + distBC / B - distAB / (A * B));
        } /*else if (params->ls_var_type == PAUPLIN) {
            // TODO: Chua test bao gio
            outError("Paulin formula not supported yet");
            lsBranch = 0.5 * (distAC + distBC) - 0.5 * distAB;
        }*/ else {
        // weighted least square
            lsBranch = 0.5 * (distAC / weightAC + distBC / weightBC - distAB / weightAB);
        }
    }
    else { // Both node are internal node
        FOR_NEIGHBOR_IT(node1, node2, it) {
            if (A == 0) {
                A = getNumTaxa((*it)->node, node1);
                nodeA = (PhyloNode*)(*it)->node;
            }
            else {
                B = getNumTaxa((*it)->node, node1);
                nodeB = (PhyloNode*)(*it)->node;
            }
        }

        FOR_NEIGHBOR_IT(node2, node1, it) {
            if (C == 0) {
                C = getNumTaxa((*it)->node, node2);
                nodeC = (PhyloNode*)(*it)->node;
            }
            else {
                D = getNumTaxa((*it)->node, node2);
                nodeD = (PhyloNode*)(*it)->node;
            }
        }

        string keyAC = nodePair2String(nodeA, nodeC);
        assert(subTreeDists.count(keyAC));
        double distAC = subTreeDists[keyAC];
        double weightAC = subTreeWeights[keyAC];

        string keyBD = nodePair2String(nodeB, nodeD);
        assert(subTreeDists.count(keyBD));
        double distBD = subTreeDists[keyBD];
        double weightBD = subTreeWeights[keyBD];

        string keyBC = nodePair2String(nodeB, nodeC);
        assert(subTreeDists.count(keyBC));
        double distBC = subTreeDists[keyBC];
        double weightBC = subTreeWeights[keyBC];

        string keyAD = nodePair2String(nodeA, nodeD);
        assert(subTreeDists.count(keyAD));
        double distAD = subTreeDists[keyAD];
        double weightAD = subTreeWeights[keyAD];

        string keyAB = nodePair2String(nodeA, nodeB);
        assert(subTreeDists.count(keyAB));
        double distAB = subTreeDists[keyAB];
        double weightAB = subTreeWeights[keyAB];

        string keyCD = nodePair2String(nodeC, nodeD);
        assert(subTreeDists.count(keyCD));
        double distCD = subTreeDists[keyCD];
        double weightCD = subTreeWeights[keyCD];

        /*if (params->ls_var_type == PAUPLIN) {
            // this distance has a typo as also seen in Mihaescu & Pachter 2008
            //lsBranch = 0.25 * (distAC + distBD + distAD + distBC) - 0.5 * (distAB - distCD);
            outError("Paulin formula not supported yet");
            lsBranch = 0.25 * (distAC + distBD + distAD + distBC) - 0.5 * (distAB + distCD);
        } else*/ if (params->ls_var_type == OLS) {
            double gamma = (B * C + A * D) / ((A + B) * (C + D));
            lsBranch = 0.5 * (gamma * (distAC / (A * C) + distBD / (B * D))
                + (1 - gamma) * (distBC / (B * C) + distAD / (A * D))
                - distAB / (A * B) - distCD / (C * D));
        }
        else {
            // weighted least square
            double K = 1.0 / weightAC + 1.0 / weightBD + 1.0 / weightAD + 1.0 / weightBC;
            lsBranch =
                ((distAC / weightAC + distBD / weightBD) * (weightAD + weightBC) / (weightAD * weightBC) +
                    (distAD / weightAD + distBC / weightBC) * (weightAC + weightBD) / (weightAC * weightBD)) / K
                - distAB / weightAB - distCD / weightCD;
            lsBranch = 0.5 * lsBranch;
        }
    }
    return lsBranch;
}

void PhyloTree::updateSubtreeDists(NNIMove& nnimove) {
    assert(subTreeDistComputed);
    PhyloNode* nodeA = NULL, * nodeB = NULL, * nodeC = NULL, * nodeD = NULL;
    PhyloNode* node1 = nnimove.node1;
    PhyloNode* node2 = nnimove.node2;
    NeighborVec::iterator node1Nei_it = nnimove.node1Nei_it;
    NeighborVec::iterator node2Nei_it = nnimove.node2Nei_it;
    Neighbor* node1Nei = *(node1Nei_it);
    Neighbor* node2Nei = *(node2Nei_it);

    // ((A,C),(B,D))
    // C and D are the 2 subtree that get swapped
    FOR_NEIGHBOR_IT(node1, node2, it) {
        if ((*it)->id != node1Nei->id) {
            nodeA = (PhyloNode*)(*it)->node;
        }
        else {
            nodeC = (PhyloNode*)(*it)->node;
        }
    }

    assert(nodeA);
    assert(nodeC);

    FOR_NEIGHBOR_IT(node2, node1, it) {
        if ((*it)->id != node2Nei->id) {
            nodeB = (PhyloNode*)(*it)->node;
        }
        else {
            nodeD = (PhyloNode*)(*it)->node;
        }
    }

    assert(nodeB);
    assert(nodeD);

    NodeVector nodeListA, nodeListB, nodeListC, nodeListD;
    getAllNodesInSubtree(nodeA, node1, nodeListA);
    getAllNodesInSubtree(nodeC, node1, nodeListC);
    getAllNodesInSubtree(nodeB, node2, nodeListB);
    getAllNodesInSubtree(nodeD, node2, nodeListD);

    for (NodeVector::iterator it = nodeListA.begin(); it != nodeListA.end(); ++it) {
        string key = nodePair2String((*it), node2);
        double distB = subTreeDists.find(nodePair2String((*it), nodeB))->second;
        double distD = subTreeDists.find(nodePair2String((*it), nodeD))->second;
        double newDist = distB + distD;
        StringDoubleMap::iterator dist_it = subTreeDists.find(key);
        assert(dist_it != subTreeDists.end());
        dist_it->second = newDist;
    }

    for (NodeVector::iterator it = nodeListB.begin(); it != nodeListB.end(); ++it) {
        string key = nodePair2String((*it), node1);
        double distC = subTreeDists.find(nodePair2String((*it), nodeC))->second;
        double distA = subTreeDists.find(nodePair2String((*it), nodeA))->second;
        double newDist = distC + distA;
        StringDoubleMap::iterator dist_it = subTreeDists.find(key);
        assert(dist_it != subTreeDists.end());
        dist_it->second = newDist;
    }

    for (NodeVector::iterator it = nodeListC.begin(); it != nodeListC.end(); ++it) {
        string key = nodePair2String((*it), node2);
        double distD = subTreeDists.find(nodePair2String((*it), nodeD))->second;
        double distB = subTreeDists.find(nodePair2String((*it), nodeB))->second;
        double newDist = distD + distB;
        StringDoubleMap::iterator dist_it = subTreeDists.find(key);
        assert(dist_it != subTreeDists.end());
        dist_it->second = newDist;
    }

    for (NodeVector::iterator it = nodeListD.begin(); it != nodeListD.end(); ++it) {
        string key = nodePair2String((*it), node1);
        double distA = subTreeDists.find(nodePair2String((*it), nodeA))->second;
        double distC = subTreeDists.find(nodePair2String((*it), nodeC))->second;
        double newDist = distA + distC;
        StringDoubleMap::iterator dist_it = subTreeDists.find(key);
        assert(dist_it != subTreeDists.end());
        dist_it->second = newDist;
    }

    double distAB = subTreeDists.find(nodePair2String(nodeA, nodeB))->second;
    double distAD = subTreeDists.find(nodePair2String(nodeA, nodeD))->second;
    double distCB = subTreeDists.find(nodePair2String(nodeC, nodeB))->second;
    double distCD = subTreeDists.find(nodePair2String(nodeC, nodeD))->second;

    subTreeDists.find(nodePair2String(node1, node2))->second = distAB + distAD + distCB + distCD;

}

void PhyloTree::computeSubtreeDists() {
    PhyloNodeVector unmarkedNodes;
    subTreeDists.clear();
    subTreeWeights.clear();
    do {
        // Generate a list of unmarked node that is adjacent to exactly one unmarked nodes
        // Here we will work up the tree in a bottom up manner
        unmarkedNodes.clear();
        getUnmarkedNodes(unmarkedNodes);
        if (unmarkedNodes.size() == 0)
            break;

        for (PhyloNodeVector::iterator it = unmarkedNodes.begin(); it != unmarkedNodes.end(); ++it) {
            // if the node is an internal node then all of its child nodes should be marked
            // source_nei1 and source_nei2 are the 2 marked child node
            // nextNode is the other node, used for traversal
            PhyloNode* source_nei1 = NULL;
            PhyloNode* source_nei2 = NULL;
            PhyloNode* nextNode;
            if (!(*it)->isLeaf()) {
                // select the 2 marked child nodes
                for (NeighborVec::iterator it2 = (*it)->neighbors.begin(); it2 != (*it)->neighbors.end(); ++it2) {
                    if (markedNodeList.find((*it2)->node->id) != markedNodeList.end()) {
                        if (!source_nei1) {
                            source_nei1 = (PhyloNode*)(*it2)->node;
                        }
                        else {
                            source_nei2 = (PhyloNode*)(*it2)->node;
                        }
                    }
                    else {
                        nextNode = (PhyloNode*)(*it2)->node;
                    }
                }
                assert(source_nei1);
                assert(source_nei2);
            }
            else {
                nextNode = (PhyloNode*)(*it)->neighbors[0]->node;
            }
            // warning: 'nextNode' may be used uninitialized in this function
            computeAllSubtreeDistForOneNode((*it), source_nei1, source_nei2, (*it), nextNode);
            markedNodeList.insert(IntPhyloNodeMap::value_type((*it)->id, (*it)));
        }
    } while (true);
    markedNodeList.clear();
    subTreeDistComputed = true;
}

void PhyloTree::computeAllSubtreeDistForOneNode(PhyloNode* source, PhyloNode* source_nei1, PhyloNode* source_nei2,
    PhyloNode* node, PhyloNode* dad) {
    string key = nodePair2String(source, dad);
    double dist, weight;
    if (markedNodeList.find(dad->id) != markedNodeList.end()) {
        return;
    }
    else if (source->isLeaf() && dad->isLeaf()) {
        assert(dist_matrix);
        int nseq = aln->getNSeq();
        if (params->ls_var_type == OLS) {
            dist = dist_matrix[dad->id * nseq + source->id];
            weight = 1.0;
        }
        else {
            // this will take into account variances, also work for OLS since var = 1
            weight = 1.0 / var_matrix[dad->id * nseq + source->id];
            dist = dist_matrix[dad->id * nseq + source->id] * weight;
        }
        subTreeDists.insert(StringDoubleMap::value_type(key, dist));
        subTreeWeights.insert(StringDoubleMap::value_type(key, weight));
    }
    else if (!source->isLeaf() && dad->isLeaf()) {
        assert(source_nei1);
        assert(source_nei2);
        string key1 = nodePair2String(source_nei1, dad);
        assert(subTreeDists.find(key1) == subTreeDists.end());
        double dist1 = subTreeDists.find(key1)->second;
        double weight1 = subTreeWeights.find(key1)->second;
        string key2 = nodePair2String(source_nei2, dad);
        assert(subTreeDists.find(key2) == subTreeDists.end());
        double dist2 = subTreeDists.find(key2)->second;
        double weight2 = subTreeWeights.find(key2)->second;
        dist = dist1 + dist2;
        weight = weight1 + weight2;
        subTreeDists.insert(StringDoubleMap::value_type(key, dist));
        subTreeWeights.insert(StringDoubleMap::value_type(key, weight));
    }
    else {
        PhyloNode* dad_nei1 = NULL;
        PhyloNode* dad_nei2 = NULL;
        for (NeighborVec::iterator it = dad->neighbors.begin(); it != dad->neighbors.end(); ++it) {
            if ((*it)->node != node) {
                if (!dad_nei1) {
                    dad_nei1 = (PhyloNode*)(*it)->node;
                }
                else {
                    dad_nei2 = (PhyloNode*)(*it)->node;
                }
            }
        }
        assert(dad_nei1);
        assert(dad_nei2);
        computeAllSubtreeDistForOneNode(source, source_nei1, source_nei2, dad, dad_nei1);
        computeAllSubtreeDistForOneNode(source, source_nei1, source_nei2, dad, dad_nei2);
        string key1 = nodePair2String(source, dad_nei1);
        string key2 = nodePair2String(source, dad_nei2);
        assert(subTreeDists.find(key1) != subTreeDists.end());
        assert(subTreeDists.find(key2) != subTreeDists.end());
        double dist1 = subTreeDists.find(key1)->second;
        double weight1 = subTreeWeights.find(key1)->second;
        double dist2 = subTreeDists.find(key2)->second;
        double weight2 = subTreeWeights.find(key2)->second;
        dist = dist1 + dist2;
        weight = weight1 + weight2;
        subTreeDists.insert(StringDoubleMap::value_type(key, dist));
        subTreeWeights.insert(StringDoubleMap::value_type(key, weight));
    }
}

set<int> PhyloTree::computeNodeBranchDists(Node* node, Node* dad) {
    set<int>::iterator i, j;
    if (!nodeBranchDists) {
        cout << "nodeNum = " << nodeNum << endl;
        nodeBranchDists = new int[nodeNum * nodeNum];
    }
    if (!node) {
        memset(nodeBranchDists, 0, sizeof(int) * nodeNum * nodeNum);
        assert(root->isLeaf());
        dad = root;
        node = dad->neighbors[0]->node;
        set<int> res = computeNodeBranchDists(node, dad);
        for (i = res.begin(); i != res.end(); i++)
            nodeBranchDists[(*i) * nodeNum + dad->id] = nodeBranchDists[(dad->id) * nodeNum + (*i)] =
            nodeBranchDists[(*i) * nodeNum + node->id] + 1;
        // sanity check that all distances are filled
        for (int x = 0; x < nodeNum; x++)
            for (int y = 0; y < nodeNum; y++)
                if (x != y)
                    assert(nodeBranchDists[x * nodeNum + y] != 0);
                else
                    assert(nodeBranchDists[x * nodeNum + y] == 0);
        return res;
    }
    if (node->isLeaf()) {
        set<int> res;
        res.insert(node->id);
        return res;
    }
    assert(node->degree() == 3);
    Node* left = NULL, * right = NULL;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!left) left = (*it)->node; else right = (*it)->node;
    }
    set<int> resl = computeNodeBranchDists(left, node);
    set<int> resr = computeNodeBranchDists(right, node);
    for (i = resl.begin(); i != resl.end(); i++)
        nodeBranchDists[(*i) * nodeNum + node->id] = nodeBranchDists[(node->id) * nodeNum + (*i)] =
        nodeBranchDists[(*i) * nodeNum + left->id] + 1;
    for (i = resr.begin(); i != resr.end(); i++)
        nodeBranchDists[(*i) * nodeNum + node->id] = nodeBranchDists[(node->id) * nodeNum + (*i)] =
        nodeBranchDists[(*i) * nodeNum + right->id] + 1;
    for (i = resl.begin(); i != resl.end(); i++)
        for (j = resr.begin(); j != resr.end(); j++)
            nodeBranchDists[(*i) * nodeNum + (*j)] = nodeBranchDists[(*j) * nodeNum + (*i)] =
            nodeBranchDists[(*i) * nodeNum + node->id] + nodeBranchDists[(*j) * nodeNum + node->id];
    resl.insert(resr.begin(), resr.end());
    resl.insert(node->id);
    return resl;
}


/*
    b0: initial guess for the maximum
*/
double PhyloTree::approxOneBranch(PhyloNode* node, PhyloNode* dad, double b0) {
    double b_max, ddl, b1, b2, std, seqlen;
    double t1, t3, t5, t11, t18, t21, t26, t29, t30, t32, t44, t46, t48;
    double beps = 1 / DBL_MAX;

    /* TODO: insert call to get sequence length */
    seqlen = getAlnNSite();

    /* use a robust first order approximation to the variance */
    std = sqrt(b0 / seqlen);

    /* determine neighbour points */
    b1 = b0 - std;
    if (b1 <= 0) b1 = beps; /* only happens for b<=1 with small seq. len. */
    b2 = b0 + std;

    /* TODO: insert calls to log-likelihood function */
    PhyloNeighbor* dad_nei = (PhyloNeighbor*)(dad->findNeighbor(node));
    PhyloNeighbor* node_nei = (PhyloNeighbor*)(node->findNeighbor(dad));
    double old_len = dad_nei->length;
    dad_nei->length = node_nei->length = b0;
    double l0 = computeLikelihoodBranch(dad_nei, dad);
    dad_nei->length = node_nei->length = b1;
    double l1 = computeLikelihoodBranch(dad_nei, dad);
    dad_nei->length = node_nei->length = b2;
    double l2 = computeLikelihoodBranch(dad_nei, dad);
    dad_nei->length = node_nei->length = old_len;

    t1 = sqrt(b0);
    t3 = sqrt(b2);
    t5 = sqrt(b1);
    t11 = pow(-t1 * l2 + t3 * l0 + t5 * l2 + t1 * l1 - t5 * l0 - t3 * l1, 2.0);
    t18 = -b0 * l2 + b2 * l0 + b1 * l2 + b0 * l1 - b1 * l0 - b2 * l1;
    t21 = t1 - t5;
    t26 = -t1 * t3 + t1 * t5 + b2 - t5 * t3;
    t29 = t18 * t18;
    t30 = 1 / t11;
    t32 = sqrt(t29 * t30);
    ddl = -2.0 * t11 / t18 / t21 / t26 / t32;

    if (ddl > 0) {
        /* the analytic extremum is a minimum,
           so the maximum is at the lower bound */
        b_max = 0;
    }
    else {
        t44 = pow(-t1 * b2 + t5 * b2 - t5 * b0 + t3 * b0 - t3 * b1 + t1 * b1, 2.0);
        t46 = t21 * t21;
        t48 = t26 * t26;
        b_max = t29 * t44 / t46 / t48 * t30 / 4.0;
    }

    return(b_max);
}

void PhyloTree::approxAllBranches(PhyloNode* node, PhyloNode* dad) {
    if (!node) {
        node = (PhyloNode*)root;
    }

    if (dad) {
        PhyloNeighbor* node_dad_nei = (PhyloNeighbor*)node->findNeighbor(dad);
        PhyloNeighbor* dad_node_nei = (PhyloNeighbor*)dad->findNeighbor(node);
        double len = approxOneBranch(node, dad, dad_node_nei->length);
        node_dad_nei->length = len;
        dad_node_nei->length = len;
    }

    for (NeighborVec::iterator it = (node)->neighbors.begin(); it != (node)->neighbors.end(); it++)
        if ((*it)->node != (dad)) {
            approxAllBranches((PhyloNode*)(*it)->node, node);
        }
}

/*
 void PhyloTree::computeAllSubtreeDists(PhyloNode* node, PhyloNode* dad) {
 if (!node) {
 node = (PhyloNode*) root;
 }

 if (dad) {
 // This function compute all pairwise subtree distance between subtree rooted at dad and others

 computeSubtreeDists(node, dad);
 }

 FOR_NEIGHBOR_IT(node, dad, it) {

 computeAllSubtreeDists((PhyloNode*) (*it)->node, node);
 }
 }

 void PhyloTree::computeSubtreeDists(PhyloNode* node, PhyloNode* dad) {
 // if both nodes are leaf then it is trivial, just retrieve the values from the distance matrix
 if (dad->isLeaf() && node->isLeaf()) {
 string key = nodePair2String(dad, node);
 assert(dist_matrix);
 int nseq = aln->getNSeq();
 double dist = dist_matrix[dad->id * nseq + node->id];
 interSubtreeDistances.insert(StringDoubleMap::value_type(key, dist));
 } else if (!dad->isLeaf() && node->isLeaf()) {

 FOR_NEIGHBOR_IT(node, dad, it) {

 computeSubtreeDists(dad, (PhyloNode*) (*it)->node);
 }

 }
 }
 */

double PhyloTree::computeBayesianBranchLength(PhyloNeighbor* dad_branch, PhyloNode* dad) {
    double obsLen = 0.0;
    PhyloNode* node = (PhyloNode*)dad_branch->node;
    PhyloNeighbor* node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
    assert(node_branch);
    /*
     if (node->isLeaf() || dad->isLeaf()) {
     return -1.0;
     }*/
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(node_branch, node);
    // now combine likelihood at the branch
    int nstates = aln->num_states;
    int numCat = site_rate->getNRate();
    size_t block = numCat * nstates;
    size_t nptn = aln->size();
    size_t ptn;
    int cat, state;
    double* tmp_state_freq = new double[nstates];
    //computeLikelihoodBranchNaive(dad_branch, dad, NULL, tmp_ptn_rates);
    //double sum_rates = 0.0;
    //for (ptn = 0; ptn < nptn; ptn++)
    //    sum_rates += tmp_ptn_rates[ptn] * aln->at(ptn).frequency;
    //cout << "sum_rates = " << sum_rates << endl;

    model->getStateFrequency(tmp_state_freq);

    for (ptn = 0; ptn < nptn; ptn++) {
        // Compute the probability of each state for the current site
        double sum_prob1 = 0.0, sum_prob2 = 0.0;
        size_t offset = ptn * block;
        double* partial_lh_site = node_branch->partial_lh + (offset);
        double* partial_lh_child = dad_branch->partial_lh + (offset);
        for (state = 0; state < nstates; state++) {
            tmp_anscentral_state_prob1[state] = 0.0;
            tmp_anscentral_state_prob2[state] = 0.0;
            for (cat = 0; cat < numCat; cat++) {
                tmp_anscentral_state_prob1[state] += partial_lh_site[nstates * cat + state];
                tmp_anscentral_state_prob2[state] += partial_lh_child[nstates * cat + state];
            }
            tmp_anscentral_state_prob1[state] *= tmp_state_freq[state];
            tmp_anscentral_state_prob2[state] *= tmp_state_freq[state];
            sum_prob1 += tmp_anscentral_state_prob1[state];
            sum_prob2 += tmp_anscentral_state_prob2[state];
        }
        bool sameState = false;
        int state1 = 0, state2 = 0;
        double cutoff = 1.0 / nstates;
        for (state = 0; state < nstates; state++) {
            tmp_anscentral_state_prob1[state] /= sum_prob1;
            tmp_anscentral_state_prob2[state] /= sum_prob2;
            if (tmp_anscentral_state_prob1[state] > tmp_anscentral_state_prob1[state1])
                state1 = state;
            if (tmp_anscentral_state_prob2[state] > tmp_anscentral_state_prob2[state2])
                state2 = state;
            if (tmp_anscentral_state_prob1[state] > cutoff && tmp_anscentral_state_prob2[state] > cutoff)
                sameState = true;
        }
        sameState = sameState || (state1 == state2);
        if (!sameState) {
            obsLen += aln->at(ptn).frequency;
        }

    }
    obsLen /= getAlnNSite();
    if (obsLen < MIN_BRANCH_LEN)
        obsLen = MIN_BRANCH_LEN;
    delete[] tmp_state_freq;

    return obsLen;
}

double PhyloTree::correctBranchLengthF81(double observedBran, double alpha) {
    double H = 0.0;
    double correctedBranLen;
    for (int i = 0; i < model->num_states; i++) {
        H += model->state_freq[i] * (1 - model->state_freq[i]);
    }
    observedBran = 1.0 - observedBran / H;
    // no gamma
    if (observedBran <= 0.0)
        return MAX_BRANCH_LEN;

    if (alpha <= 0.0) {
        correctedBranLen = -H * log(observedBran);
    }
    else {
        //if (verbose_mode >= VB_MAX) cout << "alpha: " << alpha << endl;

        correctedBranLen = H * alpha * (pow(observedBran, -1 / alpha) - 1);
    }

    if (correctedBranLen < MIN_BRANCH_LEN)
        correctedBranLen = MIN_BRANCH_LEN;
    if (correctedBranLen > MAX_BRANCH_LEN)
        correctedBranLen = MAX_BRANCH_LEN;

    return correctedBranLen;
}

double PhyloTree::computeCorrectedBayesianBranchLength(PhyloNeighbor* dad_branch, PhyloNode* dad) {
    double observedBran = computeBayesianBranchLength(dad_branch, dad);
    return correctBranchLengthF81(observedBran, site_rate->getGammaShape());
}

void PhyloTree::computeAllBayesianBranchLengths(Node* node, Node* dad) {

    if (!node)
        node = root;

    FOR_NEIGHBOR_IT(node, dad, it) {
        double branch_length = computeBayesianBranchLength((PhyloNeighbor*)(*it), (PhyloNode*)node);
        (*it)->length = branch_length;
        // set the backward branch length
        (*it)->node->findNeighbor(node)->length = (*it)->length;
        computeAllBayesianBranchLengths((*it)->node, node);
    }
}

//double PhyloTree::computeLikelihoodBranchNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh, double *pattern_rate) {
double PhyloTree::computeLikelihoodBranchNaive(PhyloNeighbor* dad_branch, PhyloNode* dad, double* pattern_lh) {
    PhyloNode* node = (PhyloNode*)dad_branch->node;
    PhyloNeighbor* node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
    //assert(node_branch);
    //assert(!site_rate->isSiteSpecificRate() || !model->isSiteSpecificModel());
    if (!central_partial_lh)
        initializeAllPartialLh();
    // swap node and dad if dad is a leaf
    // NEW: swap if root_state is given
    if (node->isLeaf() || (node->name == ROOT_NAME && root_state != aln->STATE_UNKNOWN)) {
        PhyloNode* tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor* tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }

    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(node_branch, node);
    // now combine likelihood at the branch

    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    int ncat = site_rate->getNRate();
    double p_invar = site_rate->getPInvar();
    double p_var_cat = (1.0 - p_invar) / (double)ncat;
    int nstates = aln->num_states;
    size_t block = ncat * nstates;
    int trans_size = model->getTransMatrixSize();
    size_t ptn; // for big data size > 4GB memory required
    int cat, state1, state2;
    size_t nptn = aln->size() + model_factory->unobserved_ptns.size();
    size_t orig_nptn = aln->size();
    int discrete_cat = site_rate->getNDiscreteRate();
    double* trans_mat = new double[discrete_cat * trans_size];
    double* state_freq = new double[nstates];
    model->getStateFrequency(state_freq);

    if (!site_rate->isSiteSpecificRate())
        for (cat = 0; cat < discrete_cat; cat++) {
            //trans_mat[cat] = model->newTransMatrix();
            double* trans_cat = trans_mat + (cat * trans_size);
            model_factory->computeTransMatrixFreq(dad_branch->length * site_rate->getRate(cat), state_freq, trans_cat);
        }

    bool not_ptn_cat = (site_rate->getPtnCat(0) < 0);
    double prob_const = 0.0; // probability of unobserved const patterns
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh, prob_const) private(ptn, cat, state1, state2)
#endif
    for (ptn = 0; ptn < nptn; ptn++) {
        double lh_ptn = 0.0; // likelihood of the pattern
        int dad_state = 1000; // just something big enough
        int ptn_cat = site_rate->getPtnCat(ptn);
        if (dad->name == ROOT_NAME && root_state != aln->STATE_UNKNOWN) {
            dad_state = root_state;
        }
        else if (dad->isLeaf()) {
            if (ptn < orig_nptn)
                dad_state = (*aln)[ptn][dad->id];
            else
                dad_state = model_factory->unobserved_ptns[ptn - orig_nptn];
        }
        int dad_offset = dad_state * nstates;
        if (site_rate->isSiteSpecificRate()) {
            if (ptn < orig_nptn)
                model_factory->computeTransMatrixFreq(dad_branch->length * site_rate->getPtnRate(ptn), state_freq, trans_mat);
            else
                model_factory->computeTransMatrixFreq(dad_branch->length, state_freq, trans_mat);
        }
        for (cat = 0; cat < ncat; cat++) {
            double lh_cat = 0.0; // likelihood of the pattern's category
            size_t lh_offset = cat * nstates + ptn * block;
            double* partial_lh_site = node_branch->partial_lh + lh_offset;
            double* partial_lh_child = dad_branch->partial_lh + lh_offset;
            if (dad_state < nstates) { // single state
                // external node
                double* trans_state = trans_mat + ((not_ptn_cat ? cat : ptn_cat) * trans_size + dad_offset);
                if (model->isSiteSpecificModel() && ptn < nptn)
                    trans_state += (nstates * nstates * model->getPtnModelID(ptn));
                for (state2 = 0; state2 < nstates; state2++)
                    lh_cat += partial_lh_child[state2] * trans_state[state2];
            }
            else {
                // internal node, or external node but ambiguous character
                for (state1 = 0; state1 < nstates; state1++) {
                    double lh_state = 0.0; // likelihood of state1
                    double* trans_state = trans_mat + ((not_ptn_cat ? cat : ptn_cat) * trans_size + state1 * nstates);
                    if (model->isSiteSpecificModel() && ptn < nptn)
                        trans_state += (nstates * nstates * model->getPtnModelID(ptn));
                    for (state2 = 0; state2 < nstates; state2++)
                        lh_state += partial_lh_child[state2] * trans_state[state2];
                    lh_cat += lh_state * partial_lh_site[state1];
                }
            }
            lh_ptn += lh_cat;
            _pattern_lh_cat[ptn * ncat + cat] = lh_cat;
            //            if (pattern_rate)
            //                rate_ptn += lh_cat * site_rate->getRate(cat);
        }
        if (ptn < orig_nptn) {
            //			if (pattern_rate)
            //				pattern_rate[ptn] = rate_ptn / lh_ptn;
            lh_ptn *= p_var_cat;
            if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
                lh_ptn += p_invar * state_freq[(int)(*aln)[ptn][0]];
            }
            //#ifdef DEBUG
            if (lh_ptn <= 0.0)
                cout << "Negative likelihood: " << lh_ptn << " " << site_rate->getPtnRate(ptn) << endl;
            //#endif
            lh_ptn = log(lh_ptn);
            _pattern_lh[ptn] = lh_ptn;
            if (discard_saturated_site && site_rate->isSiteSpecificRate() && site_rate->getPtnRate(ptn) >= MAX_SITE_RATE)
                continue;
            tree_lh += lh_ptn * aln->at(ptn).frequency;
        }
        else {
            lh_ptn = lh_ptn * p_var_cat + p_invar * state_freq[(int)model_factory->unobserved_ptns[ptn - orig_nptn]];
            prob_const += lh_ptn;
        }

    }
    if (orig_nptn < nptn) {
        // ascertainment bias correction
        prob_const = log(1.0 - prob_const);
        for (ptn = 0; ptn < orig_nptn; ptn++)
            _pattern_lh[ptn] -= prob_const;
        tree_lh -= aln->getNSite() * prob_const;
    }
    if (pattern_lh)
        memmove(pattern_lh, _pattern_lh, aln->size() * sizeof(double));
    delete[] state_freq;
    delete[] trans_mat;
    //for (cat = ncat-1; cat >= 0; cat--)
    //delete trans_mat[cat];
    //delete state_freq;

    return tree_lh;
}

double PhyloTree::computeLikelihoodZeroBranch(PhyloNeighbor* dad_branch, PhyloNode* dad) {
    double lh_zero_branch;
    double saved_len = dad_branch->length;
    PhyloNeighbor* node_branch = (PhyloNeighbor*)dad_branch->node->findNeighbor(dad);
    dad_branch->length = 0.0;
    node_branch->length = 0.0;
    lh_zero_branch = computeLikelihoodBranch(dad_branch, dad);
    // restore branch length
    dad_branch->length = saved_len;
    node_branch->length = saved_len;

    return lh_zero_branch;
}

void PhyloTree::computePartialLikelihoodNaive(PhyloNeighbor* dad_branch, PhyloNode* dad) {
    // don't recompute the likelihood
    if (dad_branch->partial_lh_computed & 1)
        return;
    Node* node = dad_branch->node;
    size_t ptn, cat;
    int ncat = site_rate->getNRate();
    int nstates = aln->num_states;
    size_t block = nstates * site_rate->getNRate();
    int trans_size = model->getTransMatrixSize();
    size_t lh_size = (aln->size() + model_factory->unobserved_ptns.size()) * block;
    double* partial_lh_site;
    size_t nptn = aln->size() + model_factory->unobserved_ptns.size();
    size_t orig_nptn = aln->size();

    dad_branch->lh_scale_factor = 0.0;
    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));

    assert(dad_branch->partial_lh);
    //if (!dad_branch->partial_lh)
    //	dad_branch->partial_lh = newPartialLh();
    if (node->isLeaf() && dad) {
        /* external node */
        memset(dad_branch->partial_lh, 0, lh_size * sizeof(double));
        for (ptn = 0; ptn < nptn; ptn++) {
            char state;
            partial_lh_site = dad_branch->partial_lh + (ptn * block);
            if (node->name == ROOT_NAME) {
                state = aln->STATE_UNKNOWN;
            }
            else {
                assert(node->id < aln->getNSeq());
                if (ptn < orig_nptn)
                    state = (aln->at(ptn))[node->id];
                else // ascertainment bias correction
                    state = model_factory->unobserved_ptns[ptn - orig_nptn];
            }
            if (state < nstates) {
                for (cat = 0; cat < ncat; cat++)
                    partial_lh_site[cat * nstates + state] = 1.0;
            }
            else if (state == aln->STATE_UNKNOWN) {
                // fill all entries (also over rate category) with 1.0
                dad_branch->scale_num[ptn] = -1;
                for (int state2 = 0; state2 < block; state2++) {
                    partial_lh_site[state2] = 1.0;
                }
            }
            else if (aln->seq_type == SEQ_DNA) {
                // ambiguous character, for DNA, RNA
                state = state - (nstates - 1);
                for (int state2 = 0; state2 < nstates && state2 <= 6; state2++)
                    if (state & (1 << state2)) {
                        for (cat = 0; cat < ncat; cat++)
                            partial_lh_site[cat * nstates + state2] = 1.0;
                    }
            }
            else if (aln->seq_type == SEQ_PROTEIN) {
                // ambiguous character, for DNA, RNA
                state = state - (nstates);
                assert(state < 2);
                int state_map[2] = { 4 + 8,32 + 64 };
                for (int state2 = 0; state2 <= 6; state2++)
                    if (state_map[(int)state] & (1 << state2)) {
                        for (cat = 0; cat < ncat; cat++)
                            partial_lh_site[cat * nstates + state2] = 1.0;
                    }
            }
            else {
                outError("Internal error ", __func__);
            }
        }
    }
    else {
        /* internal node */
        int discrete_cat = site_rate->getNDiscreteRate();
        double* trans_mat = new double[discrete_cat * trans_size];
        //for (cat = 0; cat < discrete_cat; cat++) trans_mat[cat] = model->newTransMatrix();
        for (ptn = 0; ptn < lh_size; ptn++) {
            dad_branch->partial_lh[ptn] = 1.0;
        }
        for (ptn = 0; ptn < nptn; ptn++)
            dad_branch->scale_num[ptn] = -1;

        FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME) {
            computePartialLikelihoodNaive((PhyloNeighbor*)(*it), (PhyloNode*)node);

            dad_branch->lh_scale_factor += ((PhyloNeighbor*)(*it))->lh_scale_factor;

            if (!site_rate->isSiteSpecificRate())
                for (cat = 0; cat < discrete_cat; cat++)
                    model_factory->computeTransMatrix((*it)->length * site_rate->getRate(cat), trans_mat + (cat * trans_size));

            bool not_ptn_cat = (site_rate->getPtnCat(0) < 0);

            double sum_scale = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private(ptn, cat, partial_lh_site)
#endif
            for (ptn = 0; ptn < nptn; ptn++)
                if (((PhyloNeighbor*)(*it))->scale_num[ptn] >= 0) {
                    // avoid the case that all child partial likelihoods equal 1.0
                    //double *partial_lh_child = ((PhyloNeighbor*) (*it))->partial_lh + (ptn*block);
                    //if (partial_lh_child[0] < 0.0) continue;
                    //
                    if (dad_branch->scale_num[ptn] < 0) dad_branch->scale_num[ptn] = 0;
                    dad_branch->scale_num[ptn] += ((PhyloNeighbor*)(*it))->scale_num[ptn];
                    int ptn_cat = 0;
                    if (ptn < orig_nptn) {
                        ptn_cat = site_rate->getPtnCat(ptn);
                        if (site_rate->isSiteSpecificRate())
                            model_factory->computeTransMatrix((*it)->length * site_rate->getPtnRate(ptn), trans_mat);
                    }
                    else {
                        if (site_rate->isSiteSpecificRate())
                            model_factory->computeTransMatrix((*it)->length, trans_mat);
                    }
                    for (cat = 0; cat < ncat; cat++) {
                        size_t lh_offset = cat * nstates + ptn * block;
                        partial_lh_site = dad_branch->partial_lh + lh_offset;
                        double* partial_lh_child = ((PhyloNeighbor*)(*it))->partial_lh + lh_offset;
                        for (int state = 0; state < nstates; state++) {
                            double lh_child = 0.0;
                            double* trans_state = trans_mat + ((not_ptn_cat ? cat : ptn_cat) * trans_size + state * nstates);
                            if (model->isSiteSpecificModel() && ptn < orig_nptn)
                                trans_state += (nstates * nstates * model->getPtnModelID(ptn));
                            for (int state2 = 0; state2 < nstates; state2++)
                                lh_child += trans_state[state2] * partial_lh_child[state2];

                            if (!isfinite(lh_child))
                                outError("Numerical error with ", __func__);
                            partial_lh_site[state] *= lh_child;
                        }
                    }
                    // check if one should scale partial likelihoods
                    bool do_scale = true;
                    partial_lh_site = dad_branch->partial_lh + (ptn * block);
                    for (cat = 0; cat < block; cat++)
                        if (partial_lh_site[cat] > SCALING_THRESHOLD) {
                            do_scale = false;
                            break;
                        }
                    if (!do_scale) continue;
                    // now do the likelihood scaling
                    /*
                     double lh_max = partial_lh_site[0];
                     for (cat = 1; cat < block; cat++)
                     if (lh_max < partial_lh_site[cat]) lh_max = partial_lh_site[cat];
                     for (cat = 0; cat < block; cat++)
                     partial_lh_site[cat] /= lh_max;
                     dad_branch->lh_scale_factor += log(lh_max) * (*aln)[ptn].frequency;

                     */
                    for (cat = 0; cat < block; cat++)
                        partial_lh_site[cat] /= SCALING_THRESHOLD;
                    // unobserved const pattern will never have underflow
                    sum_scale += LOG_SCALING_THRESHOLD * (*aln)[ptn].frequency;
                    dad_branch->scale_num[ptn] += 1;

                    //                if (pattern_scale)
                    //                pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
                }
            dad_branch->lh_scale_factor += sum_scale;
        }
        delete[] trans_mat;
        //for (cat = ncat - 1; cat >= 0; cat--)
        //  delete [] trans_mat[cat];
    }
    dad_branch->partial_lh_computed |= 1;

}

double PhyloTree::computeLikelihoodDervNaive(PhyloNeighbor* dad_branch, PhyloNode* dad, double& df, double& ddf) {
    PhyloNode* node = (PhyloNode*)dad_branch->node;
    PhyloNeighbor* node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
    assert(node_branch);
    // swap node and dad if dad is a leaf
    // NEW: swap if root_state is given
    if (node->isLeaf() || (node->name == ROOT_NAME && root_state != aln->STATE_UNKNOWN)) {
        PhyloNode* tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor* tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodNaive(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodNaive(node_branch, node);

    // now combine likelihood at the branch
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    df = ddf = 0.0;
    int ncat = site_rate->getNRate();
    double p_invar = site_rate->getPInvar();
    double p_var_cat = (1.0 - p_invar) / (double)ncat;
    int nstates = aln->num_states;
    size_t block = ncat * nstates;
    int trans_size = model->getTransMatrixSize();
    size_t nptn = aln->size() + model_factory->unobserved_ptns.size();
    size_t orig_nptn = aln->size();
    size_t ptn, cat, state1, state2;

    int discrete_cat = site_rate->getNDiscreteRate();

    double* trans_mat = new double[discrete_cat * trans_size];
    double* trans_derv1 = new double[discrete_cat * trans_size];
    double* trans_derv2 = new double[discrete_cat * trans_size];
    double* state_freq = new double[nstates];
    model->getStateFrequency(state_freq);

    if (!site_rate->isSiteSpecificRate())
        for (cat = 0; cat < discrete_cat; cat++) {
            //trans_mat[cat] = model->newTransMatrix();
            double* trans_cat = trans_mat + (cat * trans_size);
            double* derv1_cat = trans_derv1 + (cat * trans_size);
            double* derv2_cat = trans_derv2 + (cat * trans_size);
            double rate_val = site_rate->getRate(cat);
            //double rate_sqr = rate_val * rate_val;
            model_factory->computeTransDervFreq(dad_branch->length, rate_val, state_freq, trans_cat, derv1_cat,
                derv2_cat);
            /*
             for (state1 = 0; state1 < nstates; state1++) {
             double *trans_mat_state = trans_cat + (state1 * nstates);
             double *trans_derv1_state = derv1_cat + (state1 * nstates);
             double *trans_derv2_state = derv2_cat + (state1 * nstates);

             for (state2 = 0; state2 < nstates; state2++) {
             trans_mat_state[state2] *= state_freq[state1];
             trans_derv1_state[state2] *= (state_freq[state1] * rate_val);
             trans_derv2_state[state2] *= (state_freq[state1] * rate_sqr);
             }
             }*/
        }

    bool not_ptn_cat = (site_rate->getPtnCat(0) < 0);
    double derv1_frac;
    double derv2_frac;

    double my_df = 0.0;
    double my_ddf = 0.0;
    double prob_const = 0.0, prob_const_derv1 = 0.0, prob_const_derv2 = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh, my_df, my_ddf, prob_const, prob_const_derv1, prob_const_derv2) private(ptn, cat, state1, state2, derv1_frac, derv2_frac)
#endif
    for (ptn = 0; ptn < nptn; ptn++) {
        int ptn_cat = site_rate->getPtnCat(ptn);
        if (discard_saturated_site && site_rate->isSiteSpecificRate() && nptn < orig_nptn && site_rate->getPtnRate(ptn) >= MAX_SITE_RATE)
            continue;
        double lh_ptn = 0.0; // likelihood of the pattern
        double lh_ptn_derv1 = 0.0;
        double lh_ptn_derv2 = 0.0;
        int dad_state = aln->STATE_UNKNOWN;

        if (dad->name == ROOT_NAME && root_state != aln->STATE_UNKNOWN)
            dad_state = root_state;
        else if (dad->isLeaf()) {
            if (ptn < orig_nptn)
                dad_state = (*aln)[ptn][dad->id];
            else
                dad_state = model_factory->unobserved_ptns[ptn - orig_nptn];
        }
        int dad_offset = dad_state * nstates;
        if (site_rate->isSiteSpecificRate()) {
            if (ptn < orig_nptn)
                model_factory->computeTransDervFreq(dad_branch->length, site_rate->getPtnRate(ptn), state_freq, trans_mat,
                    trans_derv1, trans_derv2);
            else
                model_factory->computeTransDervFreq(dad_branch->length, 1.0, state_freq, trans_mat,
                    trans_derv1, trans_derv2);
        }
        for (cat = 0; cat < ncat; cat++) {
            size_t lh_offset = cat * nstates + ptn * block;
            double* partial_lh_site = node_branch->partial_lh + lh_offset;
            double* partial_lh_child = dad_branch->partial_lh + lh_offset;
            if (dad_state < nstates) {
                // external node
                int cat2 = (not_ptn_cat ? cat : ptn_cat) * trans_size + dad_offset;
                if (model->isSiteSpecificModel() && ptn < orig_nptn)
                    cat2 += (nstates * nstates * model->getPtnModelID(ptn));
                double* trans_state = trans_mat + cat2;
                double* derv1_state = trans_derv1 + cat2;
                double* derv2_state = trans_derv2 + cat2;
                for (state2 = 0; state2 < nstates; state2++) {
                    lh_ptn += partial_lh_child[state2] * trans_state[state2];
                    lh_ptn_derv1 += partial_lh_child[state2] * derv1_state[state2];
                    lh_ptn_derv2 += partial_lh_child[state2] * derv2_state[state2];
                }
            }
            else {
                // internal node, or external node but ambiguous character
                for (state1 = 0; state1 < nstates; state1++) {
                    double lh_state = 0.0; // likelihood of state1
                    double lh_state_derv1 = 0.0;
                    double lh_state_derv2 = 0.0;
                    int cat2 = (not_ptn_cat ? cat : ptn_cat) * trans_size + state1 * nstates;
                    if (model->isSiteSpecificModel() && ptn < orig_nptn)
                        cat2 += (nstates * nstates * model->getPtnModelID(ptn));
                    double* trans_state = trans_mat + cat2;
                    double* derv1_state = trans_derv1 + cat2;
                    double* derv2_state = trans_derv2 + cat2;
                    for (state2 = 0; state2 < nstates; state2++) {
                        lh_state += partial_lh_child[state2] * trans_state[state2];
                        lh_state_derv1 += partial_lh_child[state2] * derv1_state[state2];
                        lh_state_derv2 += partial_lh_child[state2] * derv2_state[state2];
                    }
                    lh_ptn += lh_state * partial_lh_site[state1];
                    lh_ptn_derv1 += lh_state_derv1 * partial_lh_site[state1];
                    lh_ptn_derv2 += lh_state_derv2 * partial_lh_site[state1];
                }
            }
        }
        /*		if (p_invar > 0.0) {
         lh_ptn *= p_var_cat;
         lh_ptn_derv1 *= p_var_cat;
         lh_ptn_derv2 *= p_var_cat;
         if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
         lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
         }
         assert(lh_ptn > 0);
         double derv1_frac = lh_ptn_derv1 / lh_ptn;
         double derv2_frac = lh_ptn_derv2 / lh_ptn;
         tree_lh += log(lh_ptn) * (*aln)[ptn].frequency;
         df += derv1_frac * (*aln)[ptn].frequency;
         ddf += (derv2_frac - derv1_frac * derv1_frac) * (*aln)[ptn].frequency;
         } else {
         double derv1_frac = lh_ptn_derv1 / lh_ptn;
         double derv2_frac = lh_ptn_derv2 / lh_ptn;
         lh_ptn *= p_var_cat;
         assert(lh_ptn > 0);
         tree_lh += log(lh_ptn) * (*aln)[ptn].frequency;
         df += derv1_frac * (*aln)[ptn].frequency;
         ddf += (derv2_frac - derv1_frac * derv1_frac) * (*aln)[ptn].frequency;

         }
         */
         // Tung beo's trick
        if (lh_ptn <= 0) {
            cout << "Abnormal " << __func__;
            abort();
        }

        if (ptn < orig_nptn) {
            lh_ptn = lh_ptn * p_var_cat;
            if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
                lh_ptn += p_invar * state_freq[(int)(*aln)[ptn][0]];
            }
            double pad = p_var_cat / lh_ptn;
            if (std::isinf(pad)) {
                lh_ptn_derv1 *= p_var_cat;
                lh_ptn_derv2 *= p_var_cat;
                derv1_frac = lh_ptn_derv1 / lh_ptn;
                derv2_frac = lh_ptn_derv2 / lh_ptn;
            }
            else {
                derv1_frac = lh_ptn_derv1 * pad;
                derv2_frac = lh_ptn_derv2 * pad;
            }
            double freq = (*aln)[ptn].frequency;
            double tmp1 = derv1_frac * freq;
            double tmp2 = derv2_frac * freq;
            my_df += tmp1;
            my_ddf += tmp2 - tmp1 * derv1_frac;
            lh_ptn = log(lh_ptn);
            //cout << "lh_ptn = " << lh_ptn << endl;
            tree_lh += lh_ptn * freq;
            _pattern_lh[ptn] = lh_ptn;
            if (!isfinite(lh_ptn) || !isfinite(my_df) || !isfinite(my_ddf)) {
                cout << "Abnormal " << __func__;
                abort();
            }
        }
        else {
            lh_ptn = lh_ptn * p_var_cat + p_invar * state_freq[(int)model_factory->unobserved_ptns[ptn - orig_nptn]];
            prob_const += lh_ptn;
            prob_const_derv1 += lh_ptn_derv1 * p_var_cat;
            prob_const_derv2 += lh_ptn_derv2 * p_var_cat;
        }

    }
    if (orig_nptn < nptn) {
        // ascertainment bias correction
        prob_const = 1.0 - prob_const;
        derv1_frac = prob_const_derv1 / prob_const;
        derv2_frac = prob_const_derv2 / prob_const;
        int nsites = aln->getNSite();
        my_df += nsites * derv1_frac;
        my_ddf += nsites * (derv2_frac + derv1_frac * derv1_frac);
        prob_const = log(prob_const);
        tree_lh -= nsites * prob_const;
        for (ptn = 0; ptn < orig_nptn; ptn++)
            _pattern_lh[ptn] -= prob_const;
    }
    delete[] state_freq;
    delete[] trans_derv2;
    delete[] trans_derv1;
    delete[] trans_mat;
    //for (cat = ncat-1; cat >= 0; cat--)
    //delete trans_mat[cat];
    //delete state_freq;
    df = my_df;
    ddf = my_ddf;

    return tree_lh;
}

/****************************************************************************
 Branch length optimization by maximum likelihood
 ****************************************************************************/

double PhyloTree::computeFunction(double value) {
    current_it->length = value;
    current_it_back->length = value;

    return -computeLikelihoodBranch(current_it, (PhyloNode*)current_it_back->node);
}

double PhyloTree::computeFuncDerv(double value, double& df, double& ddf) {
    current_it->length = value;
    current_it_back->length = value;
    double lh;
    lh = -computeLikelihoodDerv(current_it, (PhyloNode*)current_it_back->node, df, ddf);

    df = -df;
    ddf = -ddf;

    return lh;
}

double PhyloTree::optimizeOneBranch(PhyloNode* node1, PhyloNode* node2, bool clearLH, int maxNRStep) {
    double negative_lh;
    current_it = (PhyloNeighbor*)node1->findNeighbor(node2);
    assert(current_it);
    current_it_back = (PhyloNeighbor*)node2->findNeighbor(node1);
    assert(current_it_back);
    if (params->maximum_parsimony) {
        current_it->partial_lh_computed = 0;
        current_it_back->partial_lh_computed = 0;
        double pars_score = computeParsimonyBranch((PhyloNeighbor*)current_it, (PhyloNode*)node1);
        return -pars_score;
    }
    double current_len = current_it->length;
    double ferror, optx;
    assert(current_len >= 0.0);
    theta_computed = false;
    if (optimize_by_newton) // Newton-Raphson method
        optx = minimizeNewton(MIN_BRANCH_LEN, current_len, MAX_BRANCH_LEN, TOL_BRANCH_LEN, negative_lh, maxNRStep);
    else
        // Brent method
        optx = minimizeOneDimen(MIN_BRANCH_LEN, current_len, MAX_BRANCH_LEN, TOL_BRANCH_LEN, &negative_lh, &ferror);

    current_it->length = optx;
    current_it_back->length = optx;
    //curScore = -negative_lh;

    if (clearLH && current_len != optx) {
        node1->clearReversePartialLh(node2);
        node2->clearReversePartialLh(node1);
    }

    return -negative_lh;
}

double PhyloTree::optimizeChildBranches(PhyloNode* node, PhyloNode* dad) {

    double tree_lh = 0.0;

    FOR_NEIGHBOR_IT(node, dad, it) {

        tree_lh = optimizeOneBranch((PhyloNode*)node, (PhyloNode*)(*it)->node);
    }
    return tree_lh;
}

void PhyloTree::optimizeAllBranchesLS(PhyloNode* node, PhyloNode* dad) {
    if (!node) {
        node = (PhyloNode*)root;
    }

    if (dad) {
        double lsBran = optimizeOneBranchLS(node, dad);
        PhyloNeighbor* node_dad_nei = (PhyloNeighbor*)node->findNeighbor(dad);
        PhyloNeighbor* dad_node_nei = (PhyloNeighbor*)dad->findNeighbor(node);
        node_dad_nei->length = lsBran;
        dad_node_nei->length = lsBran;
    }

    for (NeighborVec::iterator it = (node)->neighbors.begin(); it != (node)->neighbors.end(); it++)
        if ((*it)->node != (dad)) {
            optimizeAllBranchesLS((PhyloNode*)(*it)->node, node);
        }
}

double PhyloTree::optimizeAllBranches(PhyloNode* node, PhyloNode* dad, int maxNRStep) {
    //double tree_lh = optimizeChildBranches(node, dad);
    double tree_lh = -DBL_MAX;

    for (NeighborVec::iterator it = (node)->neighbors.begin(); it != (node)->neighbors.end(); it++)
        if ((*it)->node != (dad)) {
            //if (!(*it)->node->isLeaf())
            double new_tree_lh = optimizeAllBranches((PhyloNode*)(*it)->node, node, maxNRStep);
            /*
             if (new_tree_lh < tree_lh)
             cout << "Wrong " << __func__ << endl;
             */
            tree_lh = new_tree_lh;
        }
    if (dad)
        tree_lh = optimizeOneBranch(node, dad, true, maxNRStep); // BQM 2014-02-24: true was missing

    return tree_lh;
}

double PhyloTree::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {
    if (params->maximum_parsimony) {
        clearAllPartialLH();
        double score = -computeParsimony();
        return score;
    }
    if (verbose_mode >= VB_MAX)
        cout << "Optimizing branch lengths (max " << my_iterations << " loops)..." << endl;
    double tree_lh = computeLikelihood();
    if (verbose_mode >= VB_MAX) {
        cout << "Initial tree log-likelihood: " << tree_lh << endl;
    }
    //cout << tree_lh << endl;
    for (int i = 0; i < my_iterations; i++) {
        double new_tree_lh = optimizeAllBranches((PhyloNode*)root, NULL, maxNRStep);
        /*
         clearAllPartialLH();
         double new_tree_lh2 = computeLikelihood();
         if (fabs(new_tree_lh - new_tree_lh2) > TOL_LIKELIHOOD) {
         cout << "Wrong " << new_tree_lh <<" "<< new_tree_lh2 << endl;
         exit(1);
         }*/
        if (verbose_mode >= VB_MAX) {
            cout << "Likelihood after iteration " << i + 1 << " : ";
            cout << new_tree_lh << endl;
        }

        // CRITICAL BUG FIX: THIS GIVES WRONG LIKELIHOOD
        /*
        if (new_tree_lh <= tree_lh + tolerance)
            return (new_tree_lh > tree_lh) ? new_tree_lh : tree_lh;
        */
        // BQM comment for above: WHY DID I DO THIS??? this will make the loop never stop after my_iterations

//        assert(new_tree_lh >= tree_lh); // make sure that the new tree likelihood never decreases
        if (tree_lh <= new_tree_lh && new_tree_lh <= tree_lh + tolerance)
            return new_tree_lh;
        tree_lh = new_tree_lh;
    }
    return tree_lh;
}

/****************************************************************************
 Stepwise addition (greedy) by maximum likelihood
 ****************************************************************************/

double PhyloTree::addTaxonML(Node* added_node, Node*& target_node, Node*& target_dad, Node* node, Node* dad) {

    Neighbor* dad_nei = dad->findNeighbor(node);

    // now insert the new node in the middle of the branch node-dad
    double len = dad_nei->length;
    node->updateNeighbor(dad, added_node, len / 2.0);
    dad->updateNeighbor(node, added_node, len / 2.0);
    added_node->updateNeighbor((Node*)1, node, len / 2.0);
    added_node->updateNeighbor((Node*)2, dad, len / 2.0);
    // compute the likelihood
    clearAllPartialLH();
    double best_score = optimizeChildBranches((PhyloNode*)added_node);
    target_node = node;
    target_dad = dad;
    // remove the added node
    node->updateNeighbor(added_node, dad, len);
    dad->updateNeighbor(added_node, node, len);
    added_node->updateNeighbor(node, (Node*)1, len);
    added_node->updateNeighbor(dad, (Node*)2, len);

    // now tranverse the tree downwards

    FOR_NEIGHBOR_IT(node, dad, it) {
        Node* target_node2;
        Node* target_dad2;
        double score = addTaxonML(added_node, target_node2, target_dad2, (*it)->node, node);
        if (score > best_score) {

            best_score = score;
            target_node = target_node2;
            target_dad = target_dad2;
        }
    }
    return best_score;
}

void PhyloTree::growTreeML(Alignment* alignment) {

    cout << "Stepwise addition using ML..." << endl;
    aln = alignment;
    int size = aln->getNSeq();
    if (size < 3)
        outError(ERR_FEW_TAXA);

    root = newNode();
    Node* new_taxon;

    // create initial tree with 3 taxa
    for (leafNum = 0; leafNum < 3; leafNum++) {
        cout << "Add " << aln->getSeqName(leafNum) << " to the tree" << endl;
        new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
        root->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(root, 1.0);
    }
    root = findNodeID(0);
    optimizeAllBranches();

    // stepwise adding the next taxon
    for (leafNum = 3; leafNum < size; leafNum++) {

        cout << "Add " << aln->getSeqName(leafNum) << " to the tree" << endl;
        // allocate a new taxon and a new ajedcent internal node
        new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
        Node* added_node = newNode();
        added_node->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(added_node, 1.0);

        // preserve two neighbors
        added_node->addNeighbor((Node*)1, 1.0);
        added_node->addNeighbor((Node*)2, 1.0);

        Node* target_node = NULL;
        Node* target_dad = NULL;
        addTaxonML(added_node, target_node, target_dad, root->neighbors[0]->node, root);
        // now insert the new node in the middle of the branch node-dad
        double len = target_dad->findNeighbor(target_node)->length;
        target_node->updateNeighbor(target_dad, added_node, len / 2.0);
        target_dad->updateNeighbor(target_node, added_node, len / 2.0);
        added_node->updateNeighbor((Node*)1, target_node, len / 2.0);
        added_node->updateNeighbor((Node*)2, target_dad, len / 2.0);
        // compute the likelihood
        clearAllPartialLH();
        optimizeAllBranches();
        //optimizeNNI();
    }

    nodeNum = 2 * leafNum - 2;
}

/****************************************************************************
 Distance function
 ****************************************************************************/

double PhyloTree::computeDist(int seq1, int seq2, double initial_dist, double& d2l) {
    // if no model or site rate is specified, return JC distance
    if (initial_dist == 0.0) {
        if (params->compute_obs_dist)
            initial_dist = aln->computeObsDist(seq1, seq2);
        else
            initial_dist = aln->computeDist(seq1, seq2);
    }
    if (!model_factory || !site_rate)
        return initial_dist; // MANUEL: here no d2l is return

    // now optimize the distance based on the model and site rate
    AlignmentPairwise aln_pair(this, seq1, seq2);

    return aln_pair.optimizeDist(initial_dist, d2l);
}

double PhyloTree::computeDist(int seq1, int seq2, double initial_dist) {
    double var;
    return computeDist(seq1, seq2, initial_dist, var);
}

double PhyloTree::correctDist(double* dist_mat) {
    int i, j, k, pos;
    int n = aln->getNSeq();
    int nsqr = n * n;
    // use Floyd algorithm to find shortest path between all pairs of taxa
    for (k = 0; k < n; k++)
        for (i = 0, pos = 0; i < n; i++)
            for (j = 0; j < n; j++, pos++) {
                double tmp = dist_mat[i * n + k] + dist_mat[k * n + j];
                if (dist_mat[pos] > tmp)
                    dist_mat[pos] = tmp;
            }
    double longest_dist = 0.0;
    for (i = 0; i < nsqr; i++)
        if (dist_mat[i] > longest_dist)
            longest_dist = dist_mat[i];

    return longest_dist;
}

double PhyloTree::computeDist(double* dist_mat, double* var_mat) {
    int nseqs = aln->getNSeq();
    int pos = 0;
    int num_pairs = nseqs * (nseqs - 1) / 2;
    double longest_dist = 0.0;
    double d2l;
    int* row_id = new int[num_pairs];
    int* col_id = new int[num_pairs];

    row_id[0] = 0;
    col_id[0] = 1;
    for (pos = 1; pos < num_pairs; pos++) {
        row_id[pos] = row_id[pos - 1];
        col_id[pos] = col_id[pos - 1] + 1;
        if (col_id[pos] >= nseqs) {
            row_id[pos]++;
            col_id[pos] = row_id[pos] + 1;
        }
    }
    // compute the upper-triangle of distance matrix
#ifdef _OPENMP
#pragma omp parallel for private(pos)
#endif

    for (pos = 0; pos < num_pairs; pos++) {
        int seq1 = row_id[pos];
        int seq2 = col_id[pos];
        int sym_pos = seq1 * nseqs + seq2;
        dist_mat[sym_pos] = computeDist(seq1, seq2, dist_mat[sym_pos], d2l);
        if (params->ls_var_type == OLS)
            var_mat[sym_pos] = 1.0;
        else if (params->ls_var_type == WLS_PAUPLIN)
            var_mat[sym_pos] = 0.0;
        else if (params->ls_var_type == WLS_FIRST_TAYLOR)
            var_mat[sym_pos] = dist_mat[sym_pos];
        else if (params->ls_var_type == WLS_FITCH_MARGOLIASH)
            var_mat[sym_pos] = dist_mat[sym_pos] * dist_mat[sym_pos];
        else if (params->ls_var_type == WLS_SECOND_TAYLOR)
            var_mat[sym_pos] = -1.0 / d2l;
    }

    // copy upper-triangle into lower-triangle and set diagonal = 0
    for (int seq1 = 0; seq1 < nseqs; seq1++)
        for (int seq2 = 0; seq2 <= seq1; seq2++) {
            pos = seq1 * nseqs + seq2;
            if (seq1 == seq2) {
                dist_mat[pos] = 0.0;
                var_mat[pos] = 0.0;
            }
            else {
                dist_mat[pos] = dist_mat[seq2 * nseqs + seq1];
                var_mat[pos] = var_mat[seq2 * nseqs + seq1];
            }
            if (dist_mat[pos] > longest_dist)
                longest_dist = dist_mat[pos];
        }
    delete[] col_id;
    delete[] row_id;

    /*
     if (longest_dist > MAX_GENETIC_DIST * 0.99)
     outWarning("Some distances are saturated. Please check your alignment again");*/
     // NOTE: Bionj does handle long distances already (thanks Manuel)
     //return correctDist(dist_mat);
    return longest_dist;
}

double PhyloTree::computeDist(Params& params, Alignment* alignment, double*& dist_mat, double*& var_mat,
    string& dist_file) {
    this->params = &params;
    double longest_dist = 0.0;
    aln = alignment;
    dist_file = params.out_prefix;
    if (!model_factory) {
        if (params.compute_obs_dist)
            dist_file += ".obsdist";
        else
            //dist_file += ".jcdist"; // too many files, I decided to discard .jcdist
            dist_file += ".mldist";
    }
    else
        dist_file += ".mldist";

    if (!dist_mat) {
        dist_mat = new double[alignment->getNSeq() * alignment->getNSeq()];
        memset(dist_mat, 0, sizeof(double) * alignment->getNSeq() * alignment->getNSeq());
        var_mat = new double[alignment->getNSeq() * alignment->getNSeq()];
        // BUG!
        //memset(var_mat, 1, sizeof(double) * alignment->getNSeq() * alignment->getNSeq());
        int nseq = alignment->getNSeq();
        for (int i = 0; i < nseq; i++)
            for (int j = 0; j < nseq; j++)
                var_mat[i * nseq + j] = 1.0;
    }
    if (!params.dist_file) {
        longest_dist = computeDist(dist_mat, var_mat);
        alignment->printDist(dist_file.c_str(), dist_mat);
    }
    else {
        longest_dist = alignment->readDist(params.dist_file, dist_mat);
        dist_file = params.dist_file;
    }
    return longest_dist;
}

double PhyloTree::computeObsDist(double* dist_mat) {
    int nseqs = aln->getNSeq();
    int pos = 0;
    double longest_dist = 0.0;
    for (int seq1 = 0; seq1 < nseqs; seq1++)
        for (int seq2 = 0; seq2 < nseqs; seq2++, pos++) {
            if (seq1 == seq2)
                dist_mat[pos] = 0.0;
            else if (seq2 > seq1) {
                dist_mat[pos] = aln->computeObsDist(seq1, seq2);
            }
            else
                dist_mat[pos] = dist_mat[seq2 * nseqs + seq1];

            if (dist_mat[pos] > longest_dist)
                longest_dist = dist_mat[pos];
        }
    return longest_dist;
    //return correctDist(dist_mat);
}

double PhyloTree::computeObsDist(Params& params, Alignment* alignment, double*& dist_mat, string& dist_file) {
    double longest_dist = 0.0;
    aln = alignment;
    dist_file = params.out_prefix;
    dist_file += ".obsdist";

    if (!dist_mat) {
        dist_mat = new double[alignment->getNSeq() * alignment->getNSeq()];
        memset(dist_mat, 0, sizeof(double) * alignment->getNSeq() * alignment->getNSeq());
    }
    longest_dist = computeObsDist(dist_mat);
    alignment->printDist(dist_file.c_str(), dist_mat);

    return longest_dist;
}

/****************************************************************************
 compute BioNJ tree, a more accurate extension of Neighbor-Joining
 ****************************************************************************/

void PhyloTree::computeBioNJ(Params& params, Alignment* alignment, string& dist_file) {
    string bionj_file = params.out_prefix;
    bionj_file += ".bionj";
    cout << "Computing BIONJ tree..." << endl;
    BioNj bionj;
    bionj.create(dist_file.c_str(), bionj_file.c_str());
    bool my_rooted = false;
    bool non_empty_tree = (root != NULL);
    if (root)
        freeNode();
    readTree(bionj_file.c_str(), my_rooted);

    if (non_empty_tree)
        initializeAllPartialLh();
    setAlignment(alignment);
}

int PhyloTree::fixNegativeBranch(bool force, Node* node, Node* dad) {

    if (!node)
        node = root;
    int fixed = 0;

    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->length < 0.0 || force) { // negative branch length detected
            int branch_subst;
            int pars_score = computeParsimonyBranch((PhyloNeighbor*)(*it), (PhyloNode*)node, &branch_subst);
            // first compute the observed parsimony distance
            double branch_length = (branch_subst > 0) ? ((double)branch_subst / getAlnNSite()) : (1.0 / getAlnNSite());
            // now correct Juke-Cantor formula
            double z = (double)aln->num_states / (aln->num_states - 1);
            double x = 1.0 - (z * branch_length);
            if (x > 0) branch_length = -log(x) / z;
            if (branch_length < MIN_BRANCH_LEN)
                branch_length = MIN_BRANCH_LEN;
            //        if (verbose_mode >= VB_DEBUG)
            //        	cout << "Negative branch length " << (*it)->length << " was set to ";
                    //(*it)->length = fixed_length;
                    //(*it)->length = random_double()+0.1;
            (*it)->length = branch_length;
            if (verbose_mode >= VB_DEBUG)
                cout << (*it)->length << " parsimony = " << pars_score << endl;
            // set the backward branch length
            (*it)->node->findNeighbor(node)->length = (*it)->length;
            fixed++;
        }
        if ((*it)->length <= 0.0) {
            (*it)->length = MIN_BRANCH_LEN;
            (*it)->node->findNeighbor(node)->length = (*it)->length;
        }
        fixed += fixNegativeBranch(force, (*it)->node, node);
    }
    return fixed;
}

//int PhyloTree::assignRandomBranchLengths(bool force, Node *node, Node *dad) {
//
//    if (!node)
//        node = root;
//    int fixed = 0;
//
//    FOR_NEIGHBOR_IT(node, dad, it){
//		if ((*it)->length < 0.0 || force) { // negative branch length detected
//			if (verbose_mode >= VB_DEBUG)
//			cout << "Negative branch length " << (*it)->length << " was set to ";
//			(*it)->length = random_double() + 0.1;
//			if (verbose_mode >= VB_DEBUG)
//			cout << (*it)->length << endl;
//			// set the backward branch length
//			(*it)->node->findNeighbor(node)->length = (*it)->length;
//			fixed++;
//		}
//		if ((*it)->length <= 0.0) {
//			(*it)->length = 1e-6;
//			(*it)->node->findNeighbor(node)->length = (*it)->length;
//		}
//		fixed += assignRandomBranchLengths(force, (*it)->node, node);
//    }
//    return fixed;
//}

/****************************************************************************
 Nearest Neighbor Interchange by maximum likelihood
 ****************************************************************************/

void PhyloTree::doOneRandomNNI(Node* node1, Node* node2) {
    assert(!node1->isLeaf() && !node2->isLeaf());
    assert(node1->degree() == 3 && node2->degree() == 3);
    assert(node1->neighbors.size() == 3 && node2->neighbors.size() == 3);

    Neighbor* node1Nei = NULL;
    Neighbor* node2Nei = NULL;
    // randomly choose one neighbor from node1 and one neighbor from node2
    bool chooseNext = false;
    FOR_NEIGHBOR_IT(node1, node2, it) {
        if (chooseNext) {
            node1Nei = (*it);
            break;
        }
        int randNum = random_int(1);
        if (randNum == 0) {
            node1Nei = (*it);
            break;
        }
        else {
            chooseNext = true;
        }
    }
    chooseNext = false;
    FOR_NEIGHBOR_IT(node2, node1, it) {
        if (chooseNext) {
            node2Nei = (*it);
            break;
        }
        int randNum = random_int(1);
        if (randNum == 0) {
            node2Nei = (*it);
            break;
        }
        else {
            chooseNext = true;
        }
    }
    assert(node1Nei != NULL && node2Nei != NULL);

    NeighborVec::iterator node1NeiIt = node1->findNeighborIt(node1Nei->node);
    NeighborVec::iterator node2NeiIt = node2->findNeighborIt(node2Nei->node);
    assert(node1NeiIt != node1->neighbors.end());
    assert(node1NeiIt != node2->neighbors.end());

    node1->updateNeighbor(node1NeiIt, node2Nei);
    node2Nei->node->updateNeighbor(node2, node1);

    node2->updateNeighbor(node2NeiIt, node1Nei);
    node1Nei->node->updateNeighbor(node1, node2);
}

void PhyloTree::doNNI(NNIMove& move, bool clearLH) {
    PhyloNode* node1 = move.node1;
    PhyloNode* node2 = move.node2;
    NeighborVec::iterator node1Nei_it = move.node1Nei_it;
    NeighborVec::iterator node2Nei_it = move.node2Nei_it;
    Neighbor* node1Nei = *(node1Nei_it);
    Neighbor* node2Nei = *(node2Nei_it);

    // TODO MINH
    /*	Node *nodeA = node1Nei->node;
     Node *nodeB = node2Nei->node;

     NeighborVec::iterator nodeA_it = nodeA->findNeighborIt(node1);
     NeighborVec::iterator nodeB_it = nodeB->findNeighborIt(node2);
     Neighbor *nodeANei = *(nodeA_it);
     Neighbor *nodeBNei = *(nodeB_it);
     *node1Nei_it = node2Nei;
     *nodeB_it = nodeANei;
     *node2Nei_it = node1Nei;
     *nodeA_it = nodeBNei;*/
     // END TODO MINH
    assert(node1->degree() == 3 && node2->degree() == 3);
    // do the NNI swap
    node1->updateNeighbor(node1Nei_it, node2Nei);
    node2Nei->node->updateNeighbor(node2, node1);

    node2->updateNeighbor(node2Nei_it, node1Nei);
    node1Nei->node->updateNeighbor(node1, node2);

    // BQM check branch ID
    /*
     if (node1->findNeighbor(nodeB)->id != nodeB->findNeighbor(node1)->id) {
     cout << node1->findNeighbor(nodeB)->id << "<->" << nodeB->findNeighbor(node1)->id << endl;
     cout << node1->id << "," << nodeB->id << endl;
     outError("Wrong ID");
     }
     if (node2->findNeighbor(nodeA)->id != nodeA->findNeighbor(node2)->id) {
     cout << node2->findNeighbor(nodeA)->id << "<->" << nodeA->findNeighbor(node2)->id << endl;
     cout << node2->id << "," << nodeA->id << endl;
     outError("Wrong ID");
     }*/

    PhyloNeighbor* node12_it = (PhyloNeighbor*)node1->findNeighbor(node2); // return neighbor of node1 which points to node 2
    PhyloNeighbor* node21_it = (PhyloNeighbor*)node2->findNeighbor(node1); // return neighbor of node2 which points to node 1

    if (clearLH) {
        // clear partial likelihood vector
        node12_it->clearPartialLh();
        node21_it->clearPartialLh();

        node2->clearReversePartialLh(node1);
        node1->clearReversePartialLh(node2);
        //if (params->nni5Branches)
        //	clearAllPartialLH();
    }

    if (params->leastSquareNNI) {
        updateSubtreeDists(move);
    }
}

void PhyloTree::changeNNIBrans(NNIMove nnimove) {
    PhyloNode* node1 = nnimove.node1;
    PhyloNode* node2 = nnimove.node2;
    PhyloNeighbor* node1_node2_nei = (PhyloNeighbor*)node1->findNeighbor(node2);
    PhyloNeighbor* node2_node1_nei = (PhyloNeighbor*)node2->findNeighbor(node1);
    node1_node2_nei->length = nnimove.newLen[0];
    node2_node1_nei->length = nnimove.newLen[0];
    if (params->nni5) {
        int i = 1;
        Neighbor* nei;
        Neighbor* nei_back;
        NeighborVec::iterator it;
        FOR_NEIGHBOR(node1, node2, it)
        {
            nei = (*it)->node->findNeighbor(node1);
            nei_back = (node1)->findNeighbor((*it)->node);
            nei->length = nnimove.newLen[i];
            nei_back->length = nnimove.newLen[i];
            i++;
        }
        FOR_NEIGHBOR(node2, node1, it)
        {
            nei = (*it)->node->findNeighbor(node2);
            nei_back = (node2)->findNeighbor((*it)->node);
            nei->length = nnimove.newLen[i];
            nei_back->length = nnimove.newLen[i];
            i++;
        }
    }
}

NNIMove PhyloTree::getBestNNIForBran(PhyloNode* node1, PhyloNode* node2, NNIMove* nniMoves) {
    assert(!node1->isLeaf() && !node2->isLeaf());
    assert(node1->degree() == 3 && node2->degree() == 3);

    NeighborVec::iterator it;
    int IT_NUM = (params->nni5) ? 6 : 2;

    NeighborVec::iterator saved_it[6];
    int id = 0;

    saved_it[id++] = node1->findNeighborIt(node2);
    saved_it[id++] = node2->findNeighborIt(node1);

    if (params->nni5) {
        FOR_NEIGHBOR(node1, node2, it)
            saved_it[id++] = (*it)->node->findNeighborIt(node1);

        FOR_NEIGHBOR(node2, node1, it)
            saved_it[id++] = (*it)->node->findNeighborIt(node2);
    }
    assert(id == IT_NUM);

    Neighbor* saved_nei[6];
    // save Neighbor and allocate new Neighbor pointer
    size_t partial_lh_size = getPartialLhBytes() / sizeof(double);
    double* new_partial_lh = new double[IT_NUM * partial_lh_size + MEM_ALIGNMENT / sizeof(double)];
    size_t mem_shift = 0;
    if (((intptr_t)new_partial_lh) % MEM_ALIGNMENT != 0)
        mem_shift = (MEM_ALIGNMENT - (((intptr_t)new_partial_lh) % MEM_ALIGNMENT)) / sizeof(double);

    for (id = 0; id < IT_NUM; id++) {
        saved_nei[id] = (*saved_it[id]);
        *saved_it[id] = new PhyloNeighbor(saved_nei[id]->node, saved_nei[id]->length);
        ((PhyloNeighbor*)(*saved_it[id]))->partial_lh = new_partial_lh + id * partial_lh_size + mem_shift;
        ((PhyloNeighbor*)(*saved_it[id]))->scale_num = newScaleNum();
        if (params->maximum_parsimony)
            ((PhyloNeighbor*)(*saved_it[id]))->partial_pars = newPartialPars();
    }

    // get the Neighbor again since it is replaced for saving purpose
    PhyloNeighbor* node12_it = (PhyloNeighbor*)node1->findNeighbor(node2);
    PhyloNeighbor* node21_it = (PhyloNeighbor*)node2->findNeighbor(node1);

    int cnt;

    //NNIMove nniMoves[2];
    bool newNNIMoves = false;
    if (!nniMoves) {
        //   Initialize the 2 NNI moves
        newNNIMoves = true;
        nniMoves = new NNIMove[2];
        nniMoves[0].ptnlh = nniMoves[1].ptnlh = NULL;
        nniMoves[0].node1 = NULL;

    }

    if (nniMoves[0].node1) {
        // assuming that node1Nei_it and node2Nei_it is defined in nniMoves structure
        for (cnt = 0; cnt < 2; cnt++) {
            // sanity check
            if (!node1->findNeighbor((*nniMoves[cnt].node1Nei_it)->node)) outError(__func__);
            if (!node2->findNeighbor((*nniMoves[cnt].node2Nei_it)->node)) outError(__func__);
        }
    }
    else {
        FOR_NEIGHBOR_IT(node1, node2, node1_it) {
            cnt = 0;
            FOR_NEIGHBOR_IT(node2, node1, node2_it) {
                //   Initialize the 2 NNI moves
                nniMoves[cnt].node1Nei_it = node1_it;
                nniMoves[cnt].node2Nei_it = node2_it;
                cnt++;
            }
            break;
        }
    }

    // Initialize node1 and node2 in nniMoves
    nniMoves[0].node1 = nniMoves[1].node1 = node1;
    nniMoves[0].node2 = nniMoves[1].node2 = node2;

    double backupScore = curScore;

    for (cnt = 0; cnt < 2; cnt++) {
        // do the NNI swap
        NeighborVec::iterator node1_it = nniMoves[cnt].node1Nei_it;
        NeighborVec::iterator node2_it = nniMoves[cnt].node2Nei_it;
        Neighbor* node1_nei = *node1_it;
        Neighbor* node2_nei = *node2_it;

        node1->updateNeighbor(node1_it, node2_nei);
        node2_nei->node->updateNeighbor(node2, node1);

        node2->updateNeighbor(node2_it, node1_nei);
        node1_nei->node->updateNeighbor(node1, node2);

        // clear partial likelihood vector
        node12_it->clearPartialLh();
        node21_it->clearPartialLh();

        // compute the score of the swapped topology
        double score = optimizeOneBranch(node1, node2, false, NNI_MAX_NR_STEP);
        nniMoves[cnt].newLen[0] = node1->findNeighbor(node2)->length;

        int i = 1;
        if (params->nni5) {
            FOR_NEIGHBOR(node1, node2, it)
            {
                ((PhyloNeighbor*)(*it)->node->findNeighbor(node1))->clearPartialLh();
                score = optimizeOneBranch(node1, (PhyloNode*)(*it)->node, false, NNI_MAX_NR_STEP);
                nniMoves[cnt].newLen[i] = node1->findNeighbor((*it)->node)->length;
                i++;
            }

            node21_it->clearPartialLh();

            FOR_NEIGHBOR(node2, node1, it)
            {
                ((PhyloNeighbor*)(*it)->node->findNeighbor(node2))->clearPartialLh();
                score = optimizeOneBranch(node2, (PhyloNode*)(*it)->node, false, NNI_MAX_NR_STEP);
                //node2_lastnei = (PhyloNeighbor*) (*it);
                nniMoves[cnt].newLen[i] = node2->findNeighbor((*it)->node)->length;
                i++;
            }
            node12_it->clearPartialLh();
        }
        nniMoves[cnt].newloglh = score;
        // compute the pattern likelihoods if wanted
        if (nniMoves[cnt].ptnlh)
            computePatternLikelihood(nniMoves[cnt].ptnlh, &score);

        if (save_all_trees == 2) {
            saveCurrentTree(score); // BQM: for new bootstrap
        }

        // else, swap back, also recover the branch lengths
        node1->updateNeighbor(node1_it, node1_nei);
        node1_nei->node->updateNeighbor(node2, node1);
        node2->updateNeighbor(node2_it, node2_nei);
        node2_nei->node->updateNeighbor(node1, node2);
    }

    // restore the Neighbor*
    for (id = IT_NUM - 1; id >= 0; id--) {
        delete[]((PhyloNeighbor*)*saved_it[id])->scale_num;
        //delete[] ((PhyloNeighbor*) *saved_it[id])->partial_lh;
        if (params->maximum_parsimony) delete[]((PhyloNeighbor*)(*saved_it[id]))->partial_pars;
        if (*saved_it[id] == current_it) current_it = (PhyloNeighbor*)saved_nei[id];
        if (*saved_it[id] == current_it_back) current_it_back = (PhyloNeighbor*)saved_nei[id];

        delete (*saved_it[id]);
        (*saved_it[id]) = saved_nei[id];
    }
    delete[] new_partial_lh;

    // restore the length of 4 branches around node1, node2
    FOR_NEIGHBOR(node1, node2, it)
        (*it)->length = (*it)->node->findNeighbor(node1)->length;
    FOR_NEIGHBOR(node2, node1, it)
        (*it)->length = (*it)->node->findNeighbor(node2)->length;

    // restore curScore
    curScore = backupScore;

    NNIMove res;
    if (nniMoves[0].newloglh > nniMoves[1].newloglh) {
        res = nniMoves[0];
    }
    else {
        res = nniMoves[1];
    }
    if (newNNIMoves) {
        delete[] nniMoves;
    }
    return res;
}


/****************************************************************************
 Subtree Pruning and Regrafting by maximum likelihood
 ****************************************************************************/

double PhyloTree::optimizeSPR_old(double cur_score, PhyloNode* node, PhyloNode* dad) {
    if (!node)
        node = (PhyloNode*)root;
    PhyloNeighbor* dad1_nei = NULL;
    PhyloNeighbor* dad2_nei = NULL;
    PhyloNode* sibling1 = NULL;
    PhyloNode* sibling2 = NULL;
    double sibling1_len = 0.0, sibling2_len = 0.0;

    if (dad && !dad->isLeaf()) {

        assert(dad->degree() == 3);
        // assign the sibling of node, with respect to dad

        FOR_NEIGHBOR_DECLARE(dad, node, it) {
            if (!sibling1) {
                dad1_nei = (PhyloNeighbor*)(*it);
                sibling1 = (PhyloNode*)(*it)->node;
                sibling1_len = (*it)->length;
            }
            else {

                dad2_nei = (PhyloNeighbor*)(*it);
                sibling2 = (PhyloNode*)(*it)->node;
                sibling2_len = (*it)->length;
            }
        }
        // remove the subtree leading to node
        double sum_len = sibling1_len + sibling2_len;
        sibling1->updateNeighbor(dad, sibling2, sum_len);
        sibling2->updateNeighbor(dad, sibling1, sum_len);
        PhyloNeighbor* sibling1_nei = (PhyloNeighbor*)sibling1->findNeighbor(sibling2);
        PhyloNeighbor* sibling2_nei = (PhyloNeighbor*)sibling2->findNeighbor(sibling1);
        sibling1_nei->clearPartialLh();
        sibling2_nei->clearPartialLh();

        // now try to move the subtree to somewhere else
        vector<PhyloNeighbor*> spr_path;

        FOR_NEIGHBOR(sibling1, sibling2, it)
        {
            spr_path.push_back(sibling1_nei);
            double score = swapSPR_old(cur_score, 1, node, dad, sibling1, sibling2, (PhyloNode*)(*it)->node, sibling1,
                spr_path);
            // if likelihood score improves, return
            if (score > cur_score)

                return score;
            spr_path.pop_back();
        }

        FOR_NEIGHBOR(sibling2, sibling1, it)
        {
            spr_path.push_back(sibling2_nei);
            double score = swapSPR_old(cur_score, 1, node, dad, sibling1, sibling2, (PhyloNode*)(*it)->node, sibling2,
                spr_path);
            // if likelihood score improves, return
            if (score > cur_score)

                return score;
            spr_path.pop_back();
        }
        // if likelihood does not imporve, swap back
        sibling1->updateNeighbor(sibling2, dad, sibling1_len);
        sibling2->updateNeighbor(sibling1, dad, sibling2_len);
        dad1_nei->node = sibling1;
        dad1_nei->length = sibling1_len;
        dad2_nei->node = sibling2;
        dad2_nei->length = sibling2_len;
        clearAllPartialLH();
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        double score = optimizeSPR_old(cur_score, (PhyloNode*)(*it)->node, node);

        if (score > cur_score) return score;
    }
    return cur_score;
}

/**
 move the subtree (dad1-node1) to the branch (dad2-node2)
 */
double PhyloTree::swapSPR_old(double cur_score, int cur_depth, PhyloNode* node1, PhyloNode* dad1, PhyloNode* orig_node1,
    PhyloNode* orig_node2, PhyloNode* node2, PhyloNode* dad2, vector<PhyloNeighbor*>& spr_path) {
    PhyloNeighbor* node1_nei = (PhyloNeighbor*)node1->findNeighbor(dad1);
    PhyloNeighbor* dad1_nei = (PhyloNeighbor*)dad1->findNeighbor(node1);
    double node1_dad1_len = node1_nei->length;
    PhyloNeighbor* node2_nei = (PhyloNeighbor*)node2->findNeighbor(dad2);

    if (dad2) {
        // now, connect (node1-dad1) to the branch (node2-dad2)

        bool first = true;
        PhyloNeighbor* node2_nei = (PhyloNeighbor*)node2->findNeighbor(dad2);
        PhyloNeighbor* dad2_nei = (PhyloNeighbor*)dad2->findNeighbor(node2);
        double len2 = node2_nei->length;

        FOR_NEIGHBOR_IT(dad1, node1, it) {
            if (first) {
                (*it)->node = dad2;
                (*it)->length = len2 / 2;
                dad2->updateNeighbor(node2, dad1, len2 / 2);
                first = false;
            }
            else {
                (*it)->node = node2;
                (*it)->length = len2 / 2;
                node2->updateNeighbor(dad2, dad1, len2 / 2);
            }
            ((PhyloNeighbor*)(*it))->clearPartialLh();
        }
        node2_nei->clearPartialLh();
        dad2_nei->clearPartialLh();
        node1_nei->clearPartialLh();
        vector<PhyloNeighbor*>::iterator it2;
        for (it2 = spr_path.begin(); it2 != spr_path.end(); it2++)
            (*it2)->clearPartialLh();
        clearAllPartialLH();
        // optimize relevant branches
        double score;

        /* testing different branch optimization */
        score = optimizeOneBranch(node1, dad1);
        //score = optimizeOneBranch(dad2, dad1);
        //score = optimizeOneBranch(node2, dad1);
        /*
         PhyloNode *cur_node = dad2;
         for (int i = spr_path.size()-1; i >= 0; i--) {
         score = optimizeOneBranch(cur_node, (PhyloNode*)spr_path[i]->node);
         cur_node = (PhyloNode*)spr_path[i]->node;
         }
         */
         //score = optimizeAllBranches(dad1);
         // if score improves, return
        if (score > cur_score)
            return score;
        // else, swap back
        node2->updateNeighbor(dad1, dad2, len2);
        dad2->updateNeighbor(dad1, node2, len2);
        node2_nei->clearPartialLh();
        dad2_nei->clearPartialLh();
        node1_nei->length = node1_dad1_len;
        dad1_nei->length = node1_dad1_len;

        // add to candiate SPR moves
        spr_moves.add(node1, dad1, node2, dad2, score);
    }
    if (cur_depth >= spr_radius)

        return cur_score;
    spr_path.push_back(node2_nei);

    FOR_NEIGHBOR_IT(node2, dad2, it) {
        double score = swapSPR(cur_score, cur_depth + 1, node1, dad1, orig_node1, orig_node2, (PhyloNode*)(*it)->node, node2, spr_path);
        if (score > cur_score) return score;
    }
    spr_path.pop_back();

    return cur_score;

}

double PhyloTree::optimizeSPR(double cur_score, PhyloNode* node, PhyloNode* dad) {
    if (!node)
        node = (PhyloNode*)root;
    PhyloNeighbor* dad1_nei = NULL;
    PhyloNeighbor* dad2_nei = NULL;
    PhyloNode* sibling1 = NULL;
    PhyloNode* sibling2 = NULL;
    double sibling1_len = 0.0, sibling2_len = 0.0;

    if (dad && !dad->isLeaf()) {

        assert(dad->degree() == 3);
        // assign the sibling of node, with respect to dad

        FOR_NEIGHBOR_DECLARE(dad, node, it) {
            if (!sibling1) {
                dad1_nei = (PhyloNeighbor*)(*it);
                sibling1 = (PhyloNode*)(*it)->node;
                sibling1_len = (*it)->length;
            }
            else {

                dad2_nei = (PhyloNeighbor*)(*it);
                sibling2 = (PhyloNode*)(*it)->node;
                sibling2_len = (*it)->length;
            }
        }
        // remove the subtree leading to node
        double sum_len = sibling1_len + sibling2_len;
        sibling1->updateNeighbor(dad, sibling2, sum_len);
        sibling2->updateNeighbor(dad, sibling1, sum_len);
        PhyloNeighbor* sibling1_nei = (PhyloNeighbor*)sibling1->findNeighbor(sibling2);
        PhyloNeighbor* sibling2_nei = (PhyloNeighbor*)sibling2->findNeighbor(sibling1);
        // save partial likelihood
        double* sibling1_partial_lh = sibling1_nei->partial_lh;
        double* sibling2_partial_lh = sibling2_nei->partial_lh;
        sibling1_nei->partial_lh = newPartialLh();
        sibling2_nei->partial_lh = newPartialLh();
        sibling1_nei->clearPartialLh();
        sibling2_nei->clearPartialLh();

        // now try to move the subtree to somewhere else
        vector<PhyloNeighbor*> spr_path;

        FOR_NEIGHBOR(sibling1, sibling2, it)
        {
            spr_path.push_back(sibling1_nei);
            double score = swapSPR(cur_score, 1, node, dad, sibling1, sibling2, (PhyloNode*)(*it)->node, sibling1,
                spr_path);
            // if likelihood score improves, return
            if (score > cur_score) {
                cout << "cur_score = " << cur_score << endl;
                cout << "Found new BETTER SCORE by SPR: " << score << endl;

                return score;
            }
            spr_path.pop_back();
        }

        FOR_NEIGHBOR(sibling2, sibling1, it)
        {
            spr_path.push_back(sibling2_nei);
            double score = swapSPR(cur_score, 1, node, dad, sibling1, sibling2, (PhyloNode*)(*it)->node, sibling2,
                spr_path);
            // if likelihood score improves, return
            if (score > cur_score) {
                cout << "cur_score = " << cur_score << endl;
                cout << "Found new BETTER SCORE by SPR: " << score << endl;

                return score;
            }
            spr_path.pop_back();
        }
        // if likelihood does not imporve, swap back
        sibling1->updateNeighbor(sibling2, dad, sibling1_len);
        sibling2->updateNeighbor(sibling1, dad, sibling2_len);
        dad1_nei->node = sibling1;
        dad1_nei->length = sibling1_len;
        dad2_nei->node = sibling2;
        dad2_nei->length = sibling2_len;
        delete[] sibling1_nei->partial_lh;
        delete[] sibling2_nei->partial_lh;
        sibling1_nei->partial_lh = sibling1_partial_lh;
        sibling2_nei->partial_lh = sibling2_partial_lh;
        //clearAllPartialLH();

    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        double score = optimizeSPR(cur_score, (PhyloNode*)(*it)->node, node);

        if (score > cur_score) return score;
    }
    return cur_score;
}

/**
 move the subtree (dad1-node1) to the branch (dad2-node2)
 */
double PhyloTree::swapSPR(double cur_score, int cur_depth, PhyloNode* node1, PhyloNode* dad1, PhyloNode* orig_node1,
    PhyloNode* orig_node2, PhyloNode* node2, PhyloNode* dad2, vector<PhyloNeighbor*>& spr_path) {

    PhyloNeighbor* node1_nei = (PhyloNeighbor*)node1->findNeighbor(dad1);
    PhyloNeighbor* dad1_nei = (PhyloNeighbor*)dad1->findNeighbor(node1);
    double node1_dad1_len = node1_nei->length;
    PhyloNeighbor* node2_nei = (PhyloNeighbor*)node2->findNeighbor(dad2);
    PhyloNeighbor* dad2_nei = (PhyloNeighbor*)dad2->findNeighbor(node2);

    //double* node1dad1_lh_save = node1_nei->partial_lh;
    //double* dad1node1_lh_save = dad1_nei->partial_lh;
    //double node1dad1_scale = node1_nei->lh_scale_factor;
    //double dad1node1_scale = dad1_nei->lh_scale_factor;

    double* node2dad2_lh_save = node2_nei->partial_lh;
    double* dad2node2_lh_save = dad2_nei->partial_lh;
    double node2dad2_scale = node2_nei->lh_scale_factor;
    double dad2node_scale = dad2_nei->lh_scale_factor;

    double len2 = node2_nei->length;
    double newLen2 = sqrt(len2);

    if (dad2 && cur_depth >= SPR_DEPTH) {
        // now, connect (node1-dad1) to the branch (node2-dad2)

        bool first = true;
        //PhyloNeighbor *node2_nei = (PhyloNeighbor*) node2->findNeighbor(dad2);
        //PhyloNeighbor *dad2_nei = (PhyloNeighbor*) dad2->findNeighbor(node2);
        //double len2 = node2_nei->length;

        FOR_NEIGHBOR_IT(dad1, node1, it) {
            // Finding new 2 neighbors for dad1 that are not node1
            if (first) {
                (*it)->node = dad2;
                //(*it)->length = len2 / 2;
                (*it)->length = newLen2;
                dad2->updateNeighbor(node2, dad1, newLen2);
                first = false;
            }
            else {
                (*it)->node = node2;
                (*it)->length = newLen2;
                node2->updateNeighbor(dad2, dad1, newLen2);
            }
            // clear all partial likelihood leading from
            // dad1 to the new neighbors
            ((PhyloNeighbor*)(*it))->clearPartialLh();
        }

        // clear partial likelihood from node2 to dad1
        node2_nei->clearPartialLh();
        // clear partial likelihood from dad2 to dad1
        dad2_nei->clearPartialLh();
        // clear partial likelihood from dad1 to node1
        node1_nei->clearPartialLh();

        // set new legnth as suggested by Alexis
        node1_nei->length = 0.9;
        dad1_nei->length = 0.9;

        //Save the partial likelihood from the removal point to the insertion point
        vector<PhyloNeighbor*>::iterator it2;
        vector<double*> saved_partial_lhs(spr_path.size());
        for (it2 = spr_path.begin(); it2 != spr_path.end(); it2++) {
            saved_partial_lhs.push_back((*it2)->partial_lh);
            (*it2)->partial_lh = newPartialLh();
            (*it2)->clearPartialLh();
        }

        // optimize relevant branches
        double score;

        /* testing different branch optimization */
        score = optimizeOneBranch(node1, dad1);
        score = optimizeOneBranch(dad2, dad1);
        score = optimizeOneBranch(node2, dad1);
        score = optimizeOneBranch(orig_node1, orig_node2);

        /*
         PhyloNode *cur_node = dad2;
         for (int i = spr_path.size()-1; i >= 0; i--) {
         score = optimizeOneBranch(cur_node, (PhyloNode*)spr_path[i]->node);
         cur_node = (PhyloNode*)spr_path[i]->node;
         }
         */
         //score = optimizeAllBranches(dad1);
         // if score improves, return
        if (score > cur_score) {
            cout << score << endl;
            return score;
        }

        // else, swap back
        node2->updateNeighbor(dad1, dad2, len2);
        dad2->updateNeighbor(dad1, node2, len2);
        //node2_nei->clearPartialLh();
        //dad2_nei->clearPartialLh();
        // restore partial likelihood vectors
        node2_nei->partial_lh = node2dad2_lh_save;
        node2_nei->lh_scale_factor = node2dad2_scale;
        dad2_nei->partial_lh = dad2node2_lh_save;
        dad2_nei->lh_scale_factor = dad2node_scale;
        node2_nei->length = len2;
        dad2_nei->length = len2;
        node1_nei->length = node1_dad1_len;
        dad1_nei->length = node1_dad1_len;
        int index = 0;
        for (it2 = spr_path.begin(); it2 != spr_path.end(); it2++) {
            delete[](*it2)->partial_lh;
            (*it2)->partial_lh = saved_partial_lhs.at(index);
            (*it2)->unclearPartialLh();
            index++;
        }

        // add to candiate SPR moves
        // Tung : why adding negative SPR move ?
        spr_moves.add(node1, dad1, node2, dad2, score);
    }
    if (cur_depth >= spr_radius)

        return cur_score;
    spr_path.push_back(node2_nei);

    FOR_NEIGHBOR_IT(node2, dad2, it) {
        double score = swapSPR(cur_score, cur_depth + 1, node1, dad1, orig_node1, orig_node2, (PhyloNode*)(*it)->node, node2, spr_path);
        if (score > cur_score) return score;
    }
    spr_path.pop_back();

    return cur_score;
}

double PhyloTree::assessSPRMove(double cur_score, const SPRMove& spr) {

    PhyloNode* dad = spr.prune_dad;
    PhyloNode* node = spr.prune_node;
    PhyloNode* dad2 = spr.regraft_dad;
    PhyloNode* node2 = spr.regraft_node;

    PhyloNeighbor* dad_nei1 = NULL;
    PhyloNeighbor* dad_nei2 = NULL;
    PhyloNode* sibling1 = NULL;
    PhyloNode* sibling2 = NULL;
    double sibling1_len = 0.0, sibling2_len = 0.0;

    PhyloNeighbor* node1_nei = (PhyloNeighbor*)node->findNeighbor(dad);
    PhyloNeighbor* dad1_nei = (PhyloNeighbor*)dad->findNeighbor(node);
    double node1_dad1_len = node1_nei->length;

    // assign the sibling of node, with respect to dad

    FOR_NEIGHBOR_DECLARE(dad, node, it) {
        if (!sibling1) {
            dad_nei1 = (PhyloNeighbor*)(*it);
            sibling1 = (PhyloNode*)(*it)->node;
            sibling1_len = (*it)->length;
        }
        else {

            dad_nei2 = (PhyloNeighbor*)(*it);
            sibling2 = (PhyloNode*)(*it)->node;
            sibling2_len = (*it)->length;
        }
    }
    // remove the subtree leading to node
    double sum_len = sibling1_len + sibling2_len;
    sibling1->updateNeighbor(dad, sibling2, sum_len);
    sibling2->updateNeighbor(dad, sibling1, sum_len);
    // now try to move the subtree to somewhere else

    bool first = true;
    PhyloNeighbor* node2_nei = (PhyloNeighbor*)node2->findNeighbor(dad2);
    //PhyloNeighbor *dad2_nei = (PhyloNeighbor*) dad2->findNeighbor(node2);
    double len2 = node2_nei->length;

    FOR_NEIGHBOR(dad, node, it)
    {
        if (first) {
            (*it)->node = dad2;
            (*it)->length = len2 / 2;
            dad2->updateNeighbor(node2, dad, len2 / 2);
            first = false;
        }
        else {
            (*it)->node = node2;
            (*it)->length = len2 / 2;
            node2->updateNeighbor(dad2, dad, len2 / 2);
        }
        ((PhyloNeighbor*)(*it))->clearPartialLh();
    }

    clearAllPartialLH();
    // optimize branches
    double score;
    score = optimizeAllBranches(dad);

    // if score improves, return
    if (score > cur_score)
        return score;
    // else, swap back
    node2->updateNeighbor(dad, dad2, len2);
    dad2->updateNeighbor(dad, node2, len2);

    node1_nei->length = node1_dad1_len;
    dad1_nei->length = node1_dad1_len;

    sibling1->updateNeighbor(sibling2, dad, sibling1_len);
    sibling2->updateNeighbor(sibling1, dad, sibling2_len);
    dad_nei1->node = sibling1;
    dad_nei1->length = sibling1_len;
    dad_nei2->node = sibling2;
    dad_nei2->length = sibling2_len;
    clearAllPartialLH();

    return cur_score;

}

double PhyloTree::optimizeSPR() {
    double cur_score = computeLikelihood();
    //spr_radius = leafNum / 5;
    spr_radius = 10;
    for (int i = 0; i < 100; i++) {
        cout << "i = " << i << endl;
        spr_moves.clear();
        double score = optimizeSPR_old(cur_score, (PhyloNode*)root->neighbors[0]->node);
        clearAllPartialLH();
        // why this?
        if (score <= cur_score) {
            for (SPRMoves::iterator it = spr_moves.begin(); it != spr_moves.end(); it++) {
                //cout << (*it).score << endl;
                score = assessSPRMove(cur_score, *it);
                // if likelihood score improves, apply to SPR
                if (score > cur_score)
                    break;
            }
            if (score <= cur_score)
                break;
        }
        else {

            cur_score = optimizeAllBranches();
            cout << "SPR " << i + 1 << " : " << cur_score << endl;
            cur_score = score;
        }
    }
    return cur_score;
    //return optimizeAllBranches();
}

double PhyloTree::optimizeSPRBranches() {
    cout << "Search with Subtree Pruning and Regrafting (SPR) using ML..." << endl;
    double cur_score = computeLikelihood();
    for (int i = 0; i < 100; i++) {
        double score = optimizeSPR();
        if (score <= cur_score + TOL_LIKELIHOOD)

            break;
        cur_score = score;
    }
    return cur_score;
}

void PhyloTree::pruneSubtree(PhyloNode* node, PhyloNode* dad, PruningInfo& info) {

    bool first = true;
    info.node = node;
    info.dad = dad;

    FOR_NEIGHBOR_IT(dad, node, it) {
        if (first) {
            info.dad_it_left = it;
            info.dad_nei_left = (*it);
            info.dad_lh_left = ((PhyloNeighbor*)(*it))->partial_lh;
            info.left_node = (*it)->node;
            info.left_len = (*it)->length;
            first = false;
        }
        else {

            info.dad_it_right = it;
            info.dad_nei_right = (*it);
            info.dad_lh_right = ((PhyloNeighbor*)(*it))->partial_lh;
            info.right_node = (*it)->node;
            info.right_len = (*it)->length;
        }
    }
    info.left_it = info.left_node->findNeighborIt(dad);
    info.right_it = info.right_node->findNeighborIt(dad);
    info.left_nei = (*info.left_it);
    info.right_nei = (*info.right_it);

    info.left_node->updateNeighbor(info.left_it, info.dad_nei_right);
    info.right_node->updateNeighbor(info.right_it, info.dad_nei_left);
    ((PhyloNeighbor*)info.dad_nei_right)->partial_lh = newPartialLh();
    ((PhyloNeighbor*)info.dad_nei_left)->partial_lh = newPartialLh();
}

void PhyloTree::regraftSubtree(PruningInfo& info, PhyloNode* in_node, PhyloNode* in_dad) {

    NeighborVec::iterator in_node_it = in_node->findNeighborIt(in_dad);
    NeighborVec::iterator in_dad_it = in_dad->findNeighborIt(in_node);
    Neighbor* in_dad_nei = (*in_dad_it);
    Neighbor* in_node_nei = (*in_node_it);
    //double in_len = in_dad_nei->length;
    info.dad->updateNeighbor(info.dad_it_right, in_dad_nei);
    info.dad->updateNeighbor(info.dad_it_left, in_node_nei);
    // SOMETHING NEED TO BE DONE
    //in_dad->updateNeighbor(in_dad_it,

}

/****************************************************************************
 Approximate Likelihood Ratio Test with SH-like interpretation
 ****************************************************************************/

 /*void PhyloTree::computeNNIPatternLh(double cur_lh, double &lh2, double *pattern_lh2, double &lh3, double *pattern_lh3,
         PhyloNode *node1, PhyloNode *node2) {

     assert(node1->degree() == 3 && node2->degree() == 3);

     // recompute pattern scaling factors if necessary
     PhyloNeighbor *node12_it = (PhyloNeighbor*) node1->findNeighbor(node2);
     PhyloNeighbor *node21_it = (PhyloNeighbor*) node2->findNeighbor(node1);
     NeighborVec::iterator it;
     const int IT_NUM = 6;

     NeighborVec::iterator saved_it[IT_NUM];
     int id = 0;

     FOR_NEIGHBOR(node1, node2, it)
     {
         saved_it[id++] = (*it)->node->findNeighborIt(node1);
     } else {

         saved_it[id++] = it;
     }

     FOR_NEIGHBOR(node2, node1, it)
     {
         saved_it[id++] = (*it)->node->findNeighborIt(node2);
     } else {
         saved_it[id++] = it;
     }
     assert(id == IT_NUM);

     Neighbor * saved_nei[IT_NUM];
     // save Neighbor and allocate new Neighbor pointer
     for (id = 0; id < IT_NUM; id++) {
         saved_nei[id] = (*saved_it[id]);
         // NOTE BUG DOWN HERE!
         *saved_it[id] = new PhyloNeighbor(saved_nei[id]->node, saved_nei[id]->length); // BUG for PhyloSuperTree!
         ((PhyloNeighbor*) (*saved_it[id]))->partial_lh = newPartialLh();
         ((PhyloNeighbor*) (*saved_it[id]))->scale_num = newScaleNum();
     }

     // get the Neighbor again since it is replaced for saving purpose
     node12_it = (PhyloNeighbor*) node1->findNeighbor(node2);
     node21_it = (PhyloNeighbor*) node2->findNeighbor(node1);

     // PhyloNeighbor *node2_lastnei = NULL;

     // save the first found neighbor of node 1 (excluding node2) in node1_it
     FOR_NEIGHBOR_DECLARE(node1, node2, node1_it)

         break;
     Neighbor *node1_nei = *node1_it;

     bool first = true;

     FOR_NEIGHBOR_IT(node2, node1, node2_it) {
         // do the NNI swap
         Neighbor *node2_nei = *node2_it;
         node1->updateNeighbor(node1_it, node2_nei);
         node2_nei->node->updateNeighbor(node2, node1);

         node2->updateNeighbor(node2_it, node1_nei);
         node1_nei->node->updateNeighbor(node1, node2);

         // re-optimize five adjacent branches
         double old_score = -INFINITY, new_score = old_score;

         // clear partial likelihood vector
         node12_it->clearPartialLh();
         node21_it->clearPartialLh();
         int i;
         for (i = 0; i < 2; i++) {

             new_score = optimizeOneBranch(node1, node2, false);

             FOR_NEIGHBOR(node1, node2, it) {
                 //for (id = 0; id < IT_NUM; id++)
                 //((PhyloNeighbor*)(*saved_it[id]))->clearPartialLh();
                 ((PhyloNeighbor*) (*it)->node->findNeighbor(node1))->clearPartialLh();
                 new_score = optimizeOneBranch(node1, (PhyloNode*) (*it)->node, false);
             }

             node21_it->clearPartialLh();

             FOR_NEIGHBOR(node2, node1, it) {
                 //for (id = 0; id < IT_NUM; id++)
                 //((PhyloNeighbor*)(*saved_it[id]))->clearPartialLh();
                 ((PhyloNeighbor*) (*it)->node->findNeighbor(node2))->clearPartialLh();
                 new_score = optimizeOneBranch(node2, (PhyloNode*) (*it)->node, false);
                 //node2_lastnei = (PhyloNeighbor*) (*it);
             }
             node12_it->clearPartialLh();
             if (new_score < old_score + TOL_LIKELIHOOD) break;
             old_score = new_score;
         }
         saveCurrentTree(new_score); // BQM: for new bootstrap

         //new_score = optimizeOneBranch(node1, node2, false);
         if (new_score > cur_lh + TOL_LIKELIHOOD)
         cout << "Alternative NNI shows better likelihood " << new_score << " > " << cur_lh << endl;
         double *result_lh;
         if (first) {
             result_lh = pattern_lh2;
             lh2 = new_score;
         } else {
             result_lh = pattern_lh3;
             lh3 = new_score;
         }
         old_score = new_score;
         computePatternLikelihood(result_lh);
         // swap back and recover the branch lengths
         node1->updateNeighbor(node1_it, node1_nei);
         node1_nei->node->updateNeighbor(node2, node1);
         node2->updateNeighbor(node2_it, node2_nei);
         node2_nei->node->updateNeighbor(node1, node2);
         first = false;
     }

 // restore the Neighbor*
     for (id = 0; id < IT_NUM; id++) {

         delete[] ((PhyloNeighbor*) *saved_it[id])->scale_num;
         delete[] ((PhyloNeighbor*) *saved_it[id])->partial_lh;
         delete (*saved_it[id]);
         (*saved_it[id]) = saved_nei[id];
     }

     // restore the length of 4 branches around node1, node2
     FOR_NEIGHBOR(node1, node2, it)
         (*it)->length = (*it)->node->findNeighbor(node1)->length;
     FOR_NEIGHBOR(node2, node1, it)
         (*it)->length = (*it)->node->findNeighbor(node2)->length;
 }*/

void PhyloTree::computeNNIPatternLh(double cur_lh, double& lh2, double* pattern_lh2, double& lh3, double* pattern_lh3,
    PhyloNode* node1, PhyloNode* node2) {
    NNIMove nniMoves[2];
    nniMoves[0].ptnlh = pattern_lh2;
    nniMoves[1].ptnlh = pattern_lh3;
    bool nni5 = params->nni5;
    params->nni5 = true; // always optimize 5 branches for accurate SH-aLRT
    nniMoves[0].node1 = nniMoves[1].node1 = NULL;
    nniMoves[0].node2 = nniMoves[1].node2 = NULL;
    getBestNNIForBran(node1, node2, nniMoves);
    params->nni5 = nni5;
    lh2 = nniMoves[0].newloglh;
    lh3 = nniMoves[1].newloglh;
    if (max(lh2, lh3) > cur_lh + TOL_LIKELIHOOD)
        cout << "Alternative NNI shows better log-likelihood " << max(lh2, lh3) << " > " << cur_lh << endl;
}

void PhyloTree::resampleLh(double** pat_lh, double* lh_new) {
    //int nsite = getAlnNSite();
    int nptn = getAlnNPattern();
    memset(lh_new, 0, sizeof(double) * 3);
    int i;
    IntVector boot_freq;
    aln->createBootstrapAlignment(boot_freq, params->bootstrap_spec);
    for (i = 0; i < nptn; i++) {

        lh_new[0] += boot_freq[i] * pat_lh[0][i];
        lh_new[1] += boot_freq[i] * pat_lh[1][i];
        lh_new[2] += boot_freq[i] * pat_lh[2][i];
    }
}

// Implementation of testBranch follows Guindon et al. (2010)

double PhyloTree::testOneBranch(double best_score, double* pattern_lh, int reps, int lbp_reps, PhyloNode* node1,
    PhyloNode* node2, double& lbp_support) {
    const int NUM_NNI = 3;
    double lh[NUM_NNI];
    double* pat_lh[NUM_NNI];
    lh[0] = best_score;
    pat_lh[0] = pattern_lh;
    pat_lh[1] = new double[getAlnNPattern()];
    pat_lh[2] = new double[getAlnNPattern()];
    computeNNIPatternLh(best_score, lh[1], pat_lh[1], lh[2], pat_lh[2], node1, node2);
    double aLRT;
    if (lh[1] > lh[2])
        aLRT = (lh[0] - lh[1]);
    else
        aLRT = (lh[0] - lh[2]);

    int support = 0;

    lbp_support = 0.0;
    int times = max(reps, lbp_reps);

    for (int i = 0; i < times; i++) {
        double lh_new[NUM_NNI];
        // resampling estimated log-likelihood (RELL)
        resampleLh(pat_lh, lh_new);
        if (lh_new[0] > lh_new[1] && lh_new[0] > lh_new[2])
            lbp_support += 1.0;
        double cs[NUM_NNI], cs_best, cs_2nd_best;
        cs[0] = lh_new[0] - lh[0];
        cs[1] = lh_new[1] - lh[1];
        cs[2] = lh_new[2] - lh[2];
        if (cs[0] >= cs[1] && cs[0] >= cs[2]) {
            cs_best = cs[0];
            if (cs[1] > cs[2])
                cs_2nd_best = cs[1];
            else
                cs_2nd_best = cs[2];
        }
        else if (cs[1] >= cs[2]) {
            cs_best = cs[1];
            if (cs[0] > cs[2])
                cs_2nd_best = cs[0];
            else
                cs_2nd_best = cs[2];
        }
        else {
            cs_best = cs[2];
            if (cs[0] > cs[1])
                cs_2nd_best = cs[0];
            else
                cs_2nd_best = cs[1];
        }
        if (aLRT > (cs_best - cs_2nd_best) + 0.05)
            support++;
    }
    delete[] pat_lh[2];
    delete[] pat_lh[1];
    lbp_support /= times;

    return ((double)support) / times;
}

int PhyloTree::testAllBranches(int threshold, double best_score, double* pattern_lh, int reps, int lbp_reps,
    PhyloNode* node, PhyloNode* dad) {
    int num_low_support = 0;
    if (!node) {
        node = (PhyloNode*)root;
        root->neighbors[0]->node->name = "";
        if (isSuperTree()) {
            int tmp = save_all_trees;
            save_all_trees = 2;
            bool nni5 = params->nni5;
            params->nni5 = true; // always optimize 5 branches for accurate SH-aLRT
            initPartitionInfo();
            params->nni5 = nni5;
            save_all_trees = tmp;
        }
    }
    if (dad && !node->isLeaf() && !dad->isLeaf()) {
        double lbp_support;
        int support = round(testOneBranch(best_score, pattern_lh, reps, lbp_reps, node, dad, lbp_support) * 100);
        node->name = convertIntToString(support);
        if (lbp_reps)
            node->name += "/" + convertIntToString(round(lbp_support * 100));
        if (support < threshold)
            num_low_support = 1;
        ((PhyloNeighbor*)node->findNeighbor(dad))->partial_pars[0] = support;
        ((PhyloNeighbor*)dad->findNeighbor(node))->partial_pars[0] = support;
    }
    FOR_NEIGHBOR_IT(node, dad, it)num_low_support += testAllBranches(threshold, best_score, pattern_lh, reps, lbp_reps, (PhyloNode*)(*it)->node, node);

    return num_low_support;
}

/****************************************************************************
 Collapse stable (highly supported) clades by one representative
 ****************************************************************************/

void PhyloTree::deleteLeaf(Node* leaf) {

    Node* near_node = leaf->neighbors[0]->node;
    assert(leaf->isLeaf() && near_node->degree() == 3);
    Node* node1 = NULL;
    Node* node2 = NULL;
    double sum_len = 0.0;

    FOR_NEIGHBOR_IT(near_node, leaf, it) {
        sum_len += (*it)->length;
        if (!node1)
            node1 = (*it)->node;

        else
            node2 = (*it)->node;
    }
    // make sure that the returned node1 and node2 are correct
    assert(node1 && node2);
    // update the neighbor
    node1->updateNeighbor(near_node, node2, sum_len);
    node2->updateNeighbor(near_node, node1, sum_len);
}

void PhyloTree::reinsertLeaf(Node* leaf, Node* node, Node* dad) {

    bool first = true;
    Node* adjacent_node = leaf->neighbors[0]->node;
    Neighbor* nei = node->findNeighbor(dad);
    //double len = nei->length;
    double len = max(nei->length, MIN_BRANCH_LEN * 2);
    // to avoid too small branch length when reinserting leaf

    FOR_NEIGHBOR_IT(adjacent_node, leaf, it) {
        if (first) {
            (*it)->node = node;
            (*it)->length = len / 2;
            node->updateNeighbor(dad, adjacent_node, len / 2);
        }
        else {
            (*it)->node = dad;
            (*it)->length = len / 2;
            dad->updateNeighbor(node, adjacent_node, len / 2);
        }
        first = false;
    }
}

bool PhyloTree::isSupportedNode(PhyloNode* node, int min_support) {
    FOR_NEIGHBOR_IT(node, NULL, it)if (!(*it)->node->isLeaf())
        if (((PhyloNeighbor*)*it)->partial_pars[0] < min_support) {

            return false;
        }
    return true;
}

int PhyloTree::collapseStableClade(int min_support, NodeVector& pruned_taxa, StrVector& linked_name,
    double*& dist_mat) {
    NodeVector taxa;
    NodeVector::iterator tax_it;
    StrVector::iterator linked_it;
    getTaxa(taxa);
    IntVector linked_taxid;
    linked_taxid.resize(leafNum, -1);
    int num_pruned_taxa; // global num of pruned taxa
    int ntaxa = leafNum;
    do {
        num_pruned_taxa = 0;
        for (tax_it = taxa.begin(); tax_it != taxa.end(); tax_it++)
            if (linked_taxid[(*tax_it)->id] < 0) {
                Node* taxon = (*tax_it);
                PhyloNode* near_node = (PhyloNode*)taxon->neighbors[0]->node;
                Node* adj_taxon = NULL;
                FOR_NEIGHBOR_DECLARE(near_node, taxon, it)
                    if ((*it)->node->isLeaf()) {
                        adj_taxon = (*it)->node;
                        break;
                    }
                // if it is not a cherry
                if (!adj_taxon)
                    continue;
                assert(linked_taxid[adj_taxon->id] < 0);
                PhyloNeighbor* near_nei = NULL;
                FOR_NEIGHBOR(near_node, taxon, it)
                    if ((*it)->node != adj_taxon) {
                        near_nei = (PhyloNeighbor*)(*it);
                        break;
                    }
                assert(near_nei);
                // continue if the cherry is not stable, or distance between two taxa is near ZERO
                if (!isSupportedNode((PhyloNode*)near_nei->node, min_support)
                    && dist_mat[taxon->id * ntaxa + adj_taxon->id] > 2e-6)
                    continue;
                // now do the taxon pruning
                Node* pruned_taxon = taxon, * stayed_taxon = adj_taxon;
                // prune the taxon that is far away
                if (adj_taxon->neighbors[0]->length > taxon->neighbors[0]->length) {
                    pruned_taxon = adj_taxon;
                    stayed_taxon = taxon;
                }
                deleteLeaf(pruned_taxon);
                linked_taxid[pruned_taxon->id] = stayed_taxon->id;
                pruned_taxa.push_back(pruned_taxon);
                linked_name.push_back(stayed_taxon->name);
                num_pruned_taxa++;
                // do not prune more than n-4 taxa
                if (pruned_taxa.size() >= ntaxa - 4)
                    break;
            }
    } while (num_pruned_taxa && pruned_taxa.size() < ntaxa - 4);

    if (pruned_taxa.empty())
        return 0;

    if (verbose_mode >= VB_MED)
        for (tax_it = pruned_taxa.begin(), linked_it = linked_name.begin(); tax_it != pruned_taxa.end();
            tax_it++, linked_it++)
            cout << "Delete " << (*tax_it)->name << " from " << (*linked_it) << endl;

    // set root to the first taxon which was not deleted
    for (tax_it = taxa.begin(); tax_it != taxa.end(); tax_it++)
        if (linked_taxid[(*tax_it)->id] < 0) {
            root = (*tax_it);
            break;
        }
    // extract the sub alignment
    IntVector stayed_id;
    int i, j;
    for (i = 0; i < taxa.size(); i++)
        if (linked_taxid[i] < 0)
            stayed_id.push_back(i);
    assert(stayed_id.size() + pruned_taxa.size() == leafNum);
    Alignment* pruned_aln = new Alignment();
    pruned_aln->extractSubAlignment(aln, stayed_id, 2); // at least 2 informative characters
    nodeNum = leafNum = stayed_id.size();
    initializeTree();
    setAlignment(pruned_aln);

    double* pruned_dist = new double[leafNum * leafNum];
    for (i = 0; i < leafNum; i++)
        for (j = 0; j < leafNum; j++)
            pruned_dist[i * leafNum + j] = dist_mat[stayed_id[i] * ntaxa + stayed_id[j]];
    dist_mat = pruned_dist;

    return pruned_taxa.size();
}

int PhyloTree::restoreStableClade(Alignment* original_aln, NodeVector& pruned_taxa, StrVector& linked_name) {
    //int num_inserted_taxa;
    NodeVector::reverse_iterator tax_it;
    StrVector::reverse_iterator linked_it;
    tax_it = pruned_taxa.rbegin();
    linked_it = linked_name.rbegin();
    for (; tax_it != pruned_taxa.rend(); tax_it++, linked_it++) {
        //cout << "Reinsert " << (*tax_it)->name << " to " << (*linked_it) << endl;
        Node* linked_taxon = findNodeName((*linked_it));
        assert(linked_taxon);
        assert(linked_taxon->isLeaf());
        leafNum++;
        reinsertLeaf((*tax_it), linked_taxon, linked_taxon->neighbors[0]->node);
    }
    assert(leafNum == original_aln->getNSeq());
    nodeNum = leafNum;
    initializeTree();
    setAlignment(original_aln);
    root = findNodeName(aln->getSeqName(0));
    //if (verbose_mode >= VB_MED) drawTree(cout);

    return 0;
}

bool PhyloTree::checkEqualScalingFactor(double& sum_scaling, PhyloNode* node, PhyloNode* dad) {
    if (!node)
        node = (PhyloNode*)root;
    if (dad) {
        double scaling = ((PhyloNeighbor*)node->findNeighbor(dad))->lh_scale_factor
            + ((PhyloNeighbor*)dad->findNeighbor(node))->lh_scale_factor;
        if (sum_scaling > 0)
            sum_scaling = scaling;
        if (fabs(sum_scaling - scaling) > 1e-6) {
            cout << sum_scaling << " " << scaling << endl;
            return false;
        }
    }
    FOR_NEIGHBOR_IT(node, dad, it)if (!checkEqualScalingFactor(sum_scaling, (PhyloNode*)(*it)->node, node)) return false;

    return true;
}

void PhyloTree::randomizeNeighbors(Node* node, Node* dad) {

    if (!node)
        node = root;
    FOR_NEIGHBOR_IT(node, dad, it)randomizeNeighbors((*it)->node, node);

    my_random_shuffle(node->neighbors.begin(), node->neighbors.end());
}

void PhyloTree::printTransMatrices(Node* node, Node* dad) {
    if (!node)
        node = root;
    int nstates = aln->num_states;

    if (dad) {
        double* trans_cat = new double[nstates * nstates];
        model_factory->computeTransMatrix(dad->findNeighbor(node)->length * site_rate->getRate(0), trans_cat);
        cout << "Transition matrix " << dad->name << " to " << node->name << endl;
        for (int i = 0; i < nstates; i++) {
            for (int j = 0; j < nstates; j++) {
                cout << "\t" << trans_cat[i * nstates + j];
            }
            cout << endl;
        }
        delete[] trans_cat;
    }
    FOR_NEIGHBOR_IT(node, dad, it)printTransMatrices((*it)->node, node);
}

// void PhyloTree::addRemainRow(const vector<string>& remainRowName, const vector<string>& remainRow, const vector<int>& perm, const vector<int>&permCol) {
//     assert(root);
//     Node *new_taxon;
//     int k = perm.size();
//     int curIdx = aln->getNSeq();

//     for(int i = 0; i < k; ++i) {
//         aln->addToAlignmentNewSeq(remainRowName[perm[i]], remainRow[perm[i]], permCol);
//     }

//     for (int i = 0; i < k; ++ i) {
//         initializeAllPartialPars();
//         clearAllPartialLH();
//         // // // allocate a new taxon and a new adjacent internal node
//         new_taxon = newNode(curIdx, remainRowName[perm[i]].c_str());
//         ++curIdx;
//         Node *added_node = newNode();
//         added_node->addNeighbor(new_taxon, -1.0);
//         new_taxon->addNeighbor(added_node, -1.0);
//         ((PhyloNeighbor*) added_node->findNeighbor(new_taxon))->partial_pars = newBitsBlock();
//         ((PhyloNeighbor*) new_taxon->findNeighbor(added_node))->partial_pars = newBitsBlock();
//         // // preserve two neighbors
//         added_node->addNeighbor((Node*) 1, -1.0);
//         added_node->addNeighbor((Node*) 2, -1.0);
//         Node *target_node = NULL;
//         Node *target_dad = NULL;
//         int score = addTaxonMPFast(added_node, target_node, target_dad, root->neighbors[0]->node, root);
//         delete[] ((PhyloNeighbor*) new_taxon->findNeighbor(added_node))->partial_pars;
//         delete[] ((PhyloNeighbor*) added_node->findNeighbor(new_taxon))->partial_pars;
//         if (verbose_mode >= VB_MAX)
//             cout << ", score = " << score << endl;
//         // now insert the new node in the middle of the branch node-dad
//         double len = target_dad->findNeighbor(target_node)->length;
//         target_node->updateNeighbor(target_dad, added_node, -1.0);
//         target_dad->updateNeighbor(target_node, added_node, -1.0);
//         added_node->updateNeighbor((Node*) 1, target_node, -1.0);
//         added_node->updateNeighbor((Node*) 2, target_dad, -1.0);
//         // compute the likelihood
//         // clearAllPartialLh();
//         // optimizeAllBranches();
//         // optimizeNNI();
//     }
// }      

void PhyloTree::computePartialMutation(UINT *states_dad, vector<int> &permCol, PhyloNeighbor *dad_branch, PhyloNode *dad)
{
    // cout << "3\n";
    Node* node = dad_branch->node;
    //assert(node->degree() <= 3);
    int ptn;
    int nstates = aln->num_states;
    int pars_size = getBitsBlockSize();
    int entry_size = getBitsEntrySize();
    int nptn = aln->size();
    int ptn_pars_start_id = pars_size - nptn - 1;

    assert(dad_branch->partial_pars);

    if (nstates == 4 && aln->seq_type == SEQ_DNA && (node->isLeaf() || node->degree() == 3)) {
        // ULTRAFAST VERSION FOR DNA, assuming that UINT is 32-bit integer
        if (node->isLeaf() && dad) {
            // external node
            PhyloNeighbor* node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
            int p = -1;
            for (ptn = 0; ptn < aln->size(); ptn += 8) {
                UINT states = 0;
                int maxi = aln->size() - ptn;
                if (maxi > 8) maxi = 8;
                for (int i = 0; i < maxi; i++) {
                    ++p;
                    UINT node_state = (dad_branch->partial_pars[ptn / 8] >> (i * 4)) & 15;
                    UINT dad_state = (states_dad[ptn/8] >> (i * 4)) & 15;
                    int c = 0;
                    for(c = 0; c < 4; ++c)
                        if(1&(dad_state >> c))
                            break;
                    
                    if((1&(node_state >> c)) == 0)
                    {
                        Mutation mut;
                        int c1 = 0;
                        for(c1 = 0; c1 < 4; ++c1)
                            if(1&(node_state >> c1))
                                break;
                        mut.position = permCol[p];
                        mut.mut_nuc = c1;
                        mut.par_nuc = c;
                        node_branch->mutations.push_back(mut);
                    }
                }
            }
        }
        else {
            // internal node
            UINT* left = NULL, * right = NULL;
            PhyloNeighbor *left_branch, *right_branch;
            int pars_steps = 0;
            UINT *save_fitch_result = new UINT[(aln->size()) + 1];
            FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME) {
                if (!left)
                    left = ((PhyloNeighbor*)(*it))->partial_pars, left_branch = (PhyloNeighbor*)(*it)->node->findNeighbor(node);
                else
                    right = ((PhyloNeighbor*)(*it))->partial_pars, right_branch = (PhyloNeighbor*)(*it)->node->findNeighbor(node);
            }
            int p = -1;
            for (ptn = 0; ptn < aln->size(); ptn += 8) {
                // cout << dad_branch->partial_pars[pars_size - 1] << ": ***\n";
                UINT left_state = left[ptn / 8];
                UINT right_state = right[ptn / 8];
                UINT dad_state = states_dad[ptn / 8];
                int maxi = aln->size() - ptn;
                if (maxi > 8) maxi = 8;
                for (int i = 0; i < maxi; i++) {
                    ++p;
                    UINT state_left = (left_state >> (i * 4)) & 15;
                    UINT state_right = (right_state >> (i * 4)) & 15;
                    UINT tmp = state_left | (state_right << 4);
                    save_fitch_result[p] = dna_fitch_step[tmp];
                    UINT state_both = (dad_state >> (i * 4)) & 15;
                    // cout << state_left << " " << state_right << " " << state_both << ": ";
                    int c = 0;
                    for(c = 0; c < 4; ++c)
                        if(1&(state_both >> c))
                            break;
                    
                    if((1&(state_left >> c)) == 0 && left_branch != NULL)
                    {
                        Mutation mut_l;
                        int c1 = 0;
                        for(c1 = 0; c1 < 4; ++c1)
                            if(1&(state_left >> c1))
                                break;
                        mut_l.position = permCol[p];
                        mut_l.mut_nuc = c1;
                        mut_l.par_nuc = c;
                        left_branch->mutations.push_back(mut_l);
                        // cout << "1+" << mut_l.position << " ";
                    }

                    if((1&(state_right>>c)) == 0 && right_branch != NULL)
                    {
                        Mutation mut_r;
                        int c1 = 0;
                        for(c1 = 0; c1 < 4; ++c1)
                            if(1&(state_right >> c1))
                                break;
                        mut_r.position = permCol[p];
                        mut_r.mut_nuc = c1;
                        mut_r.par_nuc = c;
                        right_branch->mutations.push_back(mut_r);
                        // cout << "2+" << mut_r.position << " ";
                    }
                    // cout << '\n';
                }
            }
            FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME) {
                UINT *new_states_dad = new UINT[(aln->size() + 7) / 8 + 1];
                for(int ptn = 0; ptn < aln->size(); ptn += 8)
                {
                    new_states_dad[ptn/8] = (states_dad[ptn/8] & (((PhyloNeighbor*)(*it))->partial_pars[ptn/8]));
                    int maxi = aln->size() - ptn;
                    if (maxi > 8) maxi = 8;
                    for(int i = 0; i < maxi; ++i)
                    {
                        if(save_fitch_result[ptn+i] == 1)
                        {
                            UINT cur_state = ((((PhyloNeighbor*)(*it))->partial_pars[ptn/8]) >> (i*4)) & 15;
                            new_states_dad[ptn/8] |= (cur_state << (i*4));
                        }
                    }
                }
                computePartialMutation(new_states_dad, permCol, (PhyloNeighbor*)(*it), (PhyloNode*)node);
            }
        }
        return;
    } // END OF DNA VERSION
}

void PhyloTree::computeMutationBranch(vector<int> &permCol, PhyloNeighbor* dad_branch, PhyloNode* dad, int* branch_subst) {
    // cout << "2\n";
    PhyloNode* node = (PhyloNode*)dad_branch->node;
    PhyloNeighbor* node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
    assert(node_branch);
    if (!central_partial_pars)
        initializeAllPartialPars();
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode* tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor* tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }

    int nptn = aln->size();
    if (!_pattern_pars) _pattern_pars = aligned_alloc<BootValTypePars>(nptn + VCSIZE_USHORT);
    memset(_pattern_pars, 0, sizeof(BootValTypePars) * (nptn + VCSIZE_USHORT));

    UINT *lf = new UINT[(aln->size() + 7) / 8 + 1];
    for(int ptn = 0; ptn < aln->size(); ptn += 8)
    {
        lf[ptn/8] = (save_branch_states_dad[ptn/8] & dad_branch->partial_pars[ptn/8]);
        int maxi = aln->size() - ptn;
        if (maxi > 8) maxi = 8;
        for(int i = 0; i < maxi; ++i)
        {
            if(save_branch_fitch_result[ptn+i] == 1)
            {
                UINT cur_state = (dad_branch->partial_pars[ptn/8] >> (i*4)) & 15;
                lf[ptn/8] |= (cur_state << (i*4));
            }
        }
    }
    computePartialMutation(lf, permCol, dad_branch, dad);
    UINT *rg = new UINT[(aln->size() + 7) / 8 + 1];
    for(int ptn = 0; ptn < aln->size(); ptn += 8)
    {
        rg[ptn/8] = (save_branch_states_dad[ptn/8] & node_branch->partial_pars[ptn/8]);
        int maxi = aln->size() - ptn;
        if (maxi > 8) maxi = 8;
        for(int i = 0; i < maxi; ++i)
        {
            if(save_branch_fitch_result[ptn+i] == 1)
            {
                UINT cur_state = (node_branch->partial_pars[ptn/8] >> (i*4)) & 15;
                rg[ptn/8] |= (cur_state << (i*4));
            }
        }
    }
    computePartialMutation(rg, permCol, node_branch, node);

    int i, ptn, p = -1, tree_pars = 0;

    if (aln->num_states == 4 && aln->seq_type == SEQ_DNA) {
        // ULTRAFAST VERSION FOR DNA
        for (ptn = 0; ptn < aln->size(); ptn += 8) {
            UINT states_left = node_branch->partial_pars[ptn / 8];
            UINT states_right = dad_branch->partial_pars[ptn / 8];
            UINT states_dad = save_branch_states_dad[ptn/8];
            int maxi = aln->size() - ptn;
            if (maxi > 8) maxi = 8;
            for (i = 0; i < maxi; i++) {
                ++p;
                UINT state_left = (states_left >> (i * 4)) & 15;
                UINT state_right = (states_right >> (i * 4)) & 15;
                UINT state_both = (states_dad >> (i * 4)) & 15;
                int c = 0;
                for(c = 0; c < 4; ++c)
                    if(1&(state_both >> c))
                        break;
                
                if(((1&(state_left>>c))) == 0)
                {
                    Mutation mut_l;
                    int c1 = 0;
                    for(c1 = 0; c1 < 4; ++c1)
                        if(1&(state_left >> c1))
                            break;
                    mut_l.position = permCol[p];
                    mut_l.mut_nuc = c1;
                    mut_l.par_nuc = c;
                    dad_branch->mutations.push_back(mut_l);
                    ++tree_pars;
                }

                if((1&(state_right>>c)) == 0)
                {
                    Mutation mut_r;
                    int c1 = 0;
                    for(c1 = 0; c1 < 4; ++c1)
                        if(1&(state_right >> c1))
                            break;
                    mut_r.position = permCol[p];
                    mut_r.mut_nuc = c1;
                    mut_r.par_nuc = c;
                    node_branch->mutations.push_back(mut_r);
                    ++tree_pars;
                }
            }
        }
    }

    // cout << "Mutation: " << tree_pars << '\n';
}

void PhyloTree::initMutation(vector<int> &permCol)
{
    // cout << "1\n";
    assert(root->isLeaf());
    PhyloNeighbor* nei = ((PhyloNeighbor*)root->neighbors[0]);
    current_it = nei;
    assert(current_it);
    current_it_back = (PhyloNeighbor*)nei->node->findNeighbor(root);
    assert(current_it_back);

    int nptn = aln->size();
    if (_pattern_pars == NULL) _pattern_pars = aligned_alloc<BootValTypePars>(nptn + VCSIZE_USHORT);

    computeMutationBranch(permCol, (PhyloNeighbor*)root->neighbors[0], (PhyloNode*)root);
}

int PhyloTree::computePartialParsimonyMutation(PhyloNeighbor *dad_branch, PhyloNode *dad)
{
    int par_s = 0;
    PhyloNode* node = (PhyloNode*)dad_branch->node;
    PhyloNeighbor* node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
    par_s += node_branch->mutations.size();
    FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME)
    {
        par_s += computePartialParsimonyMutation(((PhyloNeighbor*)(*it)), node);
    }
    return par_s;
}

int PhyloTree::computeParsimonyBranchMutation(PhyloNeighbor* dad_branch, PhyloNode* dad, int* branch_subst)
{
    PhyloNode* node = (PhyloNode*)dad_branch->node;
    PhyloNeighbor* node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
    assert(node_branch);
    if (!central_partial_pars)
        initializeAllPartialPars();
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode* tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor* tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }

    int par_s = 0;
    par_s += computePartialParsimonyMutation(dad_branch, dad);
    par_s += computePartialParsimonyMutation(node_branch, node);
    return par_s;
}


int PhyloTree::computeParsimonyScoreMutation()
{
    assert(root->isLeaf());
    PhyloNeighbor* nei = ((PhyloNeighbor*)root->neighbors[0]);
    current_it = nei;
    assert(current_it);
    current_it_back = (PhyloNeighbor*)nei->node->findNeighbor(root);
    assert(current_it_back);

    int parsimonyScore = 0;
    parsimonyScore += computeParsimonyBranchMutation((PhyloNeighbor*)root->neighbors[0], (PhyloNode*)root);
    return parsimonyScore;
}


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 