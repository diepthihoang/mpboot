#ifndef _MUTATION
#define _MUTATION 
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <cassert>
#include "phylonode.h"
// Forward declaration of structs from usher_graph
struct Missing_Sample;

#if SAVE_PROFILE == 1
#define TIMEIT() InstrumentationTimer timer##__LINE__(__PRETTY_FUNCTION__);
#else
#define TIMEIT()
#endif

struct Mutation
{
    std::string name, ref_name;
    int position;
    int ref_nuc;
    int par_nuc;
    int mut_nuc;
    bool is_missing;
    inline bool operator<(const Mutation &m) const
    {
        return ((*this).position < m.position);
    }
    inline Mutation copy() const
    {
        Mutation m;
        m.name = name;
        m.position = position;
        m.ref_nuc = ref_nuc;
        m.par_nuc = par_nuc;
        m.mut_nuc = mut_nuc;
        m.is_missing = is_missing;
        return m;
    }
    Mutation()
    {
        name = "";
        is_missing = false;
    }
    inline bool is_masked() const
    {
        return (position < 0);
    }
};

class MutationNode
{
public:
    size_t level;
    float branch_length;
    std::string identifier;
    std::vector<std::string> clade_annotations;
    MutationNode *parent;
    std::vector<MutationNode *> children;
    std::vector<Mutation> mutations;
    size_t dfs_idx;
    size_t dfs_end_idx;

    bool is_leaf();
    bool is_root();

    MutationNode();
    MutationNode(std::string id, float l);
    MutationNode(std::string id, MutationNode *p, float l);

    void add_mutation(Mutation mut);
    void clear_mutations();
    void clear_annotations();
};

class CandidateNode
{
    public:
        std::string missing_sample;
        PhyloNode* node;
        PhyloNeighbor *node_branch;
        std::vector<Mutation>* missing_sample_mutations;

        int* best_set_difference;
        int* set_difference;
        size_t* best_node_num_leaves;
        size_t j;
        size_t* best_j;
        size_t distance;
        size_t* best_distance;
        size_t* num_best;
        PhyloNode** best_node;
        PhyloNeighbor** best_node_branch;

        std::vector<bool>* node_has_unique;
        std::vector<size_t>* best_j_vec;

        bool* has_unique;

        std::vector<Mutation>* excess_mutations;
        std::vector<Mutation>* imputed_mutations;

        CandidateNode () {
            distance = 0;
            best_distance = &distance;
        }
};
#endif