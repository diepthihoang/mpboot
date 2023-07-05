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
// #include "phylonode.h"
// Forward declaration of structs from usher_graph

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
        std::string name;
        std::vector<Mutation> mutations;

        MutationNode();

        void add_mutation(Mutation mut);
        void clear_mutations();
};
#endif