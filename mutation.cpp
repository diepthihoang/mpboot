#include "mutation.h"
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <random>
#include <algorithm>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
// Uses one-hot encoding if base is unambiguous
// A:1,C:2,G:4,T:8
MutationNode::MutationNode()
{
    mutations.clear();
}

// Assumes mutations are added in chronological order. If a new mutation occurs
// at the same position, it should either be updated to the new allele or
// removed entirely (in case of reversal mutation)
void MutationNode::add_mutation(Mutation mut)
{
    auto iter = std::lower_bound(mutations.begin(), mutations.end(), mut);
    // check if mutation at the same position has occured before
    if ((iter != mutations.end()) && (iter->position == mut.position))
    {
        // update to new allele
        if (iter->par_nuc != mut.mut_nuc)
        {
            iter->mut_nuc = mut.mut_nuc;
        }
        // reversal mutation
        else
        {
            if (iter->mut_nuc != mut.par_nuc)
            {
                printf("ERROR: add_mutation: consecutive mutations at same position "
                                "disagree on nuc -- called out of order?\n");
                exit(1);
            }
            std::vector<Mutation> tmp;
            for (auto m : mutations)
            {
                if (m.position != iter->position)
                {
                    tmp.emplace_back(m.copy());
                }
            }
            mutations.clear();
            for (auto m : tmp)
            {
                mutations.emplace_back(m.copy());
            }
        }
    }
    // new mutation
    else
    {
        mutations.insert(iter, mut);
    }
}

void MutationNode::clear_mutations()
{
    mutations.clear();
}