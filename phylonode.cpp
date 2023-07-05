//
// C++ Implementation: phylonode
//
// Description: 
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "phylonode.h"


void PhyloNeighbor::clearForwardPartialLh(Node *dad) {
	clearPartialLh();
	for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it ++)
		if ((*it)->node != dad)
			((PhyloNeighbor*)*it)->clearForwardPartialLh(node);
}

void PhyloNeighbor::clear_mutations()
{
	mutations.clear();
}

void PhyloNeighbor::add_mutation(Mutation mut)
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

void PhyloNode::clearReversePartialLh(PhyloNode *dad) {
	PhyloNeighbor *node_nei = (PhyloNeighbor*)findNeighbor(dad);
	assert(node_nei);
	node_nei->partial_lh_computed = 0;
	for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it ++)
		if ((*it)->node != dad)
			((PhyloNode*)(*it)->node)->clearReversePartialLh(this);
}

void PhyloNode::clearAllPartialLh(PhyloNode *dad) {
	PhyloNeighbor *node_nei = (PhyloNeighbor*)findNeighbor(dad);
	node_nei->partial_lh_computed = 0;
	node_nei = (PhyloNeighbor*)dad->findNeighbor(this);
	node_nei->partial_lh_computed = 0;
	for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it ++)
		if ((*it)->node != dad)
			((PhyloNode*)(*it)->node)->clearAllPartialLh(this);
}


PhyloNode::PhyloNode()
 : Node()
{
	init();
}


PhyloNode::PhyloNode(int aid) : Node(aid)
{
	init();
}

PhyloNode::PhyloNode(int aid, int aname) : Node (aid, aname) {
	init();
}


PhyloNode::PhyloNode(int aid, const char *aname) : Node(aid, aname) {
	init();
}

void PhyloNode::init() {
	//partial_lh = NULL;
}


void PhyloNode::addNeighbor(Node *node, double length, int id) {
	neighbors.push_back(new PhyloNeighbor(node, length, id));
}
