/*
 * sprparsimony.h
 *
 *  Created on: Nov 6, 2014
 *      Author: diep
 */

#ifndef SPRPARSIMONY_H_
#define SPRPARSIMONY_H_

#include "iqtree.h"

/*
 * An alternative for pllComputeRandomizedStepwiseAdditionParsimonyTree
 * because the original one seems to have the wrong deallocation function
 */
void _pllComputeRandomizedStepwiseAdditionParsimonyTree(pllInstance * tr, partitionList * partitions, int sprDist);

void _allocateParsimonyDataStructures(pllInstance *tr, partitionList *pr);
void _pllFreeParsimonyDataStructures(pllInstance *tr, partitionList *pr);

/**
 * DTH: optimize whatever tree is stored in tr by parsimony SPR
 * @param tr: the tree instance :)
 * @param partition: the data partition :)
 * @param mintrav, maxtrav are PLL limitations for SPR radius
 * @return best parsimony score found
 */
int pllOptimizeSprParsimony(pllInstance * tr, partitionList * pr, int mintrav, int maxtrav, IQTree *iqtree);

int pllSaveCurrentTreeSprParsimony(pllInstance * tr, partitionList * pr, int cur_search_pars);

void pllComputePatternParsimony(pllInstance * tr, partitionList * pr, double *ptn_npars, double *cur_npars);
void pllComputeSiteParsimony(pllInstance * tr, partitionList * pr, int *site_pars, int *cur_pars = NULL);

// Diep: for testing site parsimony computed by PLL vs IQTree on the same tree
// this is called if params.test_site_pars == true
void testSiteParsimony(Params &params);


#endif /* SPRPARSIMONY_H_ */
