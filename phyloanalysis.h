/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef PHYLOANALYSIS_H
#define PHYLOANALYSIS_H

#include "tools.h"
#include "mexttree.h"
#include "phylotesting.h"
#include "nnisearch.h"

class PhyloTree;
class IQTree;

/**
	main function to carry out phylogenetic inference
	@param params program parameters
*/
void runPhyloAnalysis(Params &params);

void runTreeReconstruction(Params &params, string &original_model,
		IQTree &tree, vector<ModelInfo> &model_info);

/**
	take the collection of trees from input_trees, it assign support values to target_tree
	and print resulting tree to output_tree. 
	@param input_trees collection of input trees to infer split supports
	@param burnin the number trees at the beginning of input_trees to be discarded
	@param max_count max number of trees to load
	@param target_tree tree to assign support value
	@param output_tree (OUT, OVERWRITE IF EXIST) Resulting will be written to this file. If NULL,
		output_tree will be named target_tree appended with ".suptree"
*/
void assignBootstrapSupport(const char *input_trees, int burnin, int max_count, const char *target_tree, 
	bool rooted, const char *output_tree, const char *out_prefix, MExtTree &mytree, 
	const char* tree_weight_file, Params *params);

/**
 * assign branch supports from params.user_tree trees file to params.second_tree
 * @param params program parameters
 */
void assignBranchSupportNew(Params &params);

/**
	Compute the consensus tree from the collection of trees from input_trees
	and print resulting tree to output_tree. 
	@param phylo_tree used to optimize branch lengths of the consensus tree. Can be NULL
	@param input_trees collection of input trees to infer split supports
	@param burnin the number trees at the beginning of input_trees to be discarded
	@param max_count max number of trees to load
	@param cutoff only incorporate those splits that have support values more than cutoff
	@param weight_threshold minimum weight cutoff
	@param output_tree (OUT, OVERWRITE IF EXIST) Resulting consensus tree will be written to this file. If NULL,
		output_tree will be named input_trees appended with ".contree"
*/
void computeConsensusTree(const char *input_trees, int burnin, int max_count, double cutoff, double weight_threshold,
	const char *output_tree, const char *out_prefix, const char* tree_weight_file, Params *params);

/**
	Compute the consensus network from the collection of trees in input_trees.
	print consensus network to output_tree
	@param input_trees collection of input trees to infer split supports
	@param burnin the number trees at the beginning of input_trees to be discarded
	@param max_count max number of trees to load
	@param cutoff only incorporate those splits that have support values more than cutoff
	@param weight_threshold minimum weight cutoff
	@param output_tree (OUT, OVERWRITE IF EXIST) Resulting consensus tree will be written to this file. If NULL,
		output_tree will be named input_trees appended with ".connetwork"
*/
void computeConsensusNetwork(const char *input_trees, int burnin, int max_count, double cutoff,
		int weight_summary, double weight_threshold,
	const char *output_tree, const char *out_prefix, const char* tree_weight_file);

/**
 * Diep:
 * Locate the definition here so it is seen by optimizeAlignment
 */
//#define BootValTypePars int // Diep added
#define BootValTypePars unsigned short // Diep added


struct PatternComp{
	bool operator() (Pattern i, Pattern j) {
		return (i.ras_pars_score * i.frequency > j.ras_pars_score * j.frequency);
//		return (i.ras_pars_score > j.ras_pars_score);
	}
};

/**
 * Diep:
 * To print site parsimony scores for a user input tree
 */
void printSiteParsimonyUserTree(Params &params);

/**
 * Diep:
 * Based on a topology output by (PLL Random Addition Sequence method) IQTREE RAS / user tree
 * reorder pattern in the alignment by decreasing order of parsimony score
 * By this, uninformative sites will all be shifted to the end of the alignment.
 */
void optimizeAlignment(IQTree * & tree, Params & params);

void testCompConsensus(const char * infile, const char * outfile, Params *params);

string computeConsensusTreeNoFileIO(StringIntMap& input_trees, IntVector & weight, int max_count,
		double cutoff, double weight_threshold, Params *params);

#endif
