/*
 * candidateset.h
 *
 *  Created on: Jun 1, 2014
 *      Author: Tung Nguyen
 */

#ifndef CANDIDATESET_H_
#define CANDIDATESET_H_
#include "tools.h"
#include "alignment.h"
#include <stack>

struct CandidateTree {
	string tree; // with branch length
	string topology; // tree topology WITHOUT branch lengths and WITH TAXON ID (instead of taxon names) for sorting purpose
	double score; // log-likelihood under ML or parsimony score
};


/**
 * Candidate tree set, sorted in ascending order of scores, i.e. the last element is the highest scoring tree
 */
class CandidateSet : public multimap<double, CandidateTree> {

public:
    /**
     * constructor
     */
	CandidateSet(int limit, int max_candidates, Alignment *aln);

	CandidateSet();

    /**
     * return randomly one candidate tree from max_candidate
     */
    string getRandCandTree();

    /**
     * return the next parent tree for reproduction.
     * Here we always maintain a list of candidate trees which have not
     * been used for reproduction. If all candidate trees have been used, we select the
     * current best trees as the new parent trees
     */
    string getNextCandTree();

    /**
     *  Replace an existing tree in the candidate set
     *  @param tree the new tree string that will replace the existing tree
     *  @param score the score of the new tree
     *  @return true if the topology of \a tree exist in the candidate set
     */
    bool replaceTree(string tree, double score);

    /**
     *  create the parent tree set containing top trees
     */
    void initParentTrees();

    /**
     * update / insert tree into set of score is higher than lowest-scoring tree
     * @return false if the tree topology already exists
     *
     */
    bool update(string tree, double score);

    /**
     *  print score of max_candidates best trees
     *
     *  @param numScore
     *  	Number of best scores to print out starting from the highest
     */
    vector<double> getBestScores(int numBestScores = 0);

    /**
     * Return the worst score in the candidate tree set
     * @return
     */
    double getWorstScore();

    /**
     *  Return \a numTree best tree strings
     *  @param numTree number of best trees
     *  @return a list of tree
     */
    vector<string> getHighestScoringTrees(int numTree = 0);

    /**
     * get tree(s) with highest score. More than one tree is
     * returned if there are multiple optima.
     * @return a vector containing optimal trees
     */
    vector<string> getEquallyOptimalTrees();

    /**
     * destructor
     */
    virtual ~CandidateSet();

    /**
     * hard limit for number of trees (typically superset of candidate set)
     */
    int maxCandidates;

    /**
     *  best score in the set
     */
    double bestScore;

    /**
     *  maximum number of candidate trees
     */
    int popSize;

    /** index of tree topologies in set
     *
     */
    StringDoubleHashMap topologies;

    /**
     *  Trees used for reproduction
     */
    stack<string> parentTrees;

    /**
     * pointer to alignment, just to assign correct IDs for taxa
     */
    Alignment *aln;

    /**
     * check if tree topology WITHOUT branch length exist in the candidate set?
     */
    bool treeTopologyExist(string topo);

    /**
     * check if tree topology WITH branch length exist in the candidate set?
     */
    bool treeExist(string tree);

    /**
     * return a unique topology (sorted by taxon names, rooted at taxon with alphabetically smallest name) without branch lengths
     */
    string getTopology(string tree);

    /**
     *  Empty the candidate set
     */
    void clear();

    /**
     * Return a CandidateSet containing \a numTrees of current best candidate trees
     * @param numTrees
     * @return
     */
    CandidateSet getBestCandidateTrees(int numTrees);

	/*
	 * Diep added
	 * Copy all tree strings in current candidate set into vector<string> candidateTreeVec
	 */
    void copyToCandidateVec(); // Diep added to avoid bias in support values for big group

    /*
     * Diep added
     * Return a random tree stored in vector<string> candidateTreeVec
     * Somehow this works for MP better than getRandCandTree
     */
    string getRandCandVecTree();

    void syncCandidateSet();

    string getSyncTrees(int forWorker=5);

    string getCandidateSetAsSyncTrees(int limit);

    void updateSingleSyncTree(string singleTree);

    void updateSyncTrees(string syncString);
private:
	vector<string> candidateTreeVec; // Diep added to avoid bias in support values for big group
};

#endif /* CANDIDATESET_H_ */
