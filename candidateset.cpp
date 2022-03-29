/*
 * candidateset.cpp
 *
 *  Created on: Jun 1, 2014
 *      Author: Tung Nguyen
 */

#include "phylotree.h"
#include "candidateset.h"

CandidateSet::CandidateSet(int limit, int max_candidates, Alignment *aln) {
    assert(max_candidates <= limit);
    assert(aln);
    this->maxCandidates = limit;
    this->popSize = max_candidates;
    this->aln = aln;
    this->bestScore = -DBL_MAX;
}

CandidateSet::CandidateSet() {
	aln = NULL;
	maxCandidates = 0;
	popSize = 0;
	bestScore = -DBL_MAX;
}

vector<string> CandidateSet::getEquallyOptimalTrees() {
	vector<string> res;
	for (reverse_iterator rit = rbegin(); rit != rend() && rit->second.score == bestScore; rit++) {
		res.push_back(rit->second.tree);
	}
	return res;
}

string CandidateSet::getRandCandTree() {
	assert(!empty());
	if (empty())
		return "";
	int id = random_int(min(popSize, (int)size()) );
	for (reverse_iterator i = rbegin(); i != rend(); i++, id--)
		if (id == 0)
			return i->second.tree;
	assert(0);
	return "";
}

vector<string> CandidateSet::getHighestScoringTrees(int numTree) {
	assert(numTree <= maxCandidates);
	if (numTree == 0) {
		numTree = maxCandidates;
	}
	vector<string> res;
	int cnt = numTree;
	for (reverse_iterator rit = rbegin(); rit != rend() && cnt > 0; rit++, cnt--) {
		res.push_back(rit->second.tree);
	}
	return res;
}

bool CandidateSet::replaceTree(string tree, double score) {
    CandidateTree candidate;
    candidate.tree = tree;
    candidate.score = score;
    candidate.topology = getTopology(tree);
    if (treeTopologyExist(candidate.topology)) {
        topologies[candidate.topology] = score;
        for (reverse_iterator i = rbegin(); i != rend(); i++) {
            if (i->second.topology == candidate.topology) {
                erase( --(i.base()) );
                break;
            }
            insert(CandidateSet::value_type(score, candidate));
        }
    } else {
        return false;
    }
    return true;
}

string CandidateSet::getNextCandTree() {
    string tree;
    assert(!empty());
    if (parentTrees.empty()) {
        initParentTrees();
    }
    tree = parentTrees.top();
    parentTrees.pop();
    return tree;
}

void CandidateSet::initParentTrees() {
    if (parentTrees.empty()) {
        int count = this->popSize;
        for (reverse_iterator i = rbegin(); i != rend() && count >0 ; i++, count--) {
            parentTrees.push(i->second.tree);
            //cout << i->first << endl;
        }
    }
}

bool CandidateSet::update(string tree, double score) {
	bool newTree;
	CandidateTree candidate;
	candidate.tree = tree;
	candidate.score = score;
	candidate.topology = getTopology(tree);
	if (candidate.score > bestScore)
		bestScore = candidate.score;
	if (treeTopologyExist(candidate.topology)) {
	    // if tree topology already exist, we replace the old
	    // by the new one (with new branch lengths) and update the score
		if (topologies[candidate.topology] < score) {
			topologies[candidate.topology] = score;
			for (CandidateSet::iterator i = begin(); i != end(); i++)
				if (i->second.topology == candidate.topology) {
					erase(i);
					break;
				}
			// insert tree into candidate set
			insert(CandidateSet::value_type(score, candidate));
		}
		newTree = false;
	} else {
		newTree = true;
		if (size() < maxCandidates) {
			// insert tree into candidate set
			insert(CandidateSet::value_type(score, candidate));
			topologies[candidate.topology] = score;
			candidateTreeVec.push_back(tree); // Diep added
		} else if (getWorstScore() <= score){
			// remove the worst-scoring tree
			topologies.erase(begin()->second.topology);

			// Diep added
			int n = candidateTreeVec.size();
			for(int i = 0; i < n; i++)
				if(candidateTreeVec[i] == begin()->second.tree){
					candidateTreeVec[i] = tree;
					break;
				}

			erase(begin());
			// insert tree into candidate set
			insert(CandidateSet::value_type(score, candidate));
			topologies[candidate.topology] = score;
		}else
			newTree = false; // Diep added
	}
	return newTree;
}

vector<double> CandidateSet::getBestScores(int numBestScore) {
	if (numBestScore == 0 || numBestScore > size())
		numBestScore = size();
	vector<double> res;
	for (reverse_iterator rit = rbegin(); rit != rend() && numBestScore > 0; rit++, numBestScore--) {
		res.push_back(rit->first);
	}
	return res;
}

double CandidateSet::getWorstScore() {
	return begin()->first;
}

string CandidateSet::getTopology(string tree) {
	PhyloTree mtree;
	mtree.rooted = false;
	mtree.aln = aln;
	mtree.readTreeString(tree);
    mtree.root = mtree.findNodeName(aln->getSeqName(0));
	ostringstream ostr;
	mtree.printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
	return ostr.str();
}

void CandidateSet::clear() {
	multimap<double, CandidateTree>::clear();
	topologies.clear();
	bestScore = -DBL_MAX;
}

CandidateSet CandidateSet::getBestCandidateTrees(int numTrees) {
	CandidateSet res;
	if (numTrees >= size())
		numTrees = size();
	for (reverse_iterator rit = rbegin(); rit != rend() && numTrees > 0; rit++, numTrees--) {
		res.insert(*rit);
	}
	return res;
}

CandidateSet::~CandidateSet() {
}

bool CandidateSet::treeTopologyExist(string topo) {
	return topologies.find(topo) != topologies.end();
}

bool CandidateSet::treeExist(string tree) {
	return treeTopologyExist(getTopology(tree));
}

void CandidateSet::copyToCandidateVec(){
	assert(!empty());
	if (empty())
		return;
	candidateTreeVec.clear();
	for (reverse_iterator i = rbegin(); i != rend(); i++)
		candidateTreeVec.push_back(i->second.tree);
}

string CandidateSet::getRandCandVecTree(){
	assert(!empty());
	if (empty())
		return "";
	int n = candidateTreeVec.size();
	int i = random_int(n);
	return candidateTreeVec[i];
	assert(0);
	return "";
}

string CandidateSet::getCandidateSetAsSyncTrees(int limit) {
	string syncTrees;
	int sz = min((int) size(), limit);

	auto it = rbegin();
	for(int i = 0; i < sz; ++it, ++i) {
		syncTrees += to_string((int) -(it->first));
		syncTrees += ' ';
		syncTrees += (it->second).tree;
		syncTrees += '#';
	}
	return syncTrees;
}

string CandidateSet::getSyncTrees(int forWorker) {
	string syncTrees;
	int sz = min((int) size(), popSize);

	if (MPIHelper::getInstance().isWorker()) {
		sz = min(sz, forWorker);
	}

	auto it = rbegin();
	for(int i = 0; i < sz; ++it, ++i) {
		syncTrees += to_string((int) -(it->first));
		syncTrees += ' ';
		syncTrees += (it->second).tree;
		syncTrees += '#';
	}
	return syncTrees;
}

void CandidateSet::updateSingleSyncTree(string singleTree) {
	double score = 0;
	string sscore;
	string tree;
	for(int i = 0; i < singleTree.size(); ++i) {
		if (singleTree[i] == ' ') {
			tree = singleTree.substr(i + 1, singleTree.size() - 1 - i);
			break;
		} else {
			score = score * 10 + singleTree[i] - '0';
		}
	}
	score = -score;
	update(tree, score);
}


void CandidateSet::updateSyncTrees(string syncString) {
	int startPos = 0;
	for(int i = 0; i < syncString.size(); ++i) {
		if (syncString[i] == '#') {
			string curTree = syncString.substr(startPos, i - startPos);
			updateSingleSyncTree(curTree);
			startPos = i + 1;
		}
	}
}