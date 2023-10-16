#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iqtree_config.h>
#include "phylotree.h"
#include "phylosupertree.h"
#include "phylosupertreeplen.h"
#include "phyloanalysis.h"
#include "alignment.h"
#include "superalignment.h"
#include "iqtree.h"
#include "model/modelgtr.h"
#include "model/modeldna.h"
#include "myreader.h"
#include "model/rateheterogeneity.h"
#include "model/rategamma.h"
#include "model/rateinvar.h"
#include "model/rategammainvar.h"
// #include "modeltest_wrapper.h"
#include "model/modelprotein.h"
#include "model/modelbin.h"
#include "model/modelcodon.h"
#include "stoprule.h"

#include "mtreeset.h"
#include "mexttree.h"
#include "model/ratemeyerhaeseler.h"
#include "whtest_wrapper.h"
#include "model/partitionmodel.h"
#include "guidedbootstrap.h"
#include "model/modelset.h"
#include "timeutil.h"
#include "parstree.h"
#include "tinatree.h"
#include "sprparsimony.h"
#include "placement.h"
#include <algorithm>

int addMoreRowIQTree(IQTree* tree, Alignment* alignment)
{
	auto startTime = getCPUTime();
	// add more K row to this tree
	int k = alignment->remainSeq.size();
	IQTree resTree;
	vector<int> perm;
	for (int i = 0; i < k; ++i)
		perm.push_back(i);
	vector<int> permCol = alignment->findPermCol();
	int bestScore = (int)1e9 + 7;
	int timer = 0;
	do
	{
		// my_random_shuffle(perm.begin(), perm.end());
		++timer;
		if (timer > 1)
			break;
		IQTree newTree;
		char *file_name = "tree.treefile";
		tree->printTree(file_name, WT_SORT_TAXA | WT_NEWLINE);
		bool is_rooted = false;
		(newTree).readTree(file_name, is_rooted);
		// newTree.copyPhyloTree(tree);
		newTree.setAlignment(tree->aln);
		newTree.aln = new Alignment;
		newTree.aln->copyAlignment(tree->aln);

		cout << "tree parsimony score before add more k rows: " << tree->computeParsimonyScore() << " " << newTree.computeParsimonyScore() << '\n';
		newTree.addRemainRow(alignment->remainName, alignment->remainSeq, perm, permCol);
		newTree.initializeAllPartialPars();
		newTree.clearAllPartialLH();
		cout << newTree.computeParsimonyScore() << " " << newTree.computeParsimony() << '\n';
		int curScore = newTree.computeParsimonyScore();
		if (curScore < bestScore)
		{
			bestScore = curScore;
		}
		delete newTree.aln;
		newTree.aln = NULL;
	} while (next_permutation(perm.begin(), perm.end()));
	cout << "\nBest tree parsimony found after add more k rows: ";
	cout << bestScore << '\n';
	cout << "Time: " << fixed << setprecision(3) << getCPUTime() - startTime << " seconds\n";
	return bestScore;
}

const double e = 2.7182818;
int addMoreRowPLL(IQTree* tree, Alignment* alignment, Params& params)
{
	auto startTime = getCPUTime();
	int k = alignment->remainSeq.size();
	IQTree resTree;
	vector<int> permCol = alignment->findPermCol();
	int bestScore;

	if (k <= 5)
	{
		int bestScore = (int)1e9 + 7;
		vector<int> bestPerm;

		vector<int> perm;
		for (int i = 0; i < k; ++i)
			perm.push_back(i);
		do
		{
			int newScore = computeParsimonyPermutation(tree, alignment, params, permCol, perm);
			if (bestScore > newScore)
				bestPerm = perm;
			bestScore = min(bestScore, newScore);
		} while (next_permutation(perm.begin(), perm.end()));
		cout << '\n';
		cout << "Best tree parsimony found after add more k rows using PLL core: ";
		cout << bestScore << '\n';
		cout << "Time: " << getCPUTime() - startTime << " seconds\n";
		return bestScore;
	}

	int numCandidate = 5;
	vector<vector<int> > candidates;
	vector<int> candidateScore;

	// init 100 permutation
	for (int id = 1; id <= 100; ++id)
	{
		vector<int> perm;
		for (int i = 0; i < k; ++i)
			perm.push_back(i);
		my_random_shuffle(perm.begin(), perm.end());
		int maxLoop = 10;
		int maxDist = max(1, k / 20);
		int maxLen = max(2, k / 20);
		int bestScoreHit = 1;
		int curScore = computeParsimonyPermutation(tree, alignment, params, permCol, perm);
		for (int loop = 1; loop <= maxLoop; ++loop)
		{
			// rand 2..maxLen
			int len = 2 + rand() % (maxLen - 1);
			// rand 0..k-len
			int l = rand() % (k - len);
			int r = l + len - 1;
			for (int i = 0; i < k; ++i)
			{
				if (l <= i && i <= r)
					continue;
				if (min(abs(i - l), abs(i - r)) > maxDist)
					continue;
				vector<int> newPerm;
				for (int j = 0; j <= i; ++j)
				{
					if (l <= j && j <= r)
						continue;
					newPerm.push_back(perm[j]);
				}
				for (int j = l; j <= r; ++j)
					newPerm.push_back(perm[j]);
				for (int j = i + 1; j < k; ++j)
				{
					if (l <= j && j <= r)
						continue;
					newPerm.push_back(perm[j]);
				}
				assert((int)newPerm.size() == k);
				int newScore = computeParsimonyPermutation(tree, alignment, params, permCol, newPerm);

				if (newScore < curScore)
				{
					curScore = newScore;
					perm = newPerm;
					bestScoreHit = 1;
				}
				else if (newScore == curScore)
				{
					if ((double)rand() / RAND_MAX >= 1.0 / bestScoreHit)
						perm = newPerm;
					++bestScoreHit;
				}
			}
		}

		candidates.push_back(perm);
		candidateScore.push_back(curScore);
		if ((int)candidates.size() > numCandidate)
		{
			int choice = 0;
			for (int i = 0; i < (int)candidates.size(); ++i)
			{
				if (candidateScore[i] > candidateScore[choice])
					choice = i;
				else if (candidateScore[i] == candidateScore[choice])
				{
					if (rand() % 2 == 1)
						choice = i;
				}
			}
			candidates.erase(candidates.begin() + choice);
			candidateScore.erase(candidateScore.begin() + choice);
		}
	}

	bestScore = updatePermutation(tree, alignment, params, permCol, candidates, candidateScore);
	cout << '\n';
	cout << "Best tree parsimony found after add more k rows using PLL core: ";
	cout << bestScore << '\n';
	cout << "Time: " << getCPUTime() - startTime << " seconds\n";
	return bestScore;
}

int updatePermutation(IQTree* tree, Alignment* alignment, Params& params, const vector<int> permCol, vector<vector<int> >& candidates, vector<int>& candidateScore)
{
	cout << "we are here\n";
	assert((int)candidates.size() > 0);
	// main loop
	int cnt = 0;
	while (true)
	{
		int id = rand() % ((int)candidates.size());
		vector<int> perm = candidates[id];
		int curScore = candidateScore[id];
		int k = (int)perm.size();
		int bestScoreHit = 1;
		int maxDist = max(1, k / 20);
		int maxLoop = 10;
		for (int loop = 1; loop <= maxLoop; ++loop)
		{
			for (int i = 0; i < k; ++i)
			{
				for (int j = i + 1; j <= min(k - 1, i + maxDist); ++j)
				{
					swap(perm[i], perm[j]);
					int newScore = computeParsimonyPermutation(tree, alignment, params, permCol, perm);
					if (newScore < curScore)
					{
						curScore = newScore;
						bestScoreHit = 1;
					}
					else if (newScore == curScore)
					{
						++bestScoreHit;
						if ((double)rand() / RAND_MAX < 1.0 / bestScoreHit) // not accept
							swap(perm[i], perm[j]);
					}
					else
						swap(perm[i], perm[j]);
				}
			}
		}

		if (curScore == candidateScore[id])
			++cnt;
		else
			cnt = 1;

		candidates.push_back(perm);
		candidateScore.push_back(curScore);
		// cout << curScore << " " << candidateScore[id] << " " << cnt << '\n';
		int choice = 0;
		for (int i = 0; i < (int)candidates.size(); ++i)
		{
			if (candidateScore[i] > candidateScore[choice])
				choice = i;
			else if (candidateScore[i] == candidateScore[choice])
			{
				if (rand() % 2 == 1)
					choice = i;
			}
		}
		candidates.erase(candidates.begin() + choice);
		candidateScore.erase(candidateScore.begin() + choice);

		if (cnt == 10)
			break;
	}

	int bestScore = candidateScore[0];
	for (int i = 0; i < (int)candidateScore.size(); ++i)
		bestScore = min(bestScore, candidateScore[i]);
	return bestScore;
}

int computeParsimonyPermutation(IQTree* tree, Alignment* alignment, Params& params, const vector<int>& permCol, const vector<int>& perm)
{
	IQTree newTree;
	(newTree).copyPhyloTree(tree);
	newTree.aln = new Alignment;
	newTree.aln->copyAlignment(tree->aln);
	int score = newTree.addRemainRowSPR(alignment->remainName, alignment->remainSeq, perm, permCol, params);
	delete newTree.aln;
	newTree.aln = NULL;
	return score;
}

void addMoreRowMutation(Params& params)
{
    Alignment* alignment;
    IQTree* tree;

    /****************** read in alignment **********************/
    if (params.partition_file)
    {
    	// Partition model analysis
    	if (params.partition_type)
    	{
    		// since nni5 does not work yet, stop the programm
    		if (params.nni5)
    			outError("-nni5 option is unsupported yet for proportitional partition model. please use -nni1 option");
    		if (params.aLRT_replicates)
    			outError("-alrt option is unsupported yet for proportitional partition model");
    		// initialize supertree - Proportional Edges case, "-spt p" option
    		tree = new PhyloSuperTreePlen(params);
    	}
    	else
    	{
    		// initialize supertree stuff if user specifies partition file with -sp option
    		tree = new PhyloSuperTree(params);
    	}
    	// this alignment will actually be of type SuperAlignment
    	alignment = tree->aln;
    }
    else if (params.maximum_parsimony && params.sankoff_cost_file)
    {
    	alignment = new Alignment(params.aln_file, params.sequence_type, params.intype);
    	tree = new ParsTree(alignment);
    	dynamic_cast<ParsTree*>(tree)->initParsData(&params);
    }
    else
    {
    	alignment = new Alignment(params.aln_file, params.sequence_type, params.intype, params.numStartRow);
    	if (params.maximum_parsimony && !params.sankoff_cost_file && params.condense_parsimony_equiv_sites)
    	{
    		Alignment* aln = new Alignment();
    		aln->condenseParsimonyEquivalentSites(alignment);
    		delete alignment;
    		alignment = aln;
    	}
    	if (params.analyze_alignment) {
    		alignment->analyzeAlignment();
    	}
    	tree = new IQTree(alignment);
    }

    //	if(params.maximum_parsimony && (params.gbo_replicates || params.sankoff_cost_file)){
    if (params.maximum_parsimony && (params.sort_alignment || params.sankoff_cost_file))
    {
    	optimizeAlignment(tree, params); // Diep: this is to rearrange columns for better speed in REPS
    }

    string original_model = params.model_name;

    vector<ModelInfo> model_info;
    alignment->checkGappySeq();
    
    cout << "\nStart placement core\n";

    auto startTime = getCPUTime();
    IQTree newTree;
    char* file_name = params.mutation_tree_file;
    bool is_rooted = false;
    newTree.readTree(file_name, is_rooted);

    newTree.setAlignment(tree->aln);
    newTree.aln = new Alignment;
    newTree.aln->copyAlignment(tree->aln);
    newTree.aln->ungroupSitePattern();
    newTree.save_branch_states_dad = new UINT[(newTree.aln->size() + 7) / 8 + 1];
    newTree.save_branch_fitch_result = new UINT[newTree.aln->size() + 1];
    newTree.add_row = true;
    newTree.aln->missingSamples = alignment->missingSamples;
    newTree.aln->existingSamples = alignment->existingSamples;
    newTree.aln->reference_nuc = alignment->reference_nuc;

    cout << "Tree parsimony before add k rows: " << newTree.computeParsimony() << '\n';
    // cout << "Ungroup alignment: " << tree->aln->getNSite() << " " << newTree.aln->getNSite() << '\n';
    vector<int> permCol = alignment->findPermCol();
    vector<int> savePermCol = permCol;
    vector<int> pos;

    if (alignment->existingSamples.size())
    {
        for (auto& m : permCol)
            m = alignment->existingSamples[0][m].position;
    }
    else
    {
        for (auto& m : permCol) m++;
    }

    for (int i = 0; i < (int)permCol.size(); ++i)
    {
        while ((int)pos.size() <= permCol[i])
            pos.push_back(0);
        pos[permCol[i]] = i;
    }

    for (int i = 0; i < (int)alignment->missingSamples.size(); ++i)
    {
        for (auto m : alignment->missingSamples[i])
        {
            assert(newTree.aln->reference_nuc[m.position] > 0);
        }
    }
    for (int i = 0; i < (int)alignment->existingSamples.size(); ++i)
    {
        for (auto m : alignment->existingSamples[i])
        {
            assert(newTree.aln->reference_nuc[m.position] > 0);
        }
    }

    int sz = 0;
    for (int m : permCol)
        sz = max(sz, m);
    newTree.cur_missing_sample_mutations.resize(sz + 1);
    newTree.cur_ancestral_mutations.resize(sz + 1);
    newTree.visited_missing_sample_mutations.resize(sz + 1);
    newTree.visited_ancestral_mutations.resize(sz + 1);

    newTree.initMutation(permCol);

    cout << "Tree parsimony after init mutations: " << newTree.computeParsimony() << " " << newTree.computeParsimonyScoreMutation() << '\n';
    int num_sample = (int)alignment->missingSamples.size();
    vector<MutationNode> missingSamples(num_sample);
    for (int i = 0; i < (int)alignment->missingSamples.size(); ++i)
    {
        missingSamples[i].mutations = alignment->missingSamples[i];
        for (auto m : alignment->missingSamples[i])
        {
            assert((m.ref_nuc & (m.ref_nuc - 1)) == 0);
        }
        missingSamples[i].name = alignment->missingSamples[i][0].name;
    }

    // newTree.checkMutation(pos);
    // cout << "correct mutations\n\n";

    for (int i = 0; i < min((int)missingSamples.size(), params.numAddRow); ++i)
    {
        vector<pair<PhyloNode*, PhyloNeighbor*> > bfs = newTree.breadth_first_expansion();
        size_t total_nodes = (int)bfs.size();
        // Stores the excess mutations to place the sample at each
        // node of the tree in DFS order. When placement is as a
        // child, it only contains parsimony-increasing mutations in
        // the sample. When placement is as a sibling, it contains
        // parsimony-increasing mutations as well as the mutations
        // on the placed node in common with the new sample. Note
        // guaranteed to be corrrect only for optimal nodes since
        // the mapper can terminate the search early for non-optimal
        // nodes
        std::vector<std::vector<Mutation>> node_excess_mutations(total_nodes);
        // Stores the imputed mutations for ambiguous bases in the
        // sampled in order to place the sample at each node of the
        // tree in DFS order. Again, guaranteed to be corrrect only
        // for pasrimony-optimal nodes
        std::vector<std::vector<Mutation>> node_imputed_mutations(total_nodes);

        bool best_node_has_unique = false;

        std::vector<bool> node_has_unique(total_nodes, false);
        size_t best_node_num_leaves = 0;
        int best_set_difference = (int)1e9 + 7;
        int set_difference = (int)1e9 + 7;
        size_t best_j = 0;
        size_t best_distance = (size_t)1e9 + 7;

        for (int j = 0; j < (int)bfs.size(); ++j)
        {
            // cout << "here\n";
            CandidateNode inp;
            inp.node = bfs[j].first;
            inp.node_branch = bfs[j].second;
            inp.missing_sample_mutations = &missingSamples[i].mutations;
            inp.excess_mutations = &node_excess_mutations[j];
            inp.imputed_mutations = &node_imputed_mutations[j];
            inp.best_node_num_leaves = &best_node_num_leaves;
            inp.best_set_difference = &best_set_difference;
            inp.set_difference = &set_difference;
            inp.best_j = &best_j;
            inp.best_distance = &best_distance;
            inp.j = j;
            inp.has_unique = &best_node_has_unique;
            inp.node_has_unique = &(node_has_unique);

            newTree.calculatePlacementMutation(pos, inp, false, true);
        }

        newTree.addNewSample(bfs[best_j].first, bfs[best_j].second, node_excess_mutations[best_j], i, missingSamples[i].name);
        // newTree.aln->addToAlignmentNewSeq(missingSamples[i].name, alignment->remainSeq[i], savePermCol);
        // newTree.checkMutation(pos);
        // cout << newTree.computeParsimonyScoreMutation() << " " << newTree.computeParsimonyScore() << '\n';
    }
    delete newTree.aln;
    newTree.aln = NULL;
    cout << "New tree's parsimony score: " << newTree.computeParsimonyScoreMutation() << '\n';
    cout << "Time: " << fixed << setprecision(3) << (double)(getCPUTime() - startTime) << " seconds\n";
}