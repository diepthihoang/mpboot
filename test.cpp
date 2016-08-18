#include "phylotree.h"
#include "test.h"
#include "alignment.h"
#include "parstree.h"

void test(Params &params){
	testWeightedParsimony(params);
}

void testWeightedParsimony(Params &params){
	if(!params.sankoff_cost_file) return;

	// read aln
	Alignment alignment(params.aln_file, params.sequence_type, params.intype);

	// initialize an ParsTree instance connecting with the alignment
	ParsTree * ptree = new ParsTree(&alignment);

	// initialize the cost matrix
	dynamic_cast<ParsTree *>(ptree)->initParsData(&params);

	// read in a tree from user's file
	ptree->readTree(params.user_file, params.is_rooted);

	// compute score by Sankoff algorithm
	cout << "Score = " << ptree->computeParsimony() << endl;

	delete ptree;
}
