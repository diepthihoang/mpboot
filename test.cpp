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
	ptree->setAlignment(&alignment); // make BIG difference $$$$$$$$$$$$

	// compute score by Sankoff algorithm
	cout << "Score = " << ptree->computeParsimony() << endl;
	cout << "Pattern score = ";
	ptree->printPatternScore();
	cout << endl;

	cout << "mst score     = ";
	for(int i = 0; i < ptree->aln->getNPattern(); i++)
//	for(int i = 0; i < 6; i++)
		cout << ptree->findMstScore(i) << ", ";
	cout << endl;
	delete ptree;

//	params.sankoff_cost_file = NULL;
//	IQTree * iqtree = new IQTree(&alignment);
//	iqtree->readTree(params.user_file, params.is_rooted);
//	cout << "IQTree score = " << iqtree->computeParsimony() << endl;
//
//	IQTree * iqtree2 = new IQTree;
//	iqtree2->readTree(params.user_file, params.is_rooted);
//	iqtree2->setAlignment(&alignment); // make BIG difference
//	cout << "IQTree2 score = " << iqtree2->computeParsimony() << endl;

//	delete iqtree;
//	delete iqtree2;

}
