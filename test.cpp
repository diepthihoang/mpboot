#include "phylotree.h"
#include "test.h"
#include "alignment.h"
#include "parstree.h"

void test(Params &params){
	testWeightedParsimony(params);
	// testTreeConvertTaxaToID(params);
    if(params.remove_dup_seq){
        testRemoveDuplicateSeq(params);
        cout << "\nDONE testRemoveDuplicateSeq(). Check output at " << params.user_file << ".\n\n";
        cout << "NOTE: If you want to remove duplicate sequences for MPBoot search algorithm, "
            << "the correct command is: \n"
            << "./mpboot -s " << params.aln_file << " -test_mode -remove_dup_seq " << params.user_file << "\n\n";
        cout << "NOTE: If you want to remove duplicate sequences for MPBoot bootstrap approximation algorithm, "
            << "the correct command is: \n"
            << "./mpboot -s " << params.aln_file << " -test_mode -bb 1000 -remove_dup_seq " << params.user_file << "\n\n";
    }
}

// -s <alnfile> -test_mode <treefile> -cost <costfile>
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

// -s <alnfile> -test_mode <treefile>
void testTreeConvertTaxaToID(Params &params){
	// read aln
	Alignment alignment(params.aln_file, params.sequence_type, params.intype);

	// initialize an ParsTree instance connecting with the alignment
	IQTree * ptree = new IQTree(&alignment);

	// read in a tree from user's file
	ptree->readTree(params.user_file, params.is_rooted);
	ptree->setAlignment(&alignment);

    string out_file;
    out_file = params.user_file;
    out_file += ".id";

	ptree->printTree(out_file.c_str(), WT_TAXON_ID | WT_SORT_TAXA);
	cout << "Please see result in " << out_file << endl;

	delete ptree;
}


// Usage 1: -s <alnfile> -test_mode <out_aln>
// Usage 2: -s <alnfile> -test_mode <out_aln> -bb 100 
// The 2nd will keep two identical sequences
void testRemoveDuplicateSeq(Params &params){
	// read aln
	Alignment * aln = new Alignment(params.aln_file, params.sequence_type, params.intype);

    StrVector removed_seqs;
    StrVector twin_seqs;
	Alignment *new_aln;
	if (params.root)
		new_aln = aln->removeIdenticalSeq((string)params.root, params.gbo_replicates > 0, removed_seqs, twin_seqs);
	else
		new_aln = aln->removeIdenticalSeq("", params.gbo_replicates > 0, removed_seqs, twin_seqs);
	if (removed_seqs.size() > 0) {
		cout << "NOTE: " << removed_seqs.size() << " identical sequences will be ignored when print output alignment." << endl;
		if (verbose_mode >= VB_MED) {
			for (int i = 0; i < removed_seqs.size(); i++) {
				cout << removed_seqs[i] << " is identical to " << twin_seqs[i] << endl;
			}
		}
        delete aln; // REF: fix in iqtree1.6.12
		aln = new_aln;
	}
    aln->printPhylip(params.user_file);
    delete aln;
}
