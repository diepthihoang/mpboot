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

#include <algorithm>

void reportReferences(Params& params, ofstream& out, string& original_model)
{
	out << "To cite IQ-TREE please use:" << endl
		<< endl
		<< "Lam-Tung Nguyen, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh (2014)" << endl
		<< "IQ-TREE: : A fast and effective stochastic algorithm for estimating" << endl
		<< "maximum likelihood phylogenies. Mol. Biol. Evol., in press." << endl
		<< endl;

	if (params.gbo_replicates)
		out << "Since you also used ultrafast bootstrap (UFBoot) please cite: " << endl
		<< endl
		<< "Bui Quang Minh, Minh Anh Thi Nguyen, and Arndt von Haeseler (2013) Ultrafast" << endl
		<< "approximation for phylogenetic bootstrap. Mol. Biol. Evol., 30:1188-1195." << endl
		<< endl;

	/*		"*** If you use the parallel version, please cite: " << endl << endl <<
	 "Bui Quang Minh, Le Sy Vinh, Arndt von Haeseler, and Heiko A. Schmidt (2005)" << endl <<
	 "pIQPNNI - parallel reconstruction of large maximum likelihood phylogenies." << endl <<
	 "Bioinformatics, 21:3794-3796." << endl << endl;*/

	 // 	if (original_model == "TEST" || original_model == "TESTONLY")
	 // 		out << "Since you used Modeltest please also cite Posada and Crandall (1998)" << endl << endl;
}

void reportAlignment(ofstream& out, Alignment& alignment)
{
	out << "Input data: " << alignment.getNSeq() << " sequences with "
		<< alignment.getNSite() << " "
		<< ((alignment.seq_type == SEQ_BINARY) ? "binary" : ((alignment.seq_type == SEQ_DNA) ? "nucleotide" : (alignment.seq_type == SEQ_PROTEIN) ? "amino-acid"
			: (alignment.seq_type == SEQ_CODON) ? "codon"
			: "morphological"))
		<< " sites" << endl
		<< "Number of constant sites: "
		<< round(alignment.frac_const_sites * alignment.getNSite())
		<< " (= " << alignment.frac_const_sites * 100 << "% of all sites)"
		<< endl
		<< "Number of site patterns: " << alignment.size() << endl
		<< endl;
}

void pruneModelInfo(vector<ModelInfo>& model_info, PhyloSuperTree* tree)
{
	vector<ModelInfo> res_info;
	for (vector<PartitionInfo>::iterator it = tree->part_info.begin(); it != tree->part_info.end(); it++)
	{
		for (vector<ModelInfo>::iterator mit = model_info.begin(); mit != model_info.end(); mit++)
			if (mit->set_name == it->name)
				res_info.push_back(*mit);
	}
	model_info = res_info;
}

void reportModelSelection(ofstream& out, Params& params, vector<ModelInfo>& model_info, bool is_partitioned)
{
	out << "Best-fit model according to "
		<< ((params.model_test_criterion == MTC_BIC) ? "BIC" : ((params.model_test_criterion == MTC_AIC) ? "AIC" : "AICc")) << ": ";
	vector<ModelInfo>::iterator it;
	if (is_partitioned)
	{
		string set_name = "";
		for (it = model_info.begin(); it != model_info.end(); it++)
		{
			if (it->set_name != set_name)
			{
				if (set_name != "")
					out << ",";
				out << it->name << ":" << it->set_name;
				set_name = it->set_name;
			}
		}
	}
	else
	{
		out << model_info[0].name;
	}

	if (is_partitioned)
	{
		out << endl
			<< endl
			<< "List of best-fit models per partition:" << endl
			<< endl;
	}
	else
	{
		out << endl
			<< endl
			<< "List of models sorted by "
			<< ((params.model_test_criterion == MTC_BIC) ? "BIC" : ((params.model_test_criterion == MTC_AIC) ? "AIC" : "AICc"))
			<< " scores: " << endl
			<< endl;
	}
	if (is_partitioned)
		out << "  ID  ";
	out << "Model             LogL          AIC      w-AIC      AICc     w-AICc       BIC      w-BIC" << endl;
	/*
	if (is_partitioned)
		out << "----------";

	out << "----------------------------------------------------------------------------------------" << endl;
	*/
	int setid = 1;
	for (it = model_info.begin(); it != model_info.end(); it++)
	{
		if (it->AIC_score == DBL_MAX)
			continue;
		if (it != model_info.begin() && it->set_name != (it - 1)->set_name)
			setid++;
		if (is_partitioned && it != model_info.begin() && it->set_name == (it - 1)->set_name)
			continue;
		if (is_partitioned)
		{
			out.width(4);
			out << right << setid << "  ";
		}
		out.width(13);
		out << left << it->name << " ";
		out.width(11);
		out << right << it->logl << " ";
		out.width(11);
		out << it->AIC_score << ((it->AIC_conf) ? " + " : " - ") << it->AIC_weight << " ";
		out.width(11);
		out << it->AICc_score << ((it->AICc_conf) ? " + " : " - ") << it->AICc_weight << " ";
		out.width(11);
		out << it->BIC_score << ((it->BIC_conf) ? " + " : " - ") << it->BIC_weight;
		out << endl;
	}
	out << endl;
	out << "AIC, w-AIC   : Akaike information criterion scores and weights." << endl
		<< "AICc, w-AICc : Corrected AIC scores and weights." << endl
		<< "BIC, w-BIC   : Bayesian information criterion scores and weights." << endl
		<< endl

		<< "Plus signs denote the 95% confidence sets." << endl
		<< "Minus signs denote significant exclusion." << endl;
	out << endl;
}

void reportModel(ofstream& out, PhyloTree& tree)
{
	int i, j, k;
	out << "Model of substitution: " << tree.getModelName() << endl
		<< endl;

	if (tree.aln->num_states <= 4)
	{
		out << "Rate parameter R:" << endl
			<< endl;

		double* rate_mat = new double[tree.aln->num_states * tree.aln->num_states];
		if (!tree.getModel()->isSiteSpecificModel())
			tree.getModel()->getRateMatrix(rate_mat);
		else
			((ModelSet*)tree.getModel())->front()->getRateMatrix(rate_mat);
		if (tree.aln->num_states > 4)
			out << fixed;
		if (tree.getModel()->isReversible())
		{
			for (i = 0, k = 0; i < tree.aln->num_states - 1; i++)
				for (j = i + 1; j < tree.aln->num_states; j++, k++)
				{
					out << "  " << tree.aln->convertStateBackStr(i) << "-" << tree.aln->convertStateBackStr(j) << ": "
						<< rate_mat[k];
					if (tree.aln->num_states <= 4)
						out << endl;
					else if (k % 5 == 4)
						out << endl;
				}
		}
		else
		{ // non-reversible model
			for (i = 0, k = 0; i < tree.aln->num_states; i++)
				for (j = 0; j < tree.aln->num_states; j++)
					if (i != j)
					{
						out << "  " << tree.aln->convertStateBackStr(i) << "-" << tree.aln->convertStateBackStr(j)
							<< ": " << rate_mat[k];
						if (tree.aln->num_states <= 4)
							out << endl;
						else if (k % 5 == 4)
							out << endl;
						k++;
					}
		}

		// if (tree.aln->num_states > 4)
		out << endl;
		out.unsetf(ios_base::fixed);
		delete[] rate_mat;
	}
	out << "State frequencies: ";
	if (tree.getModel()->isSiteSpecificModel())
		out << "(site specific frequencies)" << endl
		<< endl;
	else
	{
		if (!tree.getModel()->isReversible())
			out << "(inferred from Q matrix)" << endl;
		else
			switch (tree.getModel()->getFreqType())
			{
			case FREQ_EMPIRICAL:
				out << "(empirical counts from alignment)" << endl;
				break;
			case FREQ_ESTIMATE:
				out << "(estimated with maximum likelihood)" << endl;
				break;
			case FREQ_USER_DEFINED:
				out << "(user-defined)" << endl;
				break;
			case FREQ_EQUAL:
				out << "(equal frequencies)" << endl;
				break;
			default:
				break;
			}
		out << endl;

		double* state_freqs = new double[tree.aln->num_states];
		tree.getModel()->getStateFrequency(state_freqs);
		for (i = 0; i < tree.aln->num_states; i++)
			out << "  pi(" << tree.aln->convertStateBackStr(i) << ") = "
			<< state_freqs[i] << endl;
		delete[] state_freqs;
		out << endl;
		if (tree.aln->num_states <= 4)
		{
			// report Q matrix
			double* q_mat = new double[tree.aln->num_states * tree.aln->num_states];
			tree.getModel()->getQMatrix(q_mat);

			out << "Rate matrix Q:" << endl
				<< endl;
			for (i = 0, k = 0; i < tree.aln->num_states; i++)
			{
				out << "  " << tree.aln->convertStateBackStr(i);
				for (j = 0; j < tree.aln->num_states; j++, k++)
				{
					out << "  ";
					out.width(8);
					out << q_mat[k];
				}
				out << endl;
			}
			out << endl;
			delete[] q_mat;
		}
	}
}

void reportRate(ofstream& out, PhyloTree& tree)
{
	int i;
	RateHeterogeneity* rate_model = tree.getRate();
	out << "Model of rate heterogeneity: " << rate_model->full_name << endl;
	rate_model->writeInfo(out);

	if (rate_model->getNDiscreteRate() > 1 || rate_model->getPInvar() > 0.0)
	{
		out << endl
			<< " Category  Relative_rate  Proportion" << endl;
		if (rate_model->getPInvar() > 0.0)
			out << "  0         0              " << rate_model->getPInvar()
			<< endl;
		int cats = rate_model->getNDiscreteRate();
		DoubleVector prop;
		if (rate_model->getGammaShape() > 0 || rate_model->getPtnCat(0) < 0)
		{
			//			prop.resize(cats, (1.0 - rate_model->getPInvar()) / rate_model->getNRate());
			prop.resize(cats);
			for (i = 0; i < cats; i++)
				prop[i] = rate_model->getProp(i);
		}
		else
		{
			prop.resize(cats, 0.0);
			for (i = 0; i < tree.aln->getNPattern(); i++)
				prop[rate_model->getPtnCat(i)] += tree.aln->at(i).frequency;
			for (i = 0; i < cats; i++)
				prop[i] /= tree.aln->getNSite();
		}
		for (i = 0; i < cats; i++)
		{
			out << "  " << i + 1 << "         ";
			out.width(14);
			out << left << rate_model->getRate(i) << " " << prop[i];
			out << endl;
		}
		if (rate_model->getGammaShape() > 0)
		{
			out << "Relative rates are computed as " << ((dynamic_cast<RateGamma*>(rate_model)->isCutMedian()) ? "MEDIAN" : "MEAN") << " of the portion of the Gamma distribution falling in the category." << endl;
		}
	}
	/*
	 if (rate_model->getNDiscreteRate() > 1 || rate_model->isSiteSpecificRate())
	 out << endl << "See file " << rate_file << " for site-specific rates and categories" << endl;*/
	out << endl;
}

void reportTree(ofstream& out, Params& params, PhyloTree& tree, double tree_lh,
	double lh_variance)
{
	double epsilon = 1.0 / tree.getAlnNSite();
	double totalLen = tree.treeLength();
	if (!params.maximum_parsimony)
	{
		out << "Total tree length (sum of branch lengths): " << totalLen << endl;
		double totalLenInternal = tree.treeLengthInternal(epsilon);
		out << "Sum of internal branch lengths: " << totalLenInternal << endl;
		out << "Sum of internal branch lengths divided by total tree length: "
			<< totalLenInternal / totalLen << endl;
		out << endl;
		// out << "ZERO BRANCH EPSILON = " << epsilon << endl;
		int zero_internal_branches = tree.countZeroInternalBranches(NULL, NULL, epsilon);
		if (zero_internal_branches > 0)
		{
			// int zero_internal_branches = tree.countZeroInternalBranches(NULL, NULL, epsilon);
			/*
			out << "WARNING: " << zero_branches
					<< " branches of near-zero lengths (<" << epsilon << ") and should be treated with caution!"
					<< endl;
			*/
			out << "WARNING: " << zero_internal_branches
				<< " near-zero internal branches (<" << epsilon << ") should be treated with caution"
				<< endl;
			/*
			cout << endl << "WARNING: " << zero_branches
					<< " branches of near-zero lengths (<" << epsilon << ") and should be treated with caution!"
					<< endl;
			*/
			out << "         Such branches are denoted by '**' in the figure"
				<< endl
				<< endl;
		}
		int long_branches = tree.countLongBranches(NULL, NULL, MAX_BRANCH_LEN - 0.2);
		if (long_branches > 0)
		{
			// stringstream sstr;
			out << "WARNING: " << long_branches
				<< " too long branches (>" << MAX_BRANCH_LEN - 0.2 << ") should be treated with caution!"
				<< endl;
			// out << sstr.str();
			// cout << sstr.str();
		}
	}
	tree.sortTaxa();
	// tree.setExtendedFigChar();
	tree.drawTree(out, WT_BR_SCALE, epsilon);
	int df = tree.getModelFactory()->getNParameters();
	int ssize = tree.getAlnNSite();
	double AIC_score, AICc_score, BIC_score;
	computeInformationScores(tree_lh, df, ssize, AIC_score, AICc_score, BIC_score);

	if (params.maximum_parsimony)
	{
		out << "Parsimony score of the tree: " << fixed << int(-tree_lh) << endl;
	}
	else
	{
		out << "Log-likelihood of the tree: " << fixed << tree_lh << " (s.e. "
			<< sqrt(lh_variance) << ")" << endl
			<< "Number of free parameters: " << df << endl
			<< "Akaike information criterion (AIC) score: " << AIC_score << endl
			<< "Corrected Akaike information criterion (AICc) score: " << AICc_score << endl
			<< "Bayesian information criterion (BIC) score: " << BIC_score << endl
			<< "Unconstrained log-likelihood (without tree): "
			<< tree.aln->computeUnconstrainedLogL() << endl
			<< endl;
		//<< "Total tree length: " << tree.treeLength() << endl << endl
	}
	out << "Tree in newick format:" << endl
		<< endl;

	//	tree.printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA); // ML
	tree.printTree(out, WT_SORT_TAXA); // Diep: for MP

	out << endl
		<< endl;
}

void reportCredits(ofstream& out)
{
	out << "CREDITS" << endl
		<< "-------" << endl
		<< endl
		<< "Some parts of the code were taken from the following packages/libraries:"
		<< endl
		<< endl
		<< "Schmidt HA, Strimmer K, Vingron M, and von Haeseler A (2002)" << endl
		<< "TREE-PUZZLE: maximum likelihood phylogenetic analysis using quartets" << endl
		<< "and parallel computing. Bioinformatics, 18(3):502-504." << endl
		<< endl

		//<< "The source code to construct the BIONJ tree were taken from BIONJ software:"
		//<< endl << endl
		<< "Gascuel O (1997) BIONJ: an improved version of the NJ algorithm" << endl
		<< "based on a simple model of sequence data. Mol. Bio. Evol., 14:685-695." << endl
		<< endl

		//<< "The Nexus file parser was taken from the Nexus Class Library:"
		//<< endl << endl
		<< "Paul O. Lewis (2003) NCL: a C++ class library for interpreting data files in" << endl
		<< "NEXUS format. Bioinformatics, 19(17):2330-2331." << endl
		<< endl

		<< "Mascagni M and Srinivasan A (2000) Algorithm 806: SPRNG: A Scalable Library" << endl
		<< "for Pseudorandom Number Generation. ACM Transactions on Mathematical Software," << endl
		<< "26: 436-461." << endl
		<< endl

		<< "Guennebaud G, Jacob B, et al. (2010) Eigen v3. http://eigen.tuxfamily.org" << endl
		<< endl;
	/*
	<< "The Modeltest 3.7 source codes were taken from:" << endl << endl
	<< "David Posada and Keith A. Crandall (1998) MODELTEST: testing the model of"
	<< endl << "DNA substitution. Bioinformatics, 14(9):817-8." << endl
	*/
}

/***********************************************************
 * CREATE REPORT FILE
 ***********************************************************/
extern StringIntMap pllTreeCounter;
void reportPhyloAnalysis(Params& params, string& original_model,
	Alignment& alignment, IQTree& tree, vector<ModelInfo>& model_info,
	StrVector& removed_seqs, StrVector& twin_seqs)
{
	if (params.count_trees)
	{
		// addon: print #distinct trees
		cout << endl
			<< "NOTE: " << pllTreeCounter.size() << " distinct trees evaluated during whole tree search" << endl;

		IntVector counts;
		for (StringIntMap::iterator i = pllTreeCounter.begin(); i != pllTreeCounter.end(); i++)
		{
			if (i->second > counts.size())
				counts.resize(i->second + 1, 0);
			counts[i->second]++;
		}
		for (IntVector::iterator i2 = counts.begin(); i2 != counts.end(); i2++)
		{
			if (*i2 != 0)
			{
				cout << "#Trees occuring " << (i2 - counts.begin()) << " times: " << *i2 << endl;
			}
		}
	}
	string outfile = params.out_prefix;

	//	outfile += ".iqtree";
	outfile += ".mpboot";
	try
	{
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(outfile.c_str());
		out << "MPBoot " << iqtree_VERSION_MAJOR << "." << iqtree_VERSION_MINOR
			<< "." << iqtree_VERSION_PATCH << " built " << __DATE__ << endl
			<< endl;
		if (params.partition_file)
			out << "Partition file name: " << params.partition_file << endl;
		if (params.aln_file)
			out << "Input file name: " << params.aln_file << endl;

		if (params.user_file)
			out << "User tree file name: " << params.user_file << endl;
		out << "Type of analysis: ";
		if (params.maximum_parsimony)
		{
			out << "Parsimony tree reconstruction";
			if (params.num_bootstrap_samples > 0)
			{
				out << " + ";
				out << "non-parametric bootstrap (" << params.num_bootstrap_samples
					<< " replicates)";
			}
			if (params.gbo_replicates > 0)
			{
				out << " + ";
				out << "MPBoot bootstrap approximation (" << params.gbo_replicates
					<< " replicates)";
			}
		}
		else
		{
			if (params.compute_ml_tree)
				out << "tree reconstruction";
			if (params.num_bootstrap_samples > 0)
			{
				if (params.compute_ml_tree)
					out << " + ";
				out << "non-parametric bootstrap (" << params.num_bootstrap_samples
					<< " replicates)";
			}
		}
		out << endl;
		out << "Random seed number: " << params.ran_seed << endl
			<< endl;
		// Diep: the two lines below do not need for UFBoot-MP
		//		out << "REFERENCES" << endl << "----------" << endl << endl;
		//		reportReferences(params, out, original_model);

		out << "SEQUENCE ALIGNMENT" << endl
			<< "------------------" << endl
			<< endl;
		if (tree.isSuperTree())
		{
			out << "Input data: " << alignment.getNSeq() << " taxa with "
				<< alignment.getNSite() << " partitions and "
				<< tree.getAlnNSite() << " total sites ("
				<< ((SuperAlignment*)tree.aln)->computeMissingData() * 100 << "% missing data)" << endl
				<< endl;

			PhyloSuperTree* stree = (PhyloSuperTree*)&tree;
			int namelen = stree->getMaxPartNameLength();
			int part;
			out.width(max(namelen + 6, 10));
			out << left << "  ID  Name"
				<< "  Type  #Seqs  #Sites  #Patterns  #Const_Sites" << endl;
			// out << string(namelen+54, '-') << endl;
			part = 0;
			for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++, part++)
			{
				// out << "FOR PARTITION " << stree->part_info[part].name << ":" << endl << endl;
				// reportAlignment(out, *((*it)->aln));
				out.width(4);
				out << right << part + 1 << "  ";
				out.width(max(namelen, 4));
				out << left << stree->part_info[part].name << "  ";
				out.width(6);
				switch ((*it)->aln->seq_type)
				{
				case SEQ_BINARY:
					out << "BIN";
					break;
				case SEQ_CODON:
					out << "CODON";
					break;
				case SEQ_DNA:
					out << "DNA";
					break;
				case SEQ_MORPH:
					out << "MORPH";
					break;
				case SEQ_MULTISTATE:
					out << "TINA";
					break;
				case SEQ_PROTEIN:
					out << "AA";
					break;
				case SEQ_UNKNOWN:
					out << "???";
					break;
				}
				out.width(5);
				out << right << (*it)->aln->getNSeq() << "  ";
				out.width(6);
				out << (*it)->aln->getNSite() << "  ";
				out.width(6);
				out << (*it)->aln->getNPattern() << "      ";
				out << round((*it)->aln->frac_const_sites * 100) << "%" << endl;
			}
			out << endl;
		}
		else
			reportAlignment(out, alignment);

		out.precision(4);
		out << fixed;

		if (!model_info.empty())
		{
			out << "MODEL SELECTION" << endl
				<< "---------------" << endl
				<< endl;
			if (tree.isSuperTree())
				pruneModelInfo(model_info, (PhyloSuperTree*)&tree);
			reportModelSelection(out, params, model_info, tree.isSuperTree());
		}

		if (!params.maximum_parsimony)
		{
			out << "SUBSTITUTION PROCESS" << endl
				<< "--------------------" << endl
				<< endl;
			if (tree.isSuperTree())
			{
				if (params.partition_type)
					out << "Proportional partition model with joint branch lengths and separate models between partitions" << endl
					<< endl;
				else
					out << "Full partition model with separate branch lengths and models between partitions" << endl
					<< endl;
				PhyloSuperTree* stree = (PhyloSuperTree*)&tree;
				PhyloSuperTree::iterator it;
				int part;
				if (params.partition_type)
					out << "  ID  Model          Rate   Parameters" << endl;
				else
					out << "  ID  Model          Parameters" << endl;
				// out << "-------------------------------------" << endl;
				for (it = stree->begin(), part = 0; it != stree->end(); it++, part++)
				{
					out.width(4);
					out << right << (part + 1) << "  ";
					out.width(14);
					if (params.partition_type)
						out << left << (*it)->getModelName() << " " << stree->part_info[part].part_rate << " " << (*it)->getModelNameParams() << endl;
					else
						out << left << (*it)->getModelName() << " " << (*it)->getModelNameParams() << endl;
				}
				out << endl;
				/*
				for (it = stree->begin(), part = 0; it != stree->end(); it++, part++) {
					reportModel(out, *(*it));
					reportRate(out, *(*it));
				}*/
			}
			else
			{
				reportModel(out, tree);
				reportRate(out, tree);
			}
		}

		/*
		out << "RATE HETEROGENEITY" << endl << "------------------" << endl
				<< endl;
		if (tree.isSuperTree()) {
			PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
			int part = 0;
			for (PhyloSuperTree::iterator it = stree->begin();
					it != stree->end(); it++, part++) {
				out << "FOR PARTITION " << stree->part_info[part].name << ":"
						<< endl << endl;
				reportRate(out, *(*it));
			}
		} else
			reportRate(out, tree);
		*/
		// Bootstrap analysis:
		// Display as outgroup: a

		if (original_model == "WHTEST")
		{
			out << "TEST OF MODEL HOMOGENEITY" << endl
				<< "-------------------------" << endl
				<< endl;
			out << "Delta of input data:                 "
				<< params.whtest_delta << endl;
			out << ".95 quantile of Delta distribution:  "
				<< params.whtest_delta_quantile << endl;
			out << "Number of simulations performed:     "
				<< params.whtest_simulations << endl;
			out << "P-value:                             "
				<< params.whtest_p_value << endl;
			if (params.whtest_p_value < 0.05)
			{
				out
					<< "RESULT: Homogeneity assumption is rejected (p-value cutoff 0.05)"
					<< endl;
			}
			else
			{
				out
					<< "RESULT: Homogeneity assumption is NOT rejected (p-value cutoff 0.05)"
					<< endl;
			}
			out << endl
				<< "*** For this result please cite:" << endl
				<< endl;
			out
				<< "G. Weiss and A. von Haeseler (2003) Testing substitution models"
				<< endl
				<< "within a phylogenetic tree. Mol. Biol. Evol, 20(4):572-578"
				<< endl
				<< endl;
		}
		/*
				out << "TREE SEARCH" << endl << "-----------" << endl << endl
						<< "Stopping rule: "
						<< ((params.stop_condition == SC_STOP_PREDICT) ? "Yes" : "No")
						<< endl << "Number of iterations: "
						<< tree.stop_rule.getNumIterations() << endl
						<< "Probability of deleting sequences: " << params.p_delete
						<< endl << "Number of representative leaves: "
						<< params.k_representative << endl
						<< "NNI log-likelihood cutoff: " << tree.getNNICutoff() << endl
						<< endl;
		*/
		if (params.compute_ml_tree)
		{
			out << (params.maximum_parsimony ? "MAXIMUM PARSIMONY TREE" : "MAXIMUM LIKELIHOOD TREE") << endl
				<< "-----------------------" << endl
				<< endl;

			tree.setRootNode(params.root);
			out << "NOTE: Tree is UNROOTED although outgroup taxon '" << tree.root->name << "' is drawn at root" << endl;
			if (params.partition_file)
				out << "NOTE: Branch lengths are weighted average over all partitions"
				<< endl
				<< "      (weighted by the number of sites in the partitions)"
				<< endl;
			if (params.aLRT_replicates > 0 || params.gbo_replicates || (params.num_bootstrap_samples && params.compute_ml_tree))
			{
				out << "Numbers in parentheses are ";
				if (params.aLRT_replicates > 0)
					out << "SH-aLRT supports";
				if (params.num_bootstrap_samples && params.compute_ml_tree)
				{
					if (params.aLRT_replicates > 0)
						out << " /";
					out << " standard bootstrap supports";
				}
				if (params.gbo_replicates)
				{
					if (params.aLRT_replicates > 0)
						out << " /";
					out << " ultrafast bootstrap supports";
				}
				out << " (%)" << endl;
			}
			out << endl;
			reportTree(out, params, tree, tree.getBestScore(), tree.logl_variance);

			if (tree.isSuperTree())
			{
				PhyloSuperTree* stree = (PhyloSuperTree*)&tree;
				stree->mapTrees();
				int empty_branches = stree->countEmptyBranches();
				if (empty_branches)
				{
					stringstream ss;
					ss << empty_branches << " branches in the overall tree with no phylogenetic information due to missing data!";
					outWarning(ss.str());
				}
				/*
				int part = 0;
				for (PhyloSuperTree::iterator it = stree->begin();
						it != stree->end(); it++, part++) {
					out << "FOR PARTITION " << stree->part_info[part].name
							<< ":" << endl << endl;
					string root_name;
					if (params.root)
						root_name = params.root;
					else
						root_name = (*it)->aln->getSeqName(0);
					(*it)->root = (*it)->findNodeName(root_name);
					assert((*it)->root);
					reportTree(out, params, *(*it), (*it)->computeLikelihood(),
							(*it)->computeLogLVariance());
				}*/
			}
		}
		/*
		 if (params.write_intermediate_trees) {
		 out << endl << "CONSENSUS OF INTERMEDIATE TREES" << endl << "-----------------------" << endl << endl
		 << "Number of intermediate trees: " << tree.stop_rule.getNumIterations() << endl
		 << "Split threshold: " << params.split_threshold << endl
		 << "Burn-in: " << params.tree_burnin << endl << endl;
		 }*/

		if (params.consensus_type == CT_CONSENSUS_TREE)
		{
			out << "CONSENSUS TREE" << endl
				<< "--------------" << endl
				<< endl;
			out << "Consensus tree is constructed from "
				<< (params.num_bootstrap_samples ? params.num_bootstrap_samples : params.gbo_replicates)
				<< " bootstrap trees" << endl
				<< "Branches with bootstrap support >"
				<< floor(params.split_threshold * 1000) / 10 << "% are kept";
			if (params.split_threshold == 0.0)
				out << " (extended consensus)";
			if (params.split_threshold == 0.5)
				out << " (majority-rule consensus)";
			if (params.split_threshold >= 0.99)
				out << " (strict consensus)";

			if (!params.maximum_parsimony)
				out << endl
				<< "Branch lengths are optimized by maximum likelihood on original alignment" << endl;
			else
				out << endl;
			out << "Numbers in parentheses are bootstrap supports (%)" << endl
				<< endl;

			string con_file = params.out_prefix;
			con_file += ".contree";
			bool rooted = false;
			MTree contree;
			contree.readTree(con_file.c_str(), rooted);
			contree.drawTree(out, WT_BR_SCALE);
			out << endl
				<< "Consensus tree in newick format: " << endl
				<< endl;
			//			contree.printTree(out); // for ML
			contree.printTree(out, WT_SORT_TAXA); // Diep: for MP
			out << endl
				<< endl;
			//			tree.freeNode();
			//			tree.root = NULL;
			//			tree.readTree(con_file.c_str(), rooted);
			//			if (removed_seqs.size() > 0) {
			//				tree.reinsertIdenticalSeqs(tree.aln, removed_seqs, twin_seqs);
			//			}
			//			tree.setAlignment(tree.aln);

			// bug fix
			//			if ((tree.sse == LK_EIGEN || tree.sse == LK_EIGEN_SSE) && !tree.isBifurcating()) {
			//				cout << "NOTE: Changing to old kernel as consensus tree is multifurcating" << endl;
			//				tree.changeLikelihoodKernel(LK_SSE);
			//			}

			//			tree.initializeAllPartialLh();
			//			tree.fixNegativeBranch(false);
			//			if (tree.isSuperTree())
			//				((PhyloSuperTree*) &tree)->mapTrees();
			//			tree.optimizeAllBranches();
			//			tree.printTree(con_file.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
			//			tree.sortTaxa();
			//			tree.drawTree(out, WT_BR_SCALE);
			//			out << endl << "Consensus tree in newick format: " << endl << endl;
			//			tree.printResultTree(out);
			//			out << endl << endl;
		}

		/* evaluate user trees */
		vector<TreeInfo> info;
		IntVector distinct_trees;
		if (params.treeset_file)
		{
			evaluateTrees(params, &tree, info, distinct_trees);
			out.precision(4);

			out << endl
				<< "USER TREES" << endl
				<< "----------" << endl
				<< endl;
			out << "See " << params.treeset_file << ".trees for trees with branch lengths." << endl
				<< endl;
			if (params.topotest_replicates && info.size() > 1)
			{
				if (params.do_weighted_test)
				{
					out << "Tree      logL    deltaL  bp-RELL    p-KH     p-SH    p-WKH    p-WSH    c-ELW" << endl;
					out << "-------------------------------------------------------------------------------" << endl;
				}
				else
				{
					out << "Tree      logL    deltaL  bp-RELL    p-KH     p-SH    c-ELW" << endl;
					out << "-------------------------------------------------------------" << endl;
				}
			}
			else
			{
				out << "Tree      logL    deltaL" << endl;
				out << "-------------------------" << endl;
			}
			double maxL = -DBL_MAX;
			int tid, orig_id;
			for (tid = 0; tid < info.size(); tid++)
				if (info[tid].logl > maxL)
					maxL = info[tid].logl;
			for (orig_id = 0, tid = 0; orig_id < distinct_trees.size(); orig_id++)
			{
				out.width(3);
				out << right << orig_id + 1 << " ";
				if (distinct_trees[orig_id] >= 0)
				{
					out << " = tree " << distinct_trees[orig_id] + 1 << endl;
					continue;
				}
				out.precision(3);
				out.width(12);
				out << info[tid].logl << " ";
				out.width(7);
				out << maxL - info[tid].logl;
				if (!params.topotest_replicates || info.size() <= 1)
				{
					out << endl;
					tid++;
					continue;
				}
				out.precision(4);
				out << "  ";
				out.width(6);
				out << info[tid].rell_bp;
				if (info[tid].rell_confident)
					out << " + ";
				else
					out << " - ";
				out.width(6);
				out << right << info[tid].kh_pvalue;
				if (info[tid].kh_pvalue < 0.05)
					out << " - ";
				else
					out << " + ";
				out.width(6);
				out << right << info[tid].sh_pvalue;
				if (info[tid].sh_pvalue < 0.05)
					out << " - ";
				else
					out << " + ";
				if (params.do_weighted_test)
				{
					out.width(6);
					out << right << info[tid].wkh_pvalue;
					if (info[tid].wkh_pvalue < 0.05)
						out << " - ";
					else
						out << " + ";
					out.width(6);
					out << right << info[tid].wsh_pvalue;
					if (info[tid].wsh_pvalue < 0.05)
						out << " - ";
					else
						out << " + ";
				}
				out.width(6);
				out << info[tid].elw_value;
				if (info[tid].elw_confident)
					out << " +";
				else
					out << " -";
				out << endl;
				tid++;
			}
			out << endl;

			if (params.topotest_replicates)
			{
				out << "deltaL  : logL difference from the maximal logl in the set." << endl
					<< "bp-RELL : bootstrap proportion using RELL method (Kishino et al. 1990)." << endl
					<< "p-KH    : p-value of one sided Kishino-Hasegawa test (1989)." << endl
					<< "p-SH    : p-value of Shimodaira-Hasegawa test (2000)." << endl;
				if (params.do_weighted_test)
				{
					out << "p-WKH   : p-value of weighted KH test." << endl
						<< "p-WSH   : p-value of weighted SH test." << endl;
				}
				out << "c-ELW   : Expected Likelihood Weight (Strimmer & Rambaut 2002)." << endl
					<< endl
					<< "Plus signs denote the 95% confidence sets." << endl
					<< "Minus signs denote significant exclusion." << endl
					<< "All tests performed "
					<< params.topotest_replicates << " resamplings using the RELL method." << endl;
			}
			out << endl;
		}

		time_t cur_time;
		time(&cur_time);

		char* date_str;
		date_str = ctime(&cur_time);
		out.unsetf(ios_base::fixed);
		out << "TIME STAMP" << endl
			<< "----------" << endl
			<< endl
			<< "Date and time: " << date_str << "Total CPU time used: "
			<< (double)params.run_time << " seconds (" << convert_time(params.run_time) << ")" << endl
			<< "Total wall-clock time used: " << getRealTime() - params.start_real_time
			<< " seconds (" << convert_time(getRealTime() - params.start_real_time) << ")" << endl
			<< endl;

		// reportCredits(out); // not needed, now in the manual
		out.close();
	}
	catch (ios::failure)
	{
		outError(ERR_WRITE_OUTPUT, outfile);
	}

	cout << endl
		<< "Analysis results written to: " << endl
		<< "  MPBoot report:                 " << params.out_prefix << ".mpboot"
		<< endl;
	if (params.compute_ml_tree)
	{
		cout << "  Maximum-" << (params.maximum_parsimony ? "parsimony " : "likelihood") << " tree:       " << params.out_prefix
			<< ".treefile" << endl;
		if (params.write_local_optimal_trees)
			cout << "  Candidate trees (" << tree.candidateTrees.size() << "):      " << params.out_prefix << ".localtrees" << endl;
	}
	if (!params.user_file && params.start_tree == STT_BIONJ)
	{
		cout << "  BIONJ tree:                    " << params.out_prefix << ".bionj"
			<< endl;
	}
	if (!params.dist_file)
	{
		// cout << "  Juke-Cantor distances:    " << params.out_prefix << ".jcdist" << endl;
		if (params.compute_ml_dist)
			cout << "  Likelihood distances:          " << params.out_prefix
			<< ".mldist" << endl;
		if (params.print_conaln)
			cout << "  Concatenated alignment:        " << params.out_prefix
			<< ".conaln" << endl;
	}

	if (tree.getRate()->getGammaShape() > 0 && params.print_site_rate)
		cout << "  Gamma-distributed rates:       " << params.out_prefix << ".rate"
		<< endl;

	if ((tree.getRate()->isSiteSpecificRate() || tree.getRate()->getPtnCat(0) >= 0) && params.print_site_rate)
		cout << "  Site-rates by MH model:        " << params.out_prefix << ".rate"
		<< endl;

	if (params.print_site_lh)
		cout << "  Site log-likelihoods:          " << params.out_prefix << ".sitelh"
		<< endl;

	if (params.print_site_pars)
		cout << "  Site parsimony scores:         " << params.out_prefix << ".sitepars"
		<< endl;

	if (params.write_intermediate_trees)
		cout << "  All intermediate trees:        " << params.out_prefix << ".treels"
		<< endl;

	if (params.gbo_replicates)
	{
		if (params.maximum_parsimony)
			cout << endl
			<< "Maximum parsimony bootstrap approximation results written to:" << endl
			<< "  Split support values:          " << params.out_prefix << ".splits.nex" << endl
			<< "  Consensus tree:                " << params.out_prefix << ".contree" << endl;
		else
			cout << endl
			<< "Ultrafast bootstrap approximation results written to:" << endl
			<< "  Split support values:          " << params.out_prefix << ".splits.nex" << endl
			<< "  Consensus tree:                " << params.out_prefix << ".contree" << endl;
		if (params.print_ufboot_trees)
			cout << "  UFBoot trees:                  " << params.out_prefix << ".ufboot" << endl;
	}

	if (params.treeset_file)
	{
		cout << "  Evaluated user trees:          " << params.out_prefix << ".trees" << endl;

		if (params.print_tree_lh)
		{
			cout << "  Tree log-likelihoods:          " << params.out_prefix << ".treelh" << endl;
		}
		if (params.print_site_lh)
		{
			cout << "  Site log-likelihoods:          " << params.out_prefix << ".sitelh" << endl;
		}
		if (params.print_site_pars)
		{
			cout << "  Site parsimony scores:         " << params.out_prefix << ".sitepars" << endl;
		}
	}
	cout << "  Screen log file:               " << params.out_prefix << ".log" << endl;
	/*	if (original_model == "WHTEST")
	 cout <<"  WH-TEST report:           " << params.out_prefix << ".whtest" << endl;*/
	cout << endl;
}

void checkZeroDist(Alignment* aln, double* dist)
{
	int ntaxa = aln->getNSeq();
	IntVector checked;
	checked.resize(ntaxa, 0);
	int i, j;
	for (i = 0; i < ntaxa - 1; i++)
	{
		if (checked[i])
			continue;
		string str = "";
		bool first = true;
		for (j = i + 1; j < ntaxa; j++)
			if (dist[i * ntaxa + j] <= 1e-6)
			{
				if (first)
					str = "ZERO distance between sequences " + aln->getSeqName(i);
				str += ", " + aln->getSeqName(j);
				checked[j] = 1;
				first = false;
			}
		checked[i] = 1;
		if (str != "")
			outWarning(str);
	}
}

void printAnalysisInfo(int model_df, IQTree& iqtree, Params& params)
{
	//	if (!params.raxmllib) {
	cout << "Model of evolution: ";
	if (iqtree.isSuperTree())
	{
		cout << iqtree.getModelName() << " (" << model_df << " free parameters)" << endl;
	}
	else
	{
		cout << iqtree.getModelName() << " with ";
		switch (iqtree.getModel()->getFreqType())
		{
		case FREQ_EQUAL:
			cout << "equal";
			break;
		case FREQ_EMPIRICAL:
			cout << "counted";
			break;
		case FREQ_USER_DEFINED:
			cout << "user-defined";
			break;
		case FREQ_ESTIMATE:
			cout << "optimized";
			break;
		case FREQ_CODON_1x4:
			cout << "counted 1x4";
			break;
		case FREQ_CODON_3x4:
			cout << "counted 3x4";
			break;
		case FREQ_CODON_3x4C:
			cout << "counted 3x4-corrected";
			break;
		default:
			outError("Wrong specified state frequencies");
		}
		cout << " frequencies (" << model_df << " free parameters)" << endl;
	}
	cout << "Fixed branch lengths: "
		<< ((params.fixed_branch_length) ? "Yes" : "No") << endl;

	if (params.min_iterations > 0)
	{
		cout << "Tree search algorithm: " << (params.snni ? "Stochastic nearest neighbor interchange" : "IQPNNI") << endl;
		cout << "Termination condition: ";
		if (params.stop_condition == SC_REAL_TIME)
		{
			cout << "after " << params.maxtime << " minutes" << endl;
		}
		else if (params.stop_condition == SC_UNSUCCESS_ITERATION)
		{
			cout << "after " << params.unsuccess_iteration << " unsuccessful iterations" << endl;
		}
		else if (params.stop_condition == SC_FIXED_ITERATION)
		{
			cout << params.min_iterations << " iterations" << endl;
		}
		else if (params.stop_condition == SC_WEIBULL)
		{
			cout << "predicted in [" << params.min_iterations << ","
				<< params.max_iterations << "] (confidence "
				<< params.stop_confidence << ")" << endl;
		}
		else if (params.stop_condition == SC_BOOTSTRAP_CORRELATION)
		{
			cout << "min " << params.min_correlation << " correlation coefficient" << endl;
		}

		if (!params.snni)
		{
			cout << "Number of representative leaves  : " << params.k_representative << endl;
			cout << "Probability of deleting sequences: " << iqtree.getProbDelete() << endl;
			cout << "Number of leaves to be deleted   : " << iqtree.getDelete() << endl;
			cout << "Important quartets assessed on: "
				<< ((params.iqp_assess_quartet == IQP_DISTANCE) ? "Distance" : ((params.iqp_assess_quartet == IQP_PARSIMONY) ? "Parsimony" : "Bootstrap"))
				<< endl;
		}
		cout << "NNI assessed on: " << ((params.nni5) ? "5 branches" : "1 branch") << endl;
	}
	cout << "Phylogenetic likelihood library: " << (params.pll ? "Yes" : "No") << endl;
	cout << "Branch length optimization method: "
		<< ((iqtree.optimize_by_newton) ? "Newton" : "Brent") << endl;
	cout << "Number of Newton-Raphson steps in NNI evaluation and branch length optimization: " << NNI_MAX_NR_STEP
		<< " / " << PLL_NEWZPERCYCLE << endl;
	cout << "SSE instructions: "
		<< ((iqtree.sse) ? "Yes" : "No") << endl;
	if (params.maximum_parsimony && params.spr_parsimony)
	{
		cout << "Hill-climbing using SPR radius: " << params.spr_maxtrav << endl;
	}
	cout << endl;
}

void computeMLDist(Params& params, IQTree& iqtree, string& dist_file, double begin_time, double& bestTreeScore)
{
	double longest_dist;
	stringstream best_tree_string;
	iqtree.printTree(best_tree_string, WT_BR_LEN + WT_TAXON_ID);
	cout << "Computing ML distances based on estimated model parameters...";
	double* ml_dist = NULL;
	double* ml_var = NULL;
	longest_dist = iqtree.computeDist(params, iqtree.aln, ml_dist, ml_var, dist_file);
	cout << " " << (getCPUTime() - begin_time) << " sec" << endl;
	cout << endl;
	if (longest_dist > MAX_GENETIC_DIST * 0.99)
	{
		outWarning("Some pairwise ML distances are too long (saturated)");
		// cout << "Some ML distances are too long, using old distances..." << endl;
	} // else
	{
		if (!iqtree.dist_matrix)
		{
			iqtree.dist_matrix = new double[iqtree.aln->getNSeq() * iqtree.aln->getNSeq()];
		}
		if (!iqtree.var_matrix)
		{
			iqtree.var_matrix = new double[iqtree.aln->getNSeq() * iqtree.aln->getNSeq()];
		}
		memmove(iqtree.dist_matrix, ml_dist,
			sizeof(double) * iqtree.aln->getNSeq() * iqtree.aln->getNSeq());
		memmove(iqtree.var_matrix, ml_var,
			sizeof(double) * iqtree.aln->getNSeq() * iqtree.aln->getNSeq());
	}
	delete[] ml_dist;
	delete[] ml_var;
}

void computeInitialDist(Params& params, IQTree& iqtree, string& dist_file)
{
	double longest_dist;
	if (params.dist_file)
	{
		cout << "Reading distance matrix file " << params.dist_file << " ..." << endl;
	}
	else if (params.compute_jc_dist)
	{
		cout << "Computing Juke-Cantor distances..." << endl;
	}
	else if (params.compute_obs_dist)
	{
		cout << "Computing observed distances..." << endl;
	}

	if (params.compute_jc_dist || params.compute_obs_dist || params.partition_file)
	{
		longest_dist = iqtree.computeDist(params, iqtree.aln, iqtree.dist_matrix, iqtree.var_matrix, dist_file);
		checkZeroDist(iqtree.aln, iqtree.dist_matrix);
		if (longest_dist > MAX_GENETIC_DIST * 0.99)
		{
			outWarning("Some pairwise distances are too long (saturated)");
		}
	}
}

void computeInitialTree(Params& params, IQTree& iqtree, string& dist_file, int& numInitTrees, string& initTree)
{
	double start = getCPUTime();

	string out_file = params.out_prefix;
	if (params.user_file)
	{
		// start the search with user-defined tree
		cout << endl;
		cout << "Reading input tree file " << params.user_file << " ..." << endl;
		bool myrooted = params.is_rooted;
		iqtree.readTree(params.user_file, myrooted);
		iqtree.setAlignment(iqtree.aln);
		iqtree.initializeAllPartialPars(); // 2020-08-17: Diep added to fix bug while compute score of user tree
		iqtree.clearAllPartialLH();		   // 2020-08-17: Diep added to fix bug while compute score of user tree

		numInitTrees = 1;
		params.numNNITrees = 1;
		// change to old kernel if tree is multifurcating
		if ((params.SSE == LK_EIGEN || params.SSE == LK_EIGEN_SSE) && !iqtree.isBifurcating())
		{
			cout << "NOTE: Changing to old kernel as input tree is multifurcating" << endl;
			params.SSE = LK_SSE;
		}
	}
	else
		switch (params.start_tree)
		{
		case STT_PARSIMONY:
			// Create parsimony tree using IQ-Tree kernel
			cout << endl;
			cout << "Creating initial parsimony tree by random order stepwise addition..." << endl;
			iqtree.computeParsimonyTree(params.out_prefix, iqtree.aln);
			iqtree.initializeAllPartialPars();
			iqtree.clearAllPartialLH();
			iqtree.fixNegativeBranch(true);
			numInitTrees = params.numParsTrees;
			if (numInitTrees > params.min_iterations && params.stop_condition == SC_FIXED_ITERATION)
				numInitTrees = params.min_iterations;
			break;
		case STT_PLL_PARSIMONY:
			cout << endl;
			cout << "Create initial parsimony tree by phylogenetic likelihood library (PLL)... ";
			// generate a parsimony tree for model optimization
			iqtree.pllInst->randomNumberSeed = params.ran_seed;

			if (params.maximum_parsimony)
			{
				_pllComputeRandomizedStepwiseAdditionParsimonyTree(iqtree.pllInst, iqtree.pllPartitions, params.sprDist, &iqtree);
			}
			else
				pllComputeRandomizedStepwiseAdditionParsimonyTree(iqtree.pllInst, iqtree.pllPartitions, params.sprDist);

			resetBranches(iqtree.pllInst);
			pllTreeToNewick(iqtree.pllInst->tree_string, iqtree.pllInst, iqtree.pllPartitions, iqtree.pllInst->start->back,
				PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
			iqtree.readTreeString(string(iqtree.pllInst->tree_string));
			iqtree.initializeAllPartialPars();
			iqtree.clearAllPartialLH();
			if (params.write_init_tree)
			{
				out_file += ".parstree";
				iqtree.printTree(out_file.c_str(), WT_NEWLINE);
			}
			iqtree.fixNegativeBranch(true);

			cout << getCPUTime() - start << " seconds" << endl;
			numInitTrees = params.numParsTrees;
			if (numInitTrees > params.min_iterations && params.stop_condition == SC_FIXED_ITERATION)
				numInitTrees = params.min_iterations;
			break;
		case STT_BIONJ:
			// This is the old default option: using BIONJ as starting tree
			iqtree.computeBioNJ(params, iqtree.aln, dist_file);
			cout << getCPUTime() - start << " seconds" << endl;
			numInitTrees = 1;
			break;
		}

	/* Fix if negative branch lengths detected */
	int fixed_number = iqtree.fixNegativeBranch();
	if (fixed_number)
	{
		cout << "WARNING: " << fixed_number << " undefined/negative branch lengths are initialized with parsimony" << endl;
	}
	if (params.root)
	{
		string str = params.root;
		if (!iqtree.findNodeName(str))
		{
			str = "Specified root name " + str + "not found";
			outError(str);
		}
	}
	initTree = iqtree.getTreeString();
	if (params.pll)
	{
		pllNewickTree* newick = pllNewickParseString(initTree.c_str());
		pllTreeInitTopologyNewick(iqtree.pllInst, newick, PLL_TRUE);
		pllNewickParseDestroy(&newick);
		//        pllInitModel(iqtree.pllInst, iqtree.pllPartitions, iqtree.pllAlignment);
		pllInitModel(iqtree.pllInst, iqtree.pllPartitions);
	}
}

void initializeParams(Params& params, IQTree& iqtree, vector<ModelInfo>& model_info)
{
	iqtree.curScore = -DBL_MAX;
	bool test_only = params.model_name.substr(0, 8) == "TESTONLY";
	/* initialize substitution model */
	if (params.model_name.substr(0, 4) == "TEST")
	{
		if (iqtree.isSuperTree())
			((PhyloSuperTree*)&iqtree)->mapTrees();
		uint64_t mem_size = iqtree.getMemoryRequired();
		mem_size *= (params.num_rate_cats + 1);
		cout << "NOTE: MODEL SELECTION REQUIRES AT LEAST " << ((double)mem_size * sizeof(double) / 1024.0) / 1024
			<< " MB MEMORY!" << endl;
		if (mem_size >= getMemorySize())
		{
			outError("Memory required exceeds your computer RAM size!");
		}
		params.model_name = testModel(params, &iqtree, model_info);
		cout << "CPU time for model selection: " << getCPUTime() - params.startCPUTime << " seconds." << endl;
		//        alignment = iqtree.aln;
		if (test_only)
		{
			params.min_iterations = 0;
		}
	}

	if (params.model_name == "WHTEST")
	{
		if (iqtree.aln->seq_type != SEQ_DNA)
			outError("Weiss & von Haeseler test of model homogeneity only works for DNA");
		params.model_name = "GTR+G";
	}

	assert(iqtree.aln);
	if (params.gbo_replicates)
		params.speed_conf = 1.0;

	if (iqtree.isSuperTree())
		((PhyloSuperTree*)&iqtree)->mapTrees();

	// set parameter for the current tree
	iqtree.setParams(params);
}
/*
 *  Generate the initial candidate tree set
 *  @param numInitTrees number of parsimony trees to use
 *  @return number of duplicated trees
 */
int initCandidateTreeSet(Params& params, IQTree& iqtree, int numInitTrees)
{
	int nni_count = 0;
	int nni_steps = 0;
	int numDup = 0;
	cout << "Generating " << numInitTrees - 1 << " parsimony trees... ";
	cout.flush();
	double startTime = getCPUTime();
	int numDupPars = 0;
	//    if(params.maximum_parsimony) iqtree.candidateTrees.clear(); // Diep: added this to fix the bug of sorted aln <> orig aln
	for (int treeNr = 1; treeNr < numInitTrees; treeNr++)
	{
		string curParsTree;
		if (params.start_tree == STT_PLL_PARSIMONY)
		{
			iqtree.pllInst->randomNumberSeed = params.ran_seed + treeNr * 12345;

			if (params.maximum_parsimony)
			{
				_pllComputeRandomizedStepwiseAdditionParsimonyTree(iqtree.pllInst, iqtree.pllPartitions, params.sprDist, &iqtree);
			}
			else
				pllComputeRandomizedStepwiseAdditionParsimonyTree(iqtree.pllInst, iqtree.pllPartitions, params.sprDist);

			pllTreeToNewick(iqtree.pllInst->tree_string, iqtree.pllInst, iqtree.pllPartitions,
				iqtree.pllInst->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
				PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
			curParsTree = string(iqtree.pllInst->tree_string);
		}
		else
		{
			iqtree.computeParsimonyTree(NULL, iqtree.aln);
			curParsTree = iqtree.getTreeString();
		}
		if (iqtree.candidateTrees.treeExist(curParsTree))
		{
			numDupPars++;
			continue;
		}
		else
		{
			if (params.start_tree == STT_PLL_PARSIMONY)
				iqtree.readTreeString(curParsTree);
			if (params.count_trees)
			{
				string tree = iqtree.getTopology();
				if (pllTreeCounter.find(tree) == pllTreeCounter.end())
				{
					// not found in hash_map
					pllTreeCounter[curParsTree] = 1;
				}
				else
				{
					// found in hash_map
					pllTreeCounter[curParsTree]++;
				}
			}
			// Diep added IF statement for MP doesn't need branch optimization
			if (params.maximum_parsimony)
			{
				iqtree.initializeAllPartialPars();
				iqtree.clearAllPartialLH();
				iqtree.curScore = -iqtree.computeParsimony();
				iqtree.candidateTrees.update(curParsTree, iqtree.curScore);
				if (iqtree.curScore > iqtree.bestScore)
				{
					iqtree.setBestTree(curParsTree, iqtree.curScore);
				}
			}
			else
				iqtree.candidateTrees.update(curParsTree, -DBL_MAX);
		}
	}
	double parsTime = getCPUTime() - startTime;
	cout << "(" << numDupPars << " duplicated parsimony trees)" << endl;
	cout << "CPU time: " << parsTime << endl;

	// do not do anything for parsimony because tree was already optimized by SPR
	if (params.maximum_parsimony)
	{
		return -1;
	}

	cout << "Computing log-likelihood of the parsimony trees ... " << endl;

	startTime = getCPUTime();
	vector<string> unOptParTrees = iqtree.candidateTrees.getHighestScoringTrees(numInitTrees);
	for (vector<string>::iterator it = unOptParTrees.begin() + 1; it != unOptParTrees.end(); it++)
	{
		iqtree.readTreeString(*it);
		// Initialize branch lengths for the parsimony tree
		iqtree.initializeAllPartialPars();
		iqtree.clearAllPartialLH();
		iqtree.fixNegativeBranch(true);
		//        if (iqtree.isSuperTree()) {
		//            iqtree.assignRandomBranchLengths(true);
		//            ((PhyloSuperTree*)&iqtree)->mapTrees();
		//        } else {
		//        	iqtree.fixNegativeBranch(true);
		//    	}
		iqtree.initializeAllPartialLh();
		iqtree.clearAllPartialLH();
		// Optimize the branch lengths
		string tree = iqtree.optimizeBranches(2);
		// Add tree to the candidate set
		iqtree.candidateTrees.update(tree, iqtree.curScore);
		if (iqtree.curScore > iqtree.bestScore)
		{
			iqtree.setBestTree(tree, iqtree.curScore);
		}
	}
	double loglTime = getCPUTime() - startTime;
	cout << "CPU time: " << loglTime << endl;

	CandidateSet initParsimonyTrees = iqtree.candidateTrees.getBestCandidateTrees(params.numNNITrees);
	iqtree.candidateTrees.clear();
	if (verbose_mode >= VB_MED)
	{
		for (multimap<double, CandidateTree>::iterator it = iqtree.candidateTrees.begin(); it != iqtree.candidateTrees.end(); it++)
		{
			cout << it->first << " / " << it->second.topology << endl;
		}
	}

	// iqtree.candidateTrees.clear();

	/************ END: Create a set of up to (numInitTrees - 1) unique parsimony trees **********************/

	cout << endl;
	cout << "Optimizing top " << initParsimonyTrees.size() << " parsimony trees with NNI..." << endl;
	startTime = getCPUTime();
	/*********** START: Do NNI on the best parsimony trees ************************************/
	CandidateSet::reverse_iterator rit;
	iqtree.setCurIt(1);
	for (rit = initParsimonyTrees.rbegin(); rit != initParsimonyTrees.rend(); ++rit, iqtree.setCurIt(iqtree.getCurIt() + 1))
	{
		int nniCount, nniStep;
		double initLogl, nniLogl;
		string tree;
		iqtree.readTreeString(rit->second.tree);
		// cout << rit->second.tree << endl;
		iqtree.initializeAllPartialLh();
		iqtree.clearAllPartialLH();
		iqtree.computeLogL();
		// THIS HAPPEN WHENEVER USING FULL PARTITION MODEL
		//        while (iqtree.curScore - rit->first < -1.0) {
		if (iqtree.isSuperTree() && params.partition_type == 0)
		{
			if (verbose_mode >= VB_MED)
				cout << "curScore: " << iqtree.curScore << " expected score: " << rit->first << endl;
			iqtree.optimizeBranches(2);
		}
		initLogl = iqtree.curScore;
		tree = iqtree.doNNISearch(nniCount, nniStep);
		nniLogl = iqtree.curScore;
		cout << "Iteration " << iqtree.getCurIt() << " / LogL: " << iqtree.curScore;
		if (verbose_mode >= VB_MED)
		{
			cout << " / NNI count, steps: " << nniCount << "," << nniStep;
			cout << " / Parsimony logl " << initLogl << " / NNI logl: " << nniLogl;
		}
		cout << " / Time: " << convert_time(getRealTime() - params.start_real_time) << endl;

		bool newTree = iqtree.candidateTrees.update(tree, iqtree.curScore);
		if (!newTree)
		{
			numDup++;
		}
		// Better tree is found
		if (iqtree.curScore > iqtree.bestScore && newTree)
		{
			// Re-optimize model parameters (the sNNI algorithm)
			tree = iqtree.optimizeModelParameters();
			iqtree.setBestTree(tree, iqtree.curScore);
			cout << "BETTER TREE FOUND: " << iqtree.bestScore << endl;
		}
	}
	double nniTime = getCPUTime() - startTime;
	cout << "Average time for 1 NNI search: " << nniTime / initParsimonyTrees.size() << endl;
	return numDup;
}

void pruneTaxa(Params& params, IQTree& iqtree, double* pattern_lh, NodeVector& pruned_taxa, StrVector& linked_name)
{
	int num_low_support;
	double mytime;

	if (params.aLRT_threshold <= 100 && (params.aLRT_replicates > 0 || params.localbp_replicates > 0))
	{
		mytime = getCPUTime();
		cout << "Testing tree branches by SH-like aLRT with " << params.aLRT_replicates << " replicates..." << endl;
		iqtree.setRootNode(params.root);
		iqtree.computePatternLikelihood(pattern_lh, &iqtree.curScore);
		num_low_support = iqtree.testAllBranches(params.aLRT_threshold, iqtree.curScore,
			pattern_lh, params.aLRT_replicates, params.localbp_replicates);
		iqtree.printResultTree();
		cout << "  " << getCPUTime() - mytime << " sec." << endl;
		cout << num_low_support << " branches show low support values (<= " << params.aLRT_threshold << "%)" << endl;

		// tree.drawTree(cout);
		cout << "Collapsing stable clades..." << endl;
		iqtree.collapseStableClade(params.aLRT_threshold, pruned_taxa, linked_name, iqtree.dist_matrix);
		cout << pruned_taxa.size() << " taxa were pruned from stable clades" << endl;
	}

	if (!pruned_taxa.empty())
	{
		cout << "Pruned alignment contains " << iqtree.aln->getNSeq()
			<< " sequences and " << iqtree.aln->getNSite() << " sites and "
			<< iqtree.aln->getNPattern() << " patterns" << endl;
		// tree.clearAllPartialLh();
		iqtree.initializeAllPartialLh();
		iqtree.clearAllPartialLH();
		iqtree.curScore = iqtree.optimizeAllBranches();
		// cout << "Log-likelihood	after reoptimizing model parameters: " << tree.curScore << endl;
		int nni_count, nni_steps;
		iqtree.curScore = iqtree.optimizeNNI(nni_count, nni_steps);
		cout << "Log-likelihood after optimizing partial tree: "
			<< iqtree.curScore << endl;
	}
}

void restoreTaxa(IQTree& iqtree, double* saved_dist_mat, NodeVector& pruned_taxa, StrVector& linked_name)
{
	if (!pruned_taxa.empty())
	{
		cout << "Restoring full tree..." << endl;
		iqtree.restoreStableClade(iqtree.aln, pruned_taxa, linked_name);
		delete[] iqtree.dist_matrix;
		iqtree.dist_matrix = saved_dist_mat;
		iqtree.initializeAllPartialLh();
		iqtree.clearAllPartialLH();
		iqtree.curScore = iqtree.optimizeAllBranches();
		// cout << "Log-likelihood	after reoptimizing model parameters: " << tree.curScore << endl;
		int nni_count, nni_steps;
		iqtree.curScore = iqtree.optimizeNNI(nni_count, nni_steps);
		cout << "Log-likelihood	after reoptimizing full tree: "
			<< iqtree.curScore << endl; // iqtree.setBestScore(iqtree.getModelFactory()->optimizeParameters(params.fixed_branch_length, true, params.model_eps));
	}
}
void runApproximateBranchLengths(Params& params, IQTree& iqtree)
{
	if (!params.fixed_branch_length && params.leastSquareBranch)
	{
		cout << endl
			<< "Computing Least Square branch lengths..." << endl;
		iqtree.optimizeAllBranchesLS();
		iqtree.clearAllPartialLH();
		iqtree.curScore = iqtree.computeLikelihood();
		string filename = params.out_prefix;
		filename += ".lstree";
		iqtree.printTree(filename.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
		cout << "Logl of tree with LS branch lengths: " << iqtree.curScore << endl;
		cout << "Tree with LS branch lengths written to " << filename << endl;
		if (params.print_branch_lengths)
		{
			if (params.manuel_analytic_approx)
			{
				cout << "Applying Manuel's analytic approximation.." << endl;
				iqtree.approxAllBranches();
			}
			ofstream out;
			filename = params.out_prefix;
			filename += ".lsbrlen";
			out.open(filename.c_str());
			iqtree.printBranchLengths(out);
			out.close();
			cout << "LS Branch lengths written to " << filename << endl;
		}
		cout << "Total LS tree length: " << iqtree.treeLength() << endl;
	}

	if (params.pars_branch_length)
	{
		cout << endl
			<< "Computing parsimony branch lengths..." << endl;
		iqtree.fixNegativeBranch(true);
		iqtree.clearAllPartialLH();
		iqtree.curScore = iqtree.computeLikelihood();
		string filename = params.out_prefix;
		filename += ".mptree";
		iqtree.printTree(filename.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
		cout << "Logl of tree with MP branch lengths: " << iqtree.curScore << endl;
		cout << "Tree with MP branch lengths written to " << filename << endl;
		if (params.print_branch_lengths)
		{
			ofstream out;
			filename = params.out_prefix;
			filename += ".mpbrlen";
			out.open(filename.c_str());
			iqtree.printBranchLengths(out);
			out.close();
			cout << "MP Branch lengths written to " << filename << endl;
		}
		cout << "Total MP tree length: " << iqtree.treeLength() << endl;
	}

	if (params.bayes_branch_length)
	{
		cout << endl
			<< "Computing Bayesian branch lengths..." << endl;
		iqtree.computeAllBayesianBranchLengths();
		iqtree.clearAllPartialLH();
		iqtree.curScore = iqtree.computeLikelihood();
		string filename = params.out_prefix;
		filename += ".batree";
		iqtree.printTree(filename.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
		cout << "Logl of tree with Bayesian branch lengths: " << iqtree.curScore << endl;
		cout << "Tree with Bayesian branch lengths written to " << filename << endl;
		if (params.print_branch_lengths)
		{
			ofstream out;
			filename = params.out_prefix;
			filename += ".babrlen";
			out.open(filename.c_str());
			iqtree.printBranchLengths(out);
			out.close();
			cout << "Bayesian Branch lengths written to " << filename << endl;
		}
		cout << "Total Bayesian tree length: " << iqtree.treeLength() << endl;
	}
}

// 2020-08-17: Diep modified this function to allow print site parsimony scores
void printMiscInfo(Params& params, IQTree& iqtree, double* pattern_lh)
{
	string ext = "";
	if (params.print_site_lh)
		ext = ".sitelh";
	else if (params.maximum_parsimony && params.print_site_pars)
		ext = ".sitepars";

	if ((params.print_site_lh || params.print_site_pars) && !params.pll)
	{
		string site_lh_file = params.out_prefix;
		site_lh_file += ext;
		if (params.print_site_lh == 1 || params.print_site_pars == 1)
			printSiteLh(site_lh_file.c_str(), &iqtree, pattern_lh);
		else
			printSiteLhCategory(site_lh_file.c_str(), &iqtree);
	}

	if (params.print_branch_lengths)
	{
		if (params.manuel_analytic_approx)
		{
			cout << "Applying Manuel's analytic approximation.." << endl;
			iqtree.approxAllBranches();
		}
		string brlen_file = params.out_prefix;
		brlen_file += ".brlen";
		ofstream out;
		out.open(brlen_file.c_str());
		iqtree.printBranchLengths(out);
		out.close();
		cout << "Branch lengths written to " << brlen_file << endl;
	}

	if (params.print_partition_info && iqtree.isSuperTree())
	{
		string partition_info = params.out_prefix;
		partition_info += ".partinfo.nex";
		((PhyloSuperTree*)(&iqtree))->printPartition(partition_info.c_str());
	}

	if (params.mvh_site_rate)
	{
		RateMeyerHaeseler* rate_mvh = new RateMeyerHaeseler(params.rate_file,
			&iqtree, params.rate_mh_type);
		cout << endl
			<< "Computing site-specific rates by "
			<< rate_mvh->full_name << "..." << endl;
		rate_mvh->runIterativeProc(params, iqtree);
		cout << endl
			<< "BEST SCORE FOUND : " << (params.maximum_parsimony ? -iqtree.getBestScore() : iqtree.getBestScore()) << endl;
		string mhrate_file = params.out_prefix;
		mhrate_file += ".mhrate";
		iqtree.getRate()->writeSiteRates(mhrate_file.c_str());

		if (params.print_site_lh)
		{
			string site_lh_file = params.out_prefix;
			site_lh_file += ".mhsitelh";
			printSiteLh(site_lh_file.c_str(), &iqtree);
		}
	}

	if (params.print_site_rate)
	{
		string rate_file = params.out_prefix;
		rate_file += ".rate";
		iqtree.getRate()->writeSiteRates(rate_file.c_str());
		if (iqtree.isSuperTree())
		{
			PhyloSuperTree* stree = (PhyloSuperTree*)&iqtree;
			int part = 0;
			try
			{
				ofstream out;
				out.exceptions(ios::failbit | ios::badbit);
				out.open(rate_file.c_str());
				for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++, part++)
				{
					out << "SITE RATES FOR PARTITION " << stree->part_info[part].name << ":" << endl;
					(*it)->getRate()->writeSiteRates(out);
				}
				cout << "Site rates printed to " << rate_file << endl;
				out.close();
			}
			catch (ios::failure)
			{
				outError(ERR_WRITE_OUTPUT, rate_file);
			}
		}
	}
}

void printFinalSearchInfo(Params& params, IQTree& iqtree, double search_cpu_time, double search_real_time)
{
	if (!params.maximum_parsimony)
		cout << "Total tree length: " << iqtree.treeLength() << endl;

	if (iqtree.isSuperTree())
	{
		PhyloSuperTree* stree = (PhyloSuperTree*)&iqtree;
		cout << stree->evalNNIs << " NNIs evaluated from " << stree->totalNNIs << " all possible NNIs ( " << (int)(((stree->evalNNIs + 1.0) / (stree->totalNNIs + 1.0)) * 100.0) << " %)" << endl;
		cout << "Details for subtrees:" << endl;
		for (int part = 0; part < stree->size(); part++)
		{
			cout << part + 1 << ". " << stree->part_info[part].name << ": " << stree->part_info[part].evalNNIs << " ( "
				<< (int)(((stree->part_info[part].evalNNIs + 1.0) / ((stree->totalNNIs + 1.0) / stree->size())) * 100.0)
				<< " %)" << endl;
		}
	}

	params.run_time = (getCPUTime() - params.startCPUTime);
	cout.precision(3);
	cout << endl;
	cout << "CPU time used for tree search: " << search_cpu_time
		<< " sec (" << convert_time(search_cpu_time) << ")" << endl;
	cout << "Wall-clock time used for tree search: " << search_real_time
		<< " sec (" << convert_time(search_real_time) << ")" << endl;
	cout << "Total CPU time used: " << (double)params.run_time << " sec ("
		<< convert_time((double)params.run_time) << ")" << endl;
	cout << "Total wall-clock time used: "
		<< getRealTime() - params.start_real_time << " sec ("
		<< convert_time(getRealTime() - params.start_real_time) << ")" << endl;
	cout.precision(0);
}

/************************************************************
 *  MAIN TREE RECONSTRUCTION
 ***********************************************************/
void runTreeReconstruction(Params& params, string& original_model, IQTree& iqtree, vector<ModelInfo>& model_info)
{
	string dist_file;
	params.startCPUTime = getCPUTime();
	params.start_real_time = getRealTime();

	// Make sure that no partial likelihood of IQ-TREE is initialized when PLL is used to save memory
	if (params.pll)
	{
		iqtree.deleteAllPartialLh();
	}

	// Temporary fix since PLL only supports DNA/Protein: switch to IQ-TREE parsimony kernel
	if (params.start_tree == STT_PLL_PARSIMONY)
	{
		if (iqtree.isSuperTree())
		{
			PhyloSuperTree* stree = (PhyloSuperTree*)&iqtree;
			for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++)
				if ((*it)->aln->seq_type != SEQ_DNA && (*it)->aln->seq_type != SEQ_PROTEIN)
					params.start_tree = STT_PARSIMONY;
		}
		else if (iqtree.aln->seq_type != SEQ_DNA && iqtree.aln->seq_type != SEQ_PROTEIN && iqtree.aln->seq_type != SEQ_BINARY)
			params.start_tree = STT_PARSIMONY;
	}
	/***************** Initialization for PLL and sNNI ******************/
	if (params.start_tree == STT_PLL_PARSIMONY || params.pll)
	{
		/* Initialized all data structure for PLL*/
		//        cout << "WHAT'S GOING ON HERE?" << endl;
		//        verbose_mode = VB_MAX;
		iqtree.initializePLL(params);
	}

	/********************* Compute pairwise distances *******************/
	if (params.start_tree == STT_BIONJ || params.iqp || params.leastSquareBranch)
	{
		computeInitialDist(params, iqtree, dist_file);
	}

	/********************** CREATE INITIAL TREE(S) **********************/
	int numInitTrees;
	string initTree;
	computeInitialTree(params, iqtree, dist_file, numInitTrees, initTree);
	/*************** SET UP PARAMETERS and model testing ****************/

	initializeParams(params, iqtree, model_info);

	/*********************** INITIAL MODEL OPTIMIZATION *****************/

	iqtree.initializeModel(params);

	// degree of freedom
	cout << endl;
	if (verbose_mode >= VB_MED)
	{
		cout << "ML-TREE SEARCH START WITH THE FOLLOWING PARAMETERS:" << endl;
		int model_df = iqtree.getModelFactory()->getNParameters();
		printAnalysisInfo(model_df, iqtree, params);
	}

	if (!params.pll)
	{
		uint64_t mem_size = iqtree.getMemoryRequired();
#if defined __APPLE__ || defined __MACH__
		cout << "NOTE: " << ((double)mem_size * sizeof(double) / 1024.0) / 1024 << " MB RAM is required!" << endl;
#else
		cout << "NOTE: " << ((double)mem_size * sizeof(double) / 1000.0) / 1000 << " MB RAM is required!" << endl;
#endif
		if (mem_size >= getMemorySize())
		{
			outError("Memory required exceeds your computer RAM size!");
		}
	}

	// Optimize model parameters and branch lengths using ML for the initial tree
	if (!params.maximum_parsimony)
		initTree = iqtree.optimizeModelParameters(true);
	else
	{
		iqtree.readTreeString(initTree);
		iqtree.initializeAllPartialPars();
		iqtree.clearAllPartialLH();
		iqtree.curScore = -iqtree.computeParsimony();
	}

	/****************** NOW PERFORM MAXIMUM LIKELIHOOD TREE RECONSTRUCTION ******************/
	cout.precision(0);
	// Update best tree
	iqtree.candidateTrees.clear(); // Diep added
	iqtree.setBestTree(initTree, iqtree.curScore);
	cout << "Current best tree score: " << (params.maximum_parsimony ? -iqtree.bestScore : iqtree.bestScore) << endl
		<< endl;
	iqtree.candidateTrees.update(initTree, iqtree.curScore);

	// Compute maximum likelihood distance
	// ML distance is only needed for IQP
	if ((params.snni && !params.iqp) || params.min_iterations == 0)
	{
		params.compute_ml_dist = false;
	}
	if ((!params.dist_file && params.compute_ml_dist) || params.leastSquareBranch)
	{
		computeMLDist(params, iqtree, dist_file, getCPUTime(), iqtree.bestScore);
	}

	double cputime_search_start = getCPUTime();
	double realtime_search_start = getRealTime();

	if (params.min_iterations > 0)
	{
		double initTime = getCPUTime();

		if (!params.user_file && (params.start_tree == STT_PARSIMONY || params.start_tree == STT_PLL_PARSIMONY))
		{
			int numDup = initCandidateTreeSet(params, iqtree, numInitTrees);
			assert(iqtree.candidateTrees.size() != 0);
			cout << "Finish initializing candidate tree set. ";
			cout << "Number of distinct locally optimal trees: " << iqtree.candidateTrees.size() << endl;
		}
		else
		{
			int nni_count = 0;
			int nni_steps = 0;
			cout << "Doing NNI on the initial tree ... " << endl;
			string tree = iqtree.doNNISearch(nni_count, nni_steps);
			//            if (params.pll) {
			//                iqtree.curScore = iqtree.pllOptimizeNNI(nni_count, nni_steps, iqtree.searchinfo);
			//                pllTreeToNewick(iqtree.pllInst->tree_string, iqtree.pllInst, iqtree.pllPartitions,
			//                        iqtree.pllInst->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
			//                        PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
			//                iqtree.setBestTree(string(iqtree.pllInst->tree_string), iqtree.curScore);
			//            } else {
			//                iqtree.curScore = iqtree.optimizeNNI(nni_count, nni_steps);
			//                iqtree.setBestTree(iqtree.getTreeString(), iqtree.curScore);
			//            }
			//
			//        	if (iqtree.isSuperTree())
			//        		((PhyloSuperTree*) &iqtree)->computeBranchLengths();
			iqtree.setBestTree(tree, iqtree.curScore);
			iqtree.candidateTrees.update(tree, iqtree.curScore);
		}

		cout << "Current best score: " << (params.maximum_parsimony ? -iqtree.bestScore : iqtree.bestScore) << " / CPU time: "
			<< getCPUTime() - initTime << endl
			<< endl;
	}

	if (params.leastSquareNNI)
	{
		iqtree.computeSubtreeDists();
	}
	iqtree.setRootNode(params.root); // Important for NNI below

	if (original_model == "WHTEST")
	{
		cout << endl
			<< "Testing model homogeneity by Weiss & von Haeseler (2003)..." << endl;
		WHTest(params, iqtree);
	}

	NodeVector pruned_taxa;
	StrVector linked_name;
	double* saved_dist_mat = iqtree.dist_matrix;
	double* pattern_lh = NULL;

	if (!params.maximum_parsimony)
	{ // MP doesn't need this
		pattern_lh = new double[iqtree.getAlnNPattern()];

		// prune stable taxa
		pruneTaxa(params, iqtree, pattern_lh, pruned_taxa, linked_name);
	}

	/****************** Do tree search ***************************/
	if (params.min_iterations > 1)
	{
		iqtree.readTreeString(iqtree.bestTreeString);
		iqtree.doTreeSearch();
		iqtree.setAlignment(iqtree.aln);
	}
	else
	{
		/* do SPR with likelihood function */
		if (params.tree_spr)
		{
			// tree.optimizeSPRBranches();
			cout << "Doing SPR Search" << endl;
			cout << "Start tree.optimizeSPR()" << endl;
			double spr_score = iqtree.optimizeSPR();
			cout << "Finish tree.optimizeSPR()" << endl;
			// double spr_score = tree.optimizeSPR(tree.curScore, (PhyloNode*) tree.root->neighbors[0]->node);
			if (spr_score <= iqtree.curScore)
			{
				cout << "SPR search did not found any better tree" << endl;
			}
		}
	}

	// restore pruned taxa
	if (!params.maximum_parsimony)
		restoreTaxa(iqtree, saved_dist_mat, pruned_taxa, linked_name);

	double search_cpu_time = getCPUTime() - cputime_search_start;
	double search_real_time = getRealTime() - realtime_search_start;

	if (iqtree.isSuperTree())
		((PhyloSuperTree*)&iqtree)->mapTrees();
	if (params.snni && params.min_iterations)
	{
		cout << (params.maximum_parsimony ? "Scores" : "Log-likelihoods") << " of best " << params.popSize << " trees: " << endl;
		iqtree.printBestScores(iqtree.candidateTrees.popSize);
	}

	/******** Performs final model parameters optimization ******************/
	if (params.min_iterations)
	{
		iqtree.readTreeString(iqtree.bestTreeString);
		iqtree.initializeAllPartialLh();
		iqtree.clearAllPartialLH();
		if (params.maximum_parsimony)
			iqtree.computeParsimony();
		else
			iqtree.bestTreeString = iqtree.optimizeModelParameters(true);
	}
	else
	{
		iqtree.setBestScore(iqtree.curScore);
	}

	if (iqtree.isSuperTree())
		((PhyloSuperTree*)&iqtree)->computeBranchLengths();

	cout << "BEST SCORE FOUND : " << (params.maximum_parsimony ? -iqtree.getBestScore() : iqtree.getBestScore()) << endl;

	if (params.write_local_optimal_trees)
	{
		vector<string> trees = iqtree.candidateTrees.getHighestScoringTrees();
		ofstream treesOut((string(params.out_prefix) + ".localtrees").c_str(), ofstream::out);
		for (vector<string>::iterator it = trees.begin(); it != trees.end(); it++)
			treesOut << (*it) << endl;
	}

	if (params.pll)
		iqtree.inputModelPLL2IQTree();

	/* root the tree at the first sequence */
	iqtree.root = iqtree.findLeafName(iqtree.aln->getSeqName(0));
	assert(iqtree.root);

	double myscore = 0.0;

	myscore = iqtree.getBestScore();

	if (!params.pll)
	{
		iqtree.computeLikelihood(pattern_lh);
		// compute logl variance
		iqtree.logl_variance = iqtree.computeLogLVariance();
	}

	printMiscInfo(params, iqtree, pattern_lh);

	/****** perform SH-aLRT test ******************/
	if ((params.aLRT_replicates > 0 || params.localbp_replicates > 0) && !params.pll)
	{
		double mytime = getCPUTime();
		cout << endl
			<< "Testing tree branches by SH-like aLRT with "
			<< params.aLRT_replicates << " replicates..." << endl;
		iqtree.setRootNode(params.root);
		iqtree.testAllBranches(params.aLRT_threshold, myscore,
			pattern_lh, params.aLRT_replicates, params.localbp_replicates);
		cout << "CPU Time used:  " << getCPUTime() - mytime << " sec." << endl;
	}

	if (params.gbo_replicates > 0)
	{
		if (!params.online_bootstrap)
			runGuidedBootstrap(params, iqtree.aln, iqtree);
		else
			iqtree.summarizeBootstrap(params);
	}

	printFinalSearchInfo(params, iqtree, search_cpu_time, search_real_time);

	// BUG FIX: readTreeString(bestTreeString) not needed before this line
	iqtree.printResultTree();

	// To free data structures of PLL parsimony at the very end if NOT running ratchet or NOT -opt_btree
	if (params.maximum_parsimony && params.spr_parsimony && (params.ratchet_iter < 0) && (!params.optimize_boot_trees) && (params.snni || params.pll))
	{
		_pllFreeParsimonyDataStructures(iqtree.pllInst, iqtree.pllPartitions);
	}

	if (params.out_file)
		iqtree.printTree(params.out_file);

	if (pattern_lh)
		delete[] pattern_lh;

	runApproximateBranchLengths(params, iqtree);
}

/**********************************************************
 * STANDARD NON-PARAMETRIC BOOTSTRAP
 ***********************************************************/
void runStandardBootstrap(Params& params, string& original_model, Alignment* alignment, IQTree* tree)
{
	vector<ModelInfo> model_info;
	StrVector removed_seqs, twin_seqs;

	// turn off aLRT test
	int saved_aLRT_replicates = params.aLRT_replicates;
	params.aLRT_replicates = 0;
	string treefile_name = params.out_prefix;
	treefile_name += ".treefile";
	string boottrees_name = params.out_prefix;
	boottrees_name += ".boottrees";
	string bootaln_name = params.out_prefix;
	bootaln_name += ".bootaln";
	string bootlh_name = params.out_prefix;
	bootlh_name += ".bootlh";
	// first empty the boottrees file
	try
	{
		ofstream tree_out;
		tree_out.exceptions(ios::failbit | ios::badbit);
		tree_out.open(boottrees_name.c_str());
		tree_out.close();
	}
	catch (ios::failure)
	{
		outError(ERR_WRITE_OUTPUT, boottrees_name);
	}

	// empty the bootaln file
	if (params.print_bootaln)
		try
	{
		ofstream tree_out;
		tree_out.exceptions(ios::failbit | ios::badbit);
		tree_out.open(bootaln_name.c_str());
		tree_out.close();
	}
	catch (ios::failure)
	{
		outError(ERR_WRITE_OUTPUT, bootaln_name);
	}

	double start_time = getCPUTime();

	// do bootstrap analysis
	for (int sample = 0; sample < params.num_bootstrap_samples; sample++)
	{
		cout << endl
			<< "===> START BOOTSTRAP REPLICATE NUMBER "
			<< sample + 1 << endl
			<< endl;

		Alignment* bootstrap_alignment;
		cout << "Creating bootstrap alignment..." << endl;
		if (alignment->isSuperAlignment())
			bootstrap_alignment = new SuperAlignment;
		else
			bootstrap_alignment = new Alignment;
		bootstrap_alignment->createBootstrapAlignment(alignment, NULL, params.bootstrap_spec);
		if (params.print_tree_lh)
		{
			double prob;
			bootstrap_alignment->multinomialProb(*alignment, prob);
			ofstream boot_lh;
			if (sample == 0)
				boot_lh.open(bootlh_name.c_str());
			else
				boot_lh.open(bootlh_name.c_str(), ios_base::out | ios_base::app);
			boot_lh << "0\t" << prob << endl;
			boot_lh.close();
		}
		IQTree* boot_tree;
		if (alignment->isSuperAlignment())
		{
			if (params.partition_type)
			{
				boot_tree = new PhyloSuperTreePlen((SuperAlignment*)bootstrap_alignment, (PhyloSuperTree*)tree);
			}
			else
			{
				boot_tree = new PhyloSuperTree((SuperAlignment*)bootstrap_alignment, (PhyloSuperTree*)tree);
			}
		}
		else
			boot_tree = new IQTree(bootstrap_alignment);
		if (params.print_bootaln)
			bootstrap_alignment->printPhylip(bootaln_name.c_str(), true);
		runTreeReconstruction(params, original_model, *boot_tree, model_info);
		// read in the output tree file
		string tree_str;
		try
		{
			ifstream tree_in;
			tree_in.exceptions(ios::failbit | ios::badbit);
			tree_in.open(treefile_name.c_str());
			tree_in >> tree_str;
			tree_in.close();
		}
		catch (ios::failure)
		{
			outError(ERR_READ_INPUT, treefile_name);
		}
		// write the tree into .boottrees file
		try
		{
			ofstream tree_out;
			tree_out.exceptions(ios::failbit | ios::badbit);
			tree_out.open(boottrees_name.c_str(), ios_base::out | ios_base::app);
			tree_out << tree_str << endl;
			tree_out.close();
		}
		catch (ios::failure)
		{
			outError(ERR_WRITE_OUTPUT, boottrees_name);
		}
		if (params.num_bootstrap_samples == 1)
			reportPhyloAnalysis(params, original_model, *bootstrap_alignment, *boot_tree, model_info, removed_seqs, twin_seqs);
		// WHY was the following line missing, which caused memory leak?
		delete boot_tree;
		delete bootstrap_alignment;
	}

	if (params.consensus_type == CT_CONSENSUS_TREE)
	{

		cout << endl
			<< "===> COMPUTE CONSENSUS TREE FROM "
			<< params.num_bootstrap_samples << " BOOTSTRAP TREES" << endl
			<< endl;
		computeConsensusTree(boottrees_name.c_str(), 0, 1e6, -1,
			params.split_threshold, NULL, params.out_prefix, NULL, &params);
	}

	if (params.compute_ml_tree)
	{
		cout << endl
			<< "===> START ANALYSIS ON THE ORIGINAL ALIGNMENT" << endl
			<< endl;
		params.aLRT_replicates = saved_aLRT_replicates;
		runTreeReconstruction(params, original_model, *tree, model_info);

		cout << endl
			<< "===> ASSIGN BOOTSTRAP SUPPORTS TO THE TREE FROM ORIGINAL ALIGNMENT" << endl
			<< endl;
		MExtTree ext_tree;
		assignBootstrapSupport(boottrees_name.c_str(), 0, 1e6,
			treefile_name.c_str(), false, treefile_name.c_str(),
			params.out_prefix, ext_tree, NULL, &params);
		tree->copyTree(&ext_tree);
		reportPhyloAnalysis(params, original_model, *alignment, *tree, model_info, removed_seqs, twin_seqs);
	}
	else if (params.consensus_type == CT_CONSENSUS_TREE)
	{
		int mi = params.min_iterations;
		STOP_CONDITION sc = params.stop_condition;
		params.min_iterations = 0;
		params.stop_condition = SC_FIXED_ITERATION;
		runTreeReconstruction(params, original_model, *tree, model_info);
		params.min_iterations = mi;
		params.stop_condition = sc;
		tree->stop_rule.initialize(params);
		reportPhyloAnalysis(params, original_model, *alignment, *tree, model_info, removed_seqs, twin_seqs);
	}
	else
		cout << endl;

	cout << "Total CPU time for bootstrap: " << (getCPUTime() - start_time) << " seconds." << endl
		<< endl;
	cout << "Non-parametric bootstrap results written to:" << endl;
	if (params.print_bootaln)
		cout << "  Bootstrap alignments:     " << params.out_prefix << ".bootaln" << endl;
	cout << "  Bootstrap trees:          " << params.out_prefix << ".boottrees" << endl;
	if (params.consensus_type == CT_CONSENSUS_TREE)
		cout << "  Consensus tree:           " << params.out_prefix << ".contree" << endl;
	cout << endl;
}

void convertAlignment(Params& params, IQTree* iqtree)
{
	Alignment* alignment = iqtree->aln;
	if (params.num_bootstrap_samples || params.print_bootaln)
	{
		// create bootstrap alignment
		Alignment* bootstrap_alignment;
		cout << "Creating bootstrap alignment..." << endl;
		if (alignment->isSuperAlignment())
			bootstrap_alignment = new SuperAlignment;
		else
			bootstrap_alignment = new Alignment;
		bootstrap_alignment->createBootstrapAlignment(alignment, NULL, params.bootstrap_spec);
		delete alignment;
		alignment = bootstrap_alignment;
	}
	if (alignment->isSuperAlignment())
	{
		((SuperAlignment*)alignment)->printCombinedAlignment(params.aln_output);
		if (params.print_subaln)
			((SuperAlignment*)alignment)->printSubAlignments(params, ((PhyloSuperTree*)iqtree)->part_info);
	}
	else if (params.gap_masked_aln)
	{
		Alignment out_aln;
		Alignment masked_aln(params.gap_masked_aln, params.sequence_type, params.intype);
		out_aln.createGapMaskedAlignment(&masked_aln, alignment);
		out_aln.printPhylip(params.aln_output, false, params.aln_site_list,
			params.aln_nogaps, params.aln_no_const_sites, params.ref_seq_name);
		string str = params.gap_masked_aln;
		str += ".sitegaps";
		out_aln.printSiteGaps(str.c_str());
	}
	else if (params.aln_output_format == ALN_PHYLIP)
		alignment->printPhylip(params.aln_output, false, params.aln_site_list,
			params.aln_nogaps, params.aln_no_const_sites, params.ref_seq_name);
	else if (params.aln_output_format == ALN_FASTA)
		alignment->printFasta(params.aln_output, false, params.aln_site_list,
			params.aln_nogaps, params.aln_no_const_sites, params.ref_seq_name);
}

/**********************************************************
 * TOP-LEVEL FUNCTION
 ***********************************************************/
void runPhyloAnalysis(Params& params)
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

	if (params.concatenate_aln)
	{
		Alignment aln(params.concatenate_aln, params.sequence_type, params.intype);
		cout << "Concatenating " << params.aln_file << " with " << params.concatenate_aln << " ..." << endl;
		alignment->concatenateAlignment(&aln);
	}

	if (params.aln_output)
	{
		/************ convert alignment to other format and write to output file *************/
		convertAlignment(params, tree);
	}
	else if (params.gbo_replicates > 0 && params.user_file && params.second_tree)
	{
		// run one of the UFBoot analysis
		runGuidedBootstrap(params, alignment, *tree);
	}
	else if (params.avh_test)
	{
		// run one of the wondering test for Arndt
		runAvHTest(params, alignment, *tree);
	}
	else if (params.bootlh_test)
	{
		// run Arndt's plot of tree likelihoods against bootstrap alignments
		runBootLhTest(params, alignment, *tree);
	}
	else if (params.num_bootstrap_samples == 0)
	{
		// the main Maximum likelihood tree reconstruction
		vector<ModelInfo> model_info;
		alignment->checkGappySeq();

		StrVector removed_seqs;
		StrVector twin_seqs;
		// remove identical sequences
		if (params.ignore_identical_seqs)
			tree->removeIdenticalSeqs(params, removed_seqs, twin_seqs);

		// call main tree reconstruction
		runTreeReconstruction(params, original_model, *tree, model_info);

		// addMoreRowIQTree(tree, alignment);

		// addMoreRowPLL(tree, alignment, params);

		addMoreRowMutation(tree, alignment);

		if (params.gbo_replicates && params.online_bootstrap)
		{
			if (params.print_ufboot_trees)
				tree->writeUFBootTrees(params, removed_seqs, twin_seqs);

			cout << endl
				<< "Computing bootstrap consensus tree..." << endl;
			string splitsfile = params.out_prefix;
			splitsfile += ".splits.nex";
			computeConsensusTree(splitsfile.c_str(), 0, 1e6, params.split_threshold,
				params.split_weight_threshold, NULL, params.out_prefix, NULL, &params);
			// now optimize branch lengths of the consensus tree
			string current_tree = tree->getTreeString();
			splitsfile = params.out_prefix;
			splitsfile += ".contree";
			tree->readTreeFile(splitsfile);
			// bug fix
			if ((tree->sse == LK_EIGEN || tree->sse == LK_EIGEN_SSE) && !tree->isBifurcating())
			{
				cout << "NOTE: Changing to old kernel as consensus tree is multifurcating" << endl;
				tree->changeLikelihoodKernel(LK_SSE);
			}

			tree->initializeAllPartialLh();
			tree->fixNegativeBranch(true);
			//	        if (tree->isSuperTree()) {
			//	        	if (params.partition_type == 0) {
			//	        		PhyloSuperTree *stree = (PhyloSuperTree*) tree;
			//	        		tree->clearAllPartialLH();
			//	        		// full partition model
			//	        		for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++) {
			//	        			(*it)->fixNegativeBranch(true);
			//	        		}
			//	        		tree->clearAllPartialLH();
			//	        	} else {
			//	        		// joint/prop. partition model
			//					tree->assignRandomBranchLengths(true);
			//					((PhyloSuperTree*)tree)->mapTrees();
			//	        	}
			//	        } else {
			//	        	tree->fixNegativeBranch(true);
			//	    	}

			double conScore = tree->optimizeAllBranches();
			if (params.maximum_parsimony)
				cout << "Parsimony score of consensus tree: " << -conScore << endl;
			else
				cout << "Log-likelihood of consensus tree: " << conScore << endl;
			tree->setRootNode(params.root);
			tree->insertTaxa(removed_seqs, twin_seqs);
			//			tree->printTree(splitsfile.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE); // for ML
			tree->printTree(splitsfile.c_str(), WT_SORT_TAXA | WT_NEWLINE); // Diep: for MP
			// revert the best tree
			tree->readTreeString(current_tree);
			//			if (tree->isSuperTree()) {
			//				tree->optimizeAllBranches();
			//				((PhyloSuperTree*)tree)->computeBranchLengths();
			//			}
		}
		// reinsert identical sequences
		if (removed_seqs.size() > 0)
		{
			delete tree->aln;
			tree->reinsertIdenticalSeqs(alignment, removed_seqs, twin_seqs);
			tree->printResultTree();
		}
		reportPhyloAnalysis(params, original_model, *alignment, *tree, model_info, removed_seqs, twin_seqs);
	}
	else
	{
		// the classical non-parameter bootstrap (SBS)
		runStandardBootstrap(params, original_model, alignment, tree);
	}
	delete tree;
	delete alignment;
	// delete fullTree;
}

void printSiteParsimonyUserTree(Params& params)
{
	Alignment alignment(params.aln_file, params.sequence_type, params.intype);
	IQTree* ptree;

	if (params.sankoff_cost_file)
	{
		ptree = new ParsTree(&alignment);
		dynamic_cast<ParsTree*>(ptree)->initParsData(&params);
	}
	else
		ptree = new IQTree(&alignment);

	ptree->readTree(params.user_file, params.is_rooted); // Read user tree
	ptree->setAlignment(&alignment);					 // IMPORTANT: Always call setAlignment() after readTree()
	ptree->params = &params;							 // Diep: 2020-08-17, there are two variables with identical name as 'params'

	ptree->initializeAllPartialPars();
	ptree->clearAllPartialLH();
	ptree->computeParsimony();

	params.print_site_pars = 1;
	double* pattern_lh = NULL;
	printMiscInfo(params, *ptree, pattern_lh);

	delete ptree;
	if (pattern_lh)
		delete[] pattern_lh;
}

void assignBranchSupportNew(Params& params)
{
	if (!params.user_file)
		outError("No trees file provided");
	if (!params.second_tree)
		outError("No target tree file provided");
	cout << "Reading tree " << params.second_tree << " ..." << endl;
	MTree tree(params.second_tree, params.is_rooted);
	cout << tree.leafNum << " taxa and " << tree.branchNum << " branches" << endl;
	tree.assignBranchSupport(params.user_file);
	string str = params.second_tree;
	str += ".suptree";
	tree.printTree(str.c_str());
	cout << "Tree with assigned branch supports written to " << str << endl;
	if (verbose_mode >= VB_DEBUG)
		tree.drawTree(cout);
}

/**
 * assign split occurence frequencies from a set of input trees onto a target tree
 * NOTE: input trees must have the same taxon set
 * @param input_trees file containing NEWICK tree strings
 * @param burnin number of beginning trees to discard
 * @param max_count max number of trees to read in
 * @param target_tree the target tree
 * @param rooted TRUE if trees are rooted, false for unrooted trees
 * @param output_file file name to write output tree with assigned support values
 * @param out_prefix prefix of output file
 * @param mytree (OUT) resulting tree with support values assigned from target_tree
 * @param tree_weight_file file containing INTEGER weights of input trees
 * @param params program parameters
 */
void assignBootstrapSupport(const char* input_trees, int burnin, int max_count,
	const char* target_tree, bool rooted, const char* output_tree,
	const char* out_prefix, MExtTree& mytree, const char* tree_weight_file,
	Params* params)
{
	// bool rooted = false;
	//  read the tree file
	cout << "Reading tree " << target_tree << " ..." << endl;
	mytree.init(target_tree, rooted);
	// reindex the taxa in the tree to aphabetical names
	NodeVector taxa;
	mytree.getTaxa(taxa);
	sort(taxa.begin(), taxa.end(), nodenamecmp);
	int i = 0;
	for (NodeVector::iterator it = taxa.begin(); it != taxa.end(); it++)
	{
		(*it)->id = i++;
	}

	/*
	 string filename = params.boot_trees;
	 filename += ".nolen";
	 boot_trees.printTrees(filename.c_str(), false);
	 return;
	 */
	SplitGraph sg;
	SplitIntMap hash_ss;
	// make the taxa name
	vector<string> taxname;
	taxname.resize(mytree.leafNum);
	mytree.getTaxaName(taxname);

	// read the bootstrap tree file
	double scale = 100.0;
	if (params->scaling_factor > 0)
		scale = params->scaling_factor;

	MTreeSet boot_trees;
	if (params && detectInputFile((char*)input_trees) == IN_NEXUS)
	{
		sg.init(*params);
		for (SplitGraph::iterator it = sg.begin(); it != sg.end(); it++)
			hash_ss.insertSplit((*it), (*it)->getWeight());
		StrVector sgtaxname;
		sg.getTaxaName(sgtaxname);
		i = 0;
		for (StrVector::iterator sit = sgtaxname.begin();
			sit != sgtaxname.end(); sit++, i++)
		{
			Node* leaf = mytree.findLeafName(*sit);
			if (!leaf)
				outError("Tree does not contain taxon ", *sit);
			leaf->id = i;
		}
		scale /= sg.maxWeight();
	}
	else
	{
		boot_trees.init(input_trees, rooted, burnin, max_count,
			tree_weight_file);
		boot_trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1);
		scale /= boot_trees.sumTreeWeights();
	}
	// sg.report(cout);
	cout << "Rescaling split weights by " << scale << endl;
	if (params->scaling_factor < 0)
		sg.scaleWeight(scale, true);
	else
	{
		sg.scaleWeight(scale, false, params->numeric_precision);
	}

	cout << sg.size() << " splits found" << endl;
	// compute the percentage of appearance
	//	printSplitSet(sg, hash_ss);
	// sg.report(cout);
	cout << "Creating bootstrap support values..." << endl;
	mytree.createBootstrapSupport(taxname, boot_trees, sg, hash_ss);
	// mytree.scaleLength(100.0/boot_trees.size(), true);
	string out_file;
	if (output_tree)
		out_file = output_tree;
	else
	{
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = target_tree;
		out_file += ".suptree";
	}

	mytree.printTree(out_file.c_str());
	cout << "Tree with assigned bootstrap support written to " << out_file
		<< endl;
	/*
	if (out_prefix)
		out_file = out_prefix;
	else
		out_file = target_tree;
	out_file += ".supval";
	mytree.writeInternalNodeNames(out_file);

	cout << "Support values written to " << out_file << endl;
	*/
}

void computeConsensusTree(const char* input_trees, int burnin, int max_count,
	double cutoff, double weight_threshold, const char* output_tree,
	const char* out_prefix, const char* tree_weight_file, Params* params)
{
	bool rooted = false;

	// read the bootstrap tree file
	/*
	 MTreeSet boot_trees(input_trees, rooted, burnin, tree_weight_file);
	 string first_taxname = boot_trees.front()->root->name;
	 //if (params.root) first_taxname = params.root;

	 SplitGraph sg;

	 boot_trees.convertSplits(sg, cutoff, SW_COUNT, weight_threshold);*/

	 // sg.report(cout);
	SplitGraph sg;
	SplitIntMap hash_ss;
	// make the taxa name
	// vector<string> taxname;
	// taxname.resize(mytree.leafNum);
	// mytree.getTaxaName(taxname);

	// read the bootstrap tree file
	double scale = 100.0;
	if (params->scaling_factor > 0)
		scale = params->scaling_factor;

	MTreeSet boot_trees;
	if (params && detectInputFile((char*)input_trees) == IN_NEXUS)
	{
		char* user_file = params->user_file;
		params->user_file = (char*)input_trees;
		params->split_weight_summary = SW_COUNT; // count number of splits
		sg.init(*params);
		params->user_file = user_file;
		for (SplitGraph::iterator it = sg.begin(); it != sg.end(); it++)
			hash_ss.insertSplit((*it), (*it)->getWeight());
		/*		StrVector sgtaxname;
		 sg.getTaxaName(sgtaxname);
		 i = 0;
		 for (StrVector::iterator sit = sgtaxname.begin(); sit != sgtaxname.end(); sit++, i++) {
		 Node *leaf = mytree.findLeafName(*sit);
		 if (!leaf) outError("Tree does not contain taxon ", *sit);
		 leaf->id = i;
		 }*/
		scale /= sg.maxWeight();
	}
	else
	{
		boot_trees.init(input_trees, rooted, burnin, max_count,
			tree_weight_file);
		boot_trees.convertSplits(sg, cutoff, SW_COUNT, weight_threshold);
		scale /= boot_trees.sumTreeWeights();
		cout << sg.size() << " splits found" << endl;
	}
	// sg.report(cout);
	if (verbose_mode >= VB_MED)
		cout << "Rescaling split weights by " << scale << endl;
	if (params->scaling_factor < 0)
		sg.scaleWeight(scale, true);
	else
	{
		sg.scaleWeight(scale, false, params->numeric_precision);
	}

	// cout << "Creating greedy consensus tree..." << endl;
	MTree mytree;
	SplitGraph maxsg;
	sg.findMaxCompatibleSplits(maxsg);

	if (verbose_mode >= VB_MAX)
		maxsg.saveFileStarDot(cout);
	// cout << "convert compatible split system into tree..." << endl;
	mytree.convertToTree(maxsg);
	// cout << "done" << endl;
	string taxname;
	if (params->root)
		taxname = params->root;
	else
		taxname = sg.getTaxa()->GetTaxonLabel(0);
	Node* node = mytree.findLeafName(taxname);
	if (node)
		mytree.root = node;
	// mytree.scaleLength(100.0 / boot_trees.sumTreeWeights(), true);

	// mytree.getTaxaID(maxsg.getSplitsBlock()->getCycle());
	// maxsg.saveFile(cout);

	string out_file;

	if (output_tree)
		out_file = output_tree;
	else
	{
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = input_trees;
		out_file += ".contree";
	}

	//	if (removed_seqs.size() > 0)
	//		mytree.insertTaxa(removed_seqs, twin_seqs);

	mytree.printTree(out_file.c_str(), WT_BR_CLADE);
	cout << "Consensus tree written to " << out_file << endl;

	if (output_tree)
		out_file = output_tree;
	else
	{
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = input_trees;
		out_file += ".splits";
	}

	// sg.scaleWeight(0.01, false, 4);
	if (params->print_splits_file)
	{
		sg.saveFile(out_file.c_str(), IN_OTHER, true);
		cout << "Non-trivial split supports printed to star-dot file " << out_file << endl;
	}
}

void testCompConsensus(const char* infile, const char* outfile, Params* params)
{
	StringIntMap treels;
	int index = 0;
	ifstream fin(infile);
	string tree;
	while (getline(fin, tree))
	{
		treels[tree] = index;
		index++;
	}
	cout << "Total trees read = " << index << endl;
	fin.close();
	IntVector weight;
	weight.resize(treels.size(), 1);
	string con_str = computeConsensusTreeNoFileIO(treels, weight, params->tree_max_count,
		params->split_threshold, params->split_weight_threshold, params);

	ofstream fout(outfile);
	fout << con_str;
	fout.close();

	cout << "All DONE" << endl;
}

string computeConsensusTreeNoFileIO(StringIntMap& input_trees, IntVector& weight, int max_count,
	double cutoff, double weight_threshold, Params* params)
{
	bool rooted = false;

	SplitGraph sg;
	SplitIntMap hash_ss;

	// read the bootstrap tree file
	double scale = 100.0;
	if (params->scaling_factor > 0)
		scale = params->scaling_factor;

	MTreeSet boot_trees;
	//	IntVector weight;
	//	weight.resize(input_trees.size(), 1);

	params->split_weight_summary = SW_COUNT;
	boot_trees.init(input_trees, rooted, weight);
	boot_trees.convertSplits(sg, cutoff, SW_COUNT, weight_threshold);
	scale /= boot_trees.sumTreeWeights();
	//	cout << sg.size() << " splits found" << endl;

	// sg.report(cout);
	if (verbose_mode >= VB_MED)
		cout << "Rescaling split weights by " << scale << endl;
	if (params->scaling_factor < 0)
		sg.scaleWeight(scale, true);
	else
	{
		sg.scaleWeight(scale, false, params->numeric_precision);
	}

	// cout << "Creating greedy consensus tree..." << endl;
	MTree mytree;
	SplitGraph maxsg;
	sg.findMaxCompatibleSplits(maxsg);

	mytree.convertToTree(maxsg);

	// cout << "done" << endl;
	string taxname;
	if (params->root)
		taxname = params->root;
	else
		taxname = sg.getTaxa()->GetTaxonLabel(0);
	Node* node = mytree.findLeafName(taxname);
	if (node)
		mytree.root = node;

	ostringstream ostr;
	string tree_str;
	//    mytree.printTree(ostr, WT_BR_CLADE); // for testing by testCompConsensus
	mytree.printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
	tree_str = ostr.str();
	return tree_str;
}

void computeConsensusNetwork(const char* input_trees, int burnin, int max_count,
	double cutoff, int weight_summary, double weight_threshold, const char* output_tree,
	const char* out_prefix, const char* tree_weight_file)
{
	bool rooted = false;

	// read the bootstrap tree file
	MTreeSet boot_trees(input_trees, rooted, burnin, max_count,
		tree_weight_file);

	SplitGraph sg;
	// SplitIntMap hash_ss;

	boot_trees.convertSplits(sg, cutoff, weight_summary, weight_threshold);

	string out_file;

	if (output_tree)
		out_file = output_tree;
	else
	{
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = input_trees;
		out_file += ".nex";
	}

	sg.saveFile(out_file.c_str(), IN_NEXUS);
	cout << "Consensus network printed to " << out_file << endl;

	if (output_tree)
		out_file = output_tree;
	else
	{
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = input_trees;
		out_file += ".splits";
	}
	if (verbose_mode >= VB_MED)
	{
		sg.saveFile(out_file.c_str(), IN_OTHER, true);
		cout << "Non-trivial split supports printed to star-dot file " << out_file << endl;
	}
}

/**
 * Diep:
 * This function is to check whether tree->aln has duplicate patterns in it
 * @return TRUE if there exists duplicate pattern, FALSE otherwise
 */
bool checkDuplicatePattern(IQTree*& tree)
{
	int nptn = tree->aln->getNPattern();
	bool found = false;
	for (int i = 0; i < nptn; i++)
	{
		Pattern pati = tree->aln->at(i);
		for (int j = i + 1; j < nptn; j++)
		{
			Pattern patj = tree->aln->at(j);
			int nchars = pati.size();
			int k;
			for (k = 0; k < nchars; k++)
			{
				if (patj[k] != pati[k])
					break;
			}
			if (k == nchars)
			{
				cout << "FOUND: pattern#" << i << " == pattern#" << j << endl;
				found = true;
			}
		}
	}
	return found;
}

void optimizeAlignment(IQTree*& tree, Params& params)
{
	//	if(checkDuplicatePattern(tree))
	//		cout << "FIRST CHECK: Alignment patterns are not created properly!" << endl;

	double start = getCPUTime();
	tree->params = &params; // Diep: 2020-08-17, there are two variables with identical name as 'params'

	//	tree->initTopologyByPLLRandomAdition(params); // this pll version needs further sync to work with the rest
	tree->computeParsimonyTree(params.out_prefix, tree->aln); // this iqtree version plays nicely with the rest
	// extract the vector of pattern pars of the initialized tree
	tree->initializeAllPartialPars();
	tree->clearAllPartialLH();
	//	tree->fixNegativeBranch(true);
	int pars_before = tree->computeParsimony();
	tree->curScore = pars_before;
	cout << "Time for parsimony tree construction: " << getCPUTime() - start << " seconds" << endl;
	cout << "Parsimony score: " << pars_before << endl;
	BootValTypePars* tmpPatternPars = tree->getPatternPars();
	for (int i = 0; i < tree->getAlnNPattern(); i++)
	{
		(tree->aln)->at(i).ras_pars_score = tmpPatternPars[i];
	}

	if (params.sort_alignment)
	{
		cout << "Reordering patterns in alignment by decreasing order of pattern parsimony... ";
		start = getCPUTime();
		// reordering patterns
		PatternComp pcomp;
		sort(tree->aln->begin(), tree->aln->end(), pcomp);
		tree->aln->updateSitePatternAfterOptimized();

		tree->initializeAllPartialPars();
		tree->clearAllPartialLH();
		int pars_after = tree->computeParsimony();
		if (pars_after != pars_before)
			outError("Reordering alignment has bug.");
		cout << getCPUTime() - start << " seconds" << endl;
	}
	else
	{
		tree->aln->updateSitePatternAfterOptimized();
	}

	tree->doSegmenting();

	//	if(checkDuplicatePattern(tree))
	//		cout << "SECOND CHECK: Sorted alignment patterns are duplicate!" << endl;
	//	exit(0);

	//	for(int i = 0; i < tree->getAlnNPattern(); i++){
	//		if(i == 253 || i == 254)
	//			cout << "i = " << i << ", score = " << (tree->aln)->at(i).ras_pars_score << endl;
	//	}
	//	exit(0);

	//	string tree_after = tree->getTreeString();
	//	cout << "TREE BEFORE: " << tree_before << endl;
	//	cout << "TREE AFTER: " << tree_after << endl;
	//	Alignment * sortedAln = new Alignment;
	//	*sortedAln = *(tree->aln);
	//	delete tree;
	//	tree = new IQTree(sortedAln);
	//
	//	tree->aln->printPhylip("temp.phy", false);
	//	Alignment * alignment = new Alignment("temp.phy", params.sequence_type, params.intype);
	////	tree->setAlignment(alignment);
	//	IQTree * tree1 = new IQTree(alignment);
	//	tree = tree1;
}

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

void addMoreRowMutation(IQTree* tree, Alignment* alignment)
{
	auto startTime = getCPUTime();
	IQTree newTree;
	char* file_name = "tree.treefile";
	bool is_rooted = false;
	newTree.readTree(file_name, is_rooted);

	newTree.setAlignment(tree->aln);
	newTree.aln = new Alignment;
	newTree.aln->copyAlignment(tree->aln);
	newTree.aln->ungroupSitePattern();
	newTree.save_branch_states_dad = new UINT[newTree.aln->size() + 1];
	newTree.save_branch_fitch_result = new UINT[newTree.aln->size() + 1];
	newTree.add_row = true;
	newTree.aln->missingSamples = alignment->missingSamples;
	newTree.aln->existingSamples = alignment->existingSamples;
	newTree.aln->reference_nuc = alignment->reference_nuc;

	cout << "tree parsimony before add k rows: " << tree->computeParsimony() << " " << newTree.computeParsimony() << '\n';
	cout << "ungroup alignment: " << tree->aln->getNSite() << " " << newTree.aln->getNSite() << '\n';
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

	for(int i = 0; i < (int)permCol.size(); ++i)
	{
		while((int)pos.size() <= permCol[i])
			pos.push_back(0);
		pos[permCol[i]] = i;
	}

	for(int i = 0; i < (int)alignment->missingSamples.size(); ++i)
    {
        for(auto m : alignment->missingSamples[i])
        {
            assert(newTree.aln->reference_nuc[m.position] > 0);
        }
    }
    for(int i = 0; i < (int)alignment->existingSamples.size(); ++i)
    {
        for(auto m : alignment->existingSamples[i])
        {
            assert(newTree.aln->reference_nuc[m.position] > 0);
        }
    }
	
	newTree.initMutation(permCol);

	cout << "tree parsimony after init mutations: " << newTree.computeParsimony() << " " << newTree.computeParsimonyScoreMutation() << '\n';
	int num_sample = (int)alignment->missingSamples.size();
	vector<MutationNode> missingSamples(num_sample);
	for(int i = 0; i < (int)alignment->missingSamples.size(); ++i)
	{
		missingSamples[i].mutations = alignment->missingSamples[i];
		for(auto m : alignment->missingSamples[i])
		{
			assert((m.ref_nuc & (m.ref_nuc-1)) == 0);
		}
		missingSamples[i].name = alignment->missingSamples[i][0].name;
	}

	// newTree.checkMutation(pos);
	// cout << "correct mutations\n\n";

	for(int i = 0; i < (int)missingSamples.size(); ++i)
	{
		vector<pair<PhyloNode *, PhyloNeighbor *> > bfs = newTree.breadth_first_expansion();
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
		int best_set_difference = (int)1e9+7;
		int set_difference = (int)1e9+7;
		size_t best_j = 0;
		size_t best_distance = (size_t)1e9+7;

		for(int j = 0; j < (int)bfs.size(); ++j)
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
			inp.best_j =  &best_j;
			inp.best_distance = &best_distance;
			inp.j = j;
			inp.has_unique = &best_node_has_unique;
			inp.node_has_unique = &(node_has_unique);

			newTree.calculatePlacementMutation(pos, inp, true, true);
		}
	
		newTree.addNewSample(bfs[best_j].first, bfs[best_j].second, node_excess_mutations[best_j], i, missingSamples[i].name);
		// newTree.aln->addToAlignmentNewSeq(missingSamples[i].name, alignment->remainSeq[i], savePermCol);
		// newTree.checkMutation(pos);
		// cout << newTree.computeParsimonyScoreMutation() << " " << newTree.computeParsimonyScore() << '\n';
	}
	delete newTree.aln;
	newTree.aln = NULL;
	cout << "new tree's parsimony score: " << newTree.computeParsimonyScoreMutation() << '\n';
	cout << "we still alive !\n";
	cout << "Time: " << fixed << setprecision(3) << (double)(getCPUTime() - startTime) << " seconds\n";
}