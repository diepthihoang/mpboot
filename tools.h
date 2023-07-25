/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
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

#ifndef TOOLS_H
#define TOOLS_H

#include <iqtree_config.h>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

//#include <sys/time.h>
//#include <time.h>
#include <sys/stat.h>
#include <math.h>
#include "ncl/ncl.h"
#include "msetsblock.h"

#define SPRNG
#include "sprng/sprng.h"


#define USE_HASH_MAP

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#else
#define GCC_VERSION 0
#endif

// for MSVC
#ifndef __func__
#define __func__ __FUNCTION__
#endif

#if defined(USE_HASH_MAP) && !defined(_MSC_VER)
	#if !defined(__GNUC__)
		#include <hash_map>
		#include <hash_set>
		using namespace stdext;
	#elif GCC_VERSION < 40300
		#include <ext/hash_map>
		#include <ext/hash_set>
		using namespace __gnu_cxx;
		#define unordered_map hash_map
		#define unordered_set hash_set
	#else
		#include <tr1/unordered_map>
		#include <tr1/unordered_set>
		using namespace std::tr1;
	#endif
#else
	#include <map>
	#include <set>
	#include <unordered_map>
	#include <unordered_set>
#endif

using namespace std;


#if	defined(USE_HASH_MAP) && GCC_VERSION < 40300 && !defined(_MSC_VER)
/*
        Define the hash function of Split
 */
#if !defined(__GNUC__)
namespace stdext {
#else
namespace __gnu_cxx {
#endif

    template<>
    struct hash<string> {

        size_t operator()(string str) const {
            hash<const char*> hash_str;
            return hash_str(str.c_str());
        }
    };
} // namespace
#endif // USE_HASH_MAP


class Linear {
public:

    Linear(int n, double *x, double *y) {

        // calculate the averages of arrays x and y
        double xa = 0, ya = 0;
        for (int i = 0; i < n; i++) {
            xa += x[i];
            ya += y[i];
        }
        xa /= n;
        ya /= n;

        // calculate auxiliary sums
        double xx = 0, yy = 0, xy = 0;
        for (int i = 0; i < n; i++) {
            double tmpx = x[i] - xa, tmpy = y[i] - ya;
            xx += tmpx * tmpx;
            yy += tmpy * tmpy;
            xy += tmpx * tmpy;
        }

        // calculate regression line parameters

        // make sure slope is not infinite
        assert(fabs(xx) != 0);

        m_b = xy / xx;
        m_a = ya - m_b * xa;
        m_coeff = (fabs(yy) == 0) ? 1 : xy / sqrt(xx * yy);

    }

    double getValue(double x) {
        return m_a + m_b * x;
    }

    //! Returns the slope of the regression line

    double getSlope() {
        return m_b;
    }

    //! Returns the intercept on the Y axis of the regression line

    double getIntercept() {
        return m_a;
    }

    //! Returns the linear regression coefficient

    double getCoefficient() {
        return m_coeff;
    }

private:

    double m_a, m_b, m_coeff;
};

/**
        vector of double number
 */
typedef vector<double> DoubleVector;

/**
        vector of int
 */
typedef vector<int> IntList;


/**
        vector of int
 */
typedef vector<int> IntVector;

/**
 		set of int
 */
typedef set<int> IntegerSet; // Diep added

/**
        vector of bool
 */
typedef vector<bool> BoolVector;


/**
        vector of char
 */
typedef vector<char> CharVector;

/**
        vector of string
 */
typedef vector<string> StrVector;


/**
        matrix of double number
 */
#define matrix(T) vector<vector<T> >

/**
        matrix of double
 */
/*
class DoubleMatrix {
public:
        double *value;
        int rows, cols, size;
        DoubleMatrix(int arows, int acols);
        //inline double operator() (int i, int j);
        inline double &operator() (int i, int j) { return value[i * cols + j]; };
        inline double *operator[] (int i) {	return value + (i*cols); };
        virtual ~DoubleMatrix();
        void setZero();
};
 */
typedef matrix(double) DoubleMatrix;

typedef unsigned int UINT;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
        run mode of program
 */
enum RunMode {
    DETECTED, GREEDY, PRUNING, BOTH_ALG, EXHAUSTIVE, DYNAMIC_PROGRAMMING,
    CALC_DIST, PD_USER_SET, PRINT_TAXA, PRINT_AREA, SCALE_BRANCH_LEN,
    SCALE_NODE_NAME, PD_DISTRIBUTION, LINEAR_PROGRAMMING, STATS //, GBO, MPRO
}; //STATS and GBO added by MA (STATS for some statistics on tree, GBO = guided 'bootstrap'

/**
        type of generating trees or splits graphs
 */
enum TreeGenType {
    NONE, YULE_HARDING, UNIFORM, CATERPILLAR, BALANCED,
    CIRCULAR_SPLIT_GRAPH, TAXA_SET, STAR_TREE
};

/**
        when writing tree:
                WT_BR_LEN - output branch length
                WT_BR_CLADE - put branch length into internal node name
                WT_TAXON_ID - output taxon ID
                WT_INT_NODE - for draw tree, draw the internal node ID
                WT_BR_SCALE - for draw tree, draw the branch proportional to its length
                WT_SORT_TAXA - sort the taxa s.t. subtrees with least taxon ID come first
                WT_APPEND    - append the output file
                WT_NEWLINE   - print a newline after
                WT_BR_LEN_FIXED_WIDTH - print branch length in fixed number format
 */
const int WT_BR_LEN = 1;
const int WT_BR_CLADE = 2;
const int WT_TAXON_ID = 4;
const int WT_INT_NODE = 8;
const int WT_BR_SCALE = 16;
const int WT_SORT_TAXA = 32;
const int WT_APPEND = 64;
const int WT_NEWLINE = 128;
const int WT_BR_LEN_FIXED_WIDTH = 256;
const int WT_BR_ID = 512;
const int WT_BR_LEN_ROUNDING = 1024;
const int TRUE = 1;
const int FALSE = 0;

/**
 *  Specify different ways of doing an NNI.
 *  TOPO_ONLY: only change the tree topology
 *  TOPO_UPDATE_LV: the same as above but the partial likelihoods are update in addition
 *  NNI1: optimize the central branch after changing the tree topology
 *  NNI5: optimized the 5 affected branches after changing the tree topology
 */
enum NNI_Type {
    TOPO_ONLY,
    TOPO_UPDATE_LV,
    NNI1,
    NNI5
};

/**
        when computing Robinson-Foulds distances
 */
const int RF_ADJACENT_PAIR = 1;
const int RF_ALL_PAIR = 2;
const int RF_TWO_TREE_SETS = 3;
const int RF_TWO_TREE_SETS_EXTENDED = 4; // work for trees with non-equal taxon sets

/**
        split weight summarization
 */
const int SW_COUNT = 1; // just counting the number of splits
const int SW_SUM = 2; // take the sum of all split weights
const int SW_AVG_ALL = 3; // take the split weight average over all trees
const int SW_AVG_PRESENT = 4; // take the split weight average over all trees that the split is present

/**
        search mode
 */
//enum SearchMode {EXHAUSTIVE, EXHAUSTIVE_CIRCULAR};

/**
        input type, tree or splits graph
 */
enum InputType {
    IN_NEWICK, IN_NEXUS, IN_FASTA, IN_PHYLIP, IN_OTHER
};

/**
        verbose mode, determine how verbose should the screen be printed.
 */
enum VerboseMode {
    VB_QUIET, VB_MIN, VB_MED, VB_MAX, VB_DEBUG
};

/**
        verbose level on the screen
 */
extern VerboseMode verbose_mode;

/**
        consensus reconstruction type
 */
enum ConsensusType {
    CT_NONE, CT_CONSENSUS_TREE, CT_CONSENSUS_NETWORK,
    CT_ASSIGN_SUPPORT, CT_ASSIGN_SUPPORT_EXTENDED, COMPARE
};

enum TestType {
    TEST_NONE, TEST_COMPATIBLE, TEST_CIRCULAR, TEST_WEAKLY_COMPATIBLE, TEST_K_COMPATIBLE
};

/**
        State frequency type
 */
enum StateFreqType {
    FREQ_UNKNOWN, FREQ_USER_DEFINED, FREQ_EQUAL, FREQ_EMPIRICAL, FREQ_ESTIMATE,
    FREQ_CODON_1x4, FREQ_CODON_3x4, FREQ_CODON_3x4C // special frequency for codon model
};

/**
        alignment format type
 */

enum AlnFormat {
    ALN_PHYLIP, ALN_FASTA
};

enum ModelTestCriterion {
    MTC_AIC, MTC_AICC, MTC_BIC
};

/**
        Stopping condition type
 */
enum STOP_CONDITION {
    SC_FIXED_ITERATION, SC_WEIBULL, SC_UNSUCCESS_ITERATION, SC_BOOTSTRAP_CORRELATION, SC_REAL_TIME
};

enum IQP_ASSESS_QUARTET {
    IQP_DISTANCE, IQP_PARSIMONY, IQP_BOOTSTRAP
};

enum LEAST_SQUARE_VAR {
    OLS, WLS_FIRST_TAYLOR, WLS_FITCH_MARGOLIASH, WLS_SECOND_TAYLOR, WLS_PAUPLIN
};

enum START_TREE_TYPE {
	STT_BIONJ, STT_PARSIMONY, STT_PLL_PARSIMONY
};

const int MCAT_LOG = 1; // categorize by log(rate) for Meyer & von Haeseler model
const int MCAT_MEAN = 2; // take the mean of rates for each category for Meyer & von Haeseler model
const int MCAT_PATTERN = 4; // categorize site-patterns instead of sites for Meyer & von Haeseler model

const double MAX_GENETIC_DIST = 9.0;

struct NNIInfo {
    double lh_score[4]; // tree log-likelihood of zero-branch, current tree, NNI tree 1, NNI tree 2
    double br_len[4]; // length of current branch, optimized branch, NNI branch 1, NNI branch 2
    int nni_round;
    int iqpnni_iteration;
};

enum LikelihoodKernel {
	LK_NORMAL, LK_SSE, LK_EIGEN, LK_EIGEN_SSE
};

enum LhMemSave {
	LM_DETECT, LM_ALL_BRANCH, LM_PER_NODE
};

/** maximum number of newton-raphson steps for NNI branch evaluation */
extern int NNI_MAX_NR_STEP;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
        program parameters, everything is specified here
 */
struct Params {

	/**
	 *  Number of starting parsimony trees
	 */
	int numParsTrees;

	/**
	 *  SPR distance (radius) for parsimony tree
	 */
	int sprDist;

	/**
	 *  Number of NNI trees generated from the set of parsimony trees
	 *  Default = 20 (out of 100 parsimony trees)
	 */
	int numNNITrees;

	/**
	 *  Population size
	 */
	int popSize;

	/**
	 *  maximum number of trees stored in the candidate set
	 */
	int maxCandidates;

	/**
	 *  heuristics for speeding up NNI evaluation
	 */
	bool speednni;

	bool adaptPert;

	/**
	 *  portion of NNI used for perturbing the tree
	 */
	double initPerStrength;

	/**
	 *  logl epsilon for the final model parameter optimization
	 */
	double modeps;

	/**
	 *  Carry out iterated local search using NNI only.
	 */
	bool snni;

	/**
	 *  only evaluate NNIs in affected regions
	 */
	bool fastnni;

    NNI_Type nni_type;

    /**
     *  Different type of Least Square variances
     */
	LEAST_SQUARE_VAR ls_var_type;

	/**
	 *  Threshold (likelihood difference between NNI and current tree)
	 *  to start optimizing 5 branches
	 */
	double nniThresHold;

	/**
	 *  Optimize 5 branches on NNI tree
	 */
	bool nni5;


    /**
     *  Number of smoothTree iteration carried out in Phylolib for IQP Tree
     */
    int numSmoothTree;

    /**
     *   compute least square branches for a given tree
     */
    bool leastSquareBranch;

    /** TRUE to apply Manuel's analytic approximation formulae for branch length */
    bool manuel_analytic_approx;

    /** TRUE to compute parsimony branch length of final tree */
    bool pars_branch_length;

    /** TRUE to compute bayesian branch length for the final tree */
    bool bayes_branch_length;

    /**
     *  use Least Square to evaluate NNI
     */
    bool leastSquareNNI;

    /**
     *  epsilon value used to compare log-likelihood between trees
     */
    double loglh_epsilon;

    /*
     *  reinsert leaves back to tree using parsimony
     */
    bool reinsert_par;

    /*
     *  Option to evaluate 10 different starting tree and take the best
     */
    bool bestStart;

    /**
     *  Maximum running time of the tree search in minutes
     */
    double maxtime;

    /**
     *  Turn on parsimony branch length estimation
     */
    bool parbran;

    /**
     *  option to turn on phylogenetic library
     */
    bool pll;

    /**
     *  OBSOLETE! Stopping rule for the tree search
     */
//    bool autostop;

    /**
     *  Number of maximum unsuccessful iterations after the search is stopped.
     *  Used for the automatic stopping rule
     */
    int unsuccess_iteration;

    char *binary_aln_file;

    /**
     *  the speed up heuristic will be used after
     *  speedup_iter iteration
     */
    int speedup_iter;

    /**
     *  starting CPU time of the program
     */
    double startCPUTime;

    /** starting real time of the program */
    double start_real_time;

    /**
     *		write all current best trees to file
     */
    bool write_best_trees;
    /**
     *  Number iteration = num_taxa * iteration_multiple
     */
    int iteration_multiple;
    /**
             input file name
     */
    char *user_file;

    /* type of starting tree */
    START_TREE_TYPE start_tree;

    /**
            prefix of the output file, default is the same as input file
     */
    char *out_prefix;

    /**
            alignment file name
     */
    char *aln_file;

    /**
            file containing multiple trees to evaluate at the end
     */
    char *treeset_file;

    /** number of bootstrap replicates for tree topology test */
    int topotest_replicates;

    /** true to perform weighted SH and KH test */
    bool do_weighted_test;

    /** true to do the approximately unbiased (AU) test */
    bool do_au_test;

    /**
            file specifying partition model
     */
    char *partition_file;

    /**
     * 		defines the relation between edge lengths in supertree and subtrees
     * 		0 (NULL) for separate edge length (default)
     * 		'p' for proportional edge length
     * 		'j' for joint edge length
     */
    char partition_type;

    /** remove all-gap sequences in partition model to account for terrace default: TRUE */
    bool remove_empty_seq;

    /** use terrace aware data structure for partition models, default: TRUE */
    bool terrace_aware;

    /**
            B, D, or P for Binary, DNA, or Protein sequences
     */
    char *sequence_type;

    /**
            alignment output file name
     */
    char *aln_output;

    /**
            file containing site likelihood as input for 'guided bootstrap' (added by MA)
     */
    char *siteLL_file;

    /**
            alignment where the gappy patterns will be superimposed into the input alignment
     */
    char *gap_masked_aln;

    /**
            alignment to be concatenated into the input alignment
     */
    char *concatenate_aln;

    /**
            file containing list of sites posititon to keep, format:
            pos1 pos2
            ....
     */
    char *aln_site_list;

    /**
            name of the reference sequence where aln_site_list is based on,
            NULL to take alignment positions.
     */
    char *ref_seq_name;

    /**
            alignment output format
     */
    AlnFormat aln_output_format;

    /**
            TRUE to discard all gappy positions
     */
    bool aln_nogaps;

    /**
     * TRUE to discard all constant sites
     */
    bool aln_no_const_sites;

    /**
            OBSOLETE compute parsimony score on trees
     */
//    bool parsimony;

    /**
            compute random step-wise addition parsimony tree instead of BIONJ
     */
//    bool parsimony_tree;

    /**
             output file name
     */
    char *out_file;

    /**
             size of the maximal PD-tree
     */
    int sub_size;

    /**
             min size of the maximal PD-tree
             used to calculate all PD-k trees from min_size to sub_size
     */
    int min_size;

    /**
            step_size when running from min_size to sub_size
     */
    int step_size;

    /**
            conservation proprotion, another way of input set size
     */
    double pd_proportion;

    /**
            min conservation proprotion
     */
    double min_proportion;

    /**
            step conservation proprotion
     */
    double step_proportion;

    /**
            sample size for computing PD distribution
     */
    int sample_size;


    /**
            TRUE if want to find all optimal PD-k set
            with the same maximal PD score
     */
    bool find_all;

    /**
             type of random tree to be generated
     */
    TreeGenType tree_gen;

    /**
            when generating random split graph, specify the number of
            splits here!
     */
    int num_splits;

    /**
             running mode: which algorithms to be applied
     */
    RunMode run_mode;

    /**
             real running mode if run_mode == DETECTED
     */
    RunMode detected_mode;

    /**
             parameter file
     */
    char *param_file;

    /**
            file containing taxa names to be included into the PD-tree
     */
    char *initial_file;

    /**
            file containing area names to be included into the PD set
     */
    char *initial_area_file;

    /**
            file containing a list of specific taxa sets which user wants
            to compute PD score on these sets only
     */
    char *pdtaxa_file;

    /**
            sets relation file, in form of a distance matrix file
     */
    char *areas_boundary_file;

    /**
            boundary length modifier
     */
    double boundary_modifier;

    /**
            output file to store the distance matrix
     */
    char *dist_file;

    /**
            TRUE to compute the observed distances instead of Juke-Cantor distances, default: FALSE
     */
    bool compute_obs_dist;

    /**
            TRUE to compute the Juke-Cantor distances, default: FALSE
     */
    bool compute_jc_dist;

    /**
            TRUE to compute the maximum-likelihood distances
     */
    bool compute_ml_dist;

    /**
            TRUE to compute the maximum-likelihood tree
     */
    bool compute_ml_tree;

    /**
            file containing budget information
     */
    char *budget_file;

    /**
            used when generating pair of taxa set with overlapping
     */
    int overlap;

    // private use
    /**
             number of times to repeat the algorithms
     */
    int repeated_time;

    /**
             print no tree to output
     */
    int nr_output;

    /**
            input type, tree or splits graph
     */
    InputType intype;

    /**
            total budget, for cost constrained PD problem
     */
    int budget;

    /**
            minimum budget, for cost constrained PD problem
     */
    int min_budget;

    /**
            step_budget when running from min_budget to budget
     */
    int step_budget;

    /**
            name of the root taxon
     */
    char *root;

    /**
            true if tree is forced to be rooted
     */
    bool is_rooted;


    /**
            min branch length, used to create random tree/network
     */
    double min_len;

    /**
            mean branch length, used to create random tree/network
     */
    double mean_len;

    /**
            max branch length, used to create random tree/network
     */
    double max_len;

    /**
            number of internal branches to set zero length
     */
    int num_zero_len;

    /**
            random number seed
     */
    unsigned int ran_seed;

    /**
            run time of the algorithm
     */
    double run_time;

    /**
            limit on the number of optimal PD sets
     */
    int pd_limit;

    /**
            TRUE if one wants to calculate the PD gain matrix in terms of delta_k^j = pd(PD_k \/ {j}) - pd_k
     */
    bool calc_pdgain;

    /**
            TRUE if tree file contains more than 1 tree
     */
    bool multi_tree;

    /**
            2nd user tree used in assignBootstrapSupport
     */
    char *second_tree;

    /**
            2nd alignment used in computing multinomialProb (Added by MA)
     */
    char *second_align;
    /**
            type of consensus building
     */
    ConsensusType consensus_type;

    /**
            file containing weights for every tree in the input tree file
     */
    char *tree_weight_file;

    /**
            set the TRUE if want to find the minimal PD set, instead of the default maximal PD set
     */
    bool find_pd_min;

    /**
            set TRUE to find area's endemic PD instead of regular PD
     */
    bool endemic_pd;

    /**
            set TRUE to find exclusive PD instead of regular PD
     */
    bool exclusive_pd;

    /**
            to find PD complementarity given this area
     */
    char *complement_area;

    /**
            used for likelihood mapping: for each branch, print the four cluster
     */
    int branch_cluster;

    /**
            file containing taxa order
     */
    char *taxa_order_file;

    /**
            to scale branch length or clade support with a factor
     */
    double scaling_factor;

    /**
            TRUE if always use binary linear programming
     */
    bool binary_programming;

    /**
            test the input split system in one of the TestType
     */
    TestType test_input;

    /**
            burnin value: number of beginning trees to be discarded
     */
    int tree_burnin;

    /**
            maximum number of trees to consider (for e.g. consensus tree construction)
     */
    int tree_max_count;

    /**
        threshold of split frequency, splits appear less than threshold will be discarded
     */
    double split_threshold;

    /**
            threshold of split weight, splits with weight less than or equal to threshold will be discarded
     */
    double split_weight_threshold;

    /**
            Way to summarize split weight in the consensus tree or network: SW_SUM, SW_AVG_ALL, or SW_AVG_PRESENT
     */
    double split_weight_summary;

    /**
            TRUE if use quadratic programming (for GUROBI)
     */
    bool quad_programming;

    /**
            true if one wants to optimize tree by subtree pruning and regrafting
     */
    bool tree_spr;

    /**
            true if printing out of optimal sets in NEXUS format
     */
    bool nexus_output;

    /**
            k-representative parameter, used for IQP algorithm
     */
    int k_representative;

    /**
            probability of deleting a leaf, used for IQP algorithm
     */
    double p_delete;

    /**
            min number of iqpnni iterations
     */
    int min_iterations;

    /**
            max number of iqpnni iterations
     */
    int max_iterations;

    /**
            stop condition, SC_FIXED_ITERATION or SC_STOP_PREDICT
     */
    STOP_CONDITION stop_condition;

    /**
            confidence value for stop rule
     */
    double stop_confidence;

    /**
            name of the substitution model (e.g., HKY, GTR, TN+I+G, JC+G, etc.)
     */
    string model_name;

    /** set of models for testing */
    char *model_set;

    /**
            TRUE to store transition matrix into a hash table for computation efficiency
     */
    bool store_trans_matrix;

    /**
            state frequency type
     */
    StateFreqType freq_type;


    /**
            the number of rate categories
     */
    int num_rate_cats;

    /**
            shape parameter (alpha) of the Gamma distribution for site rates
     */
    double gamma_shape;

    /**
            TRUE to use median rate for discrete categories, FALSE to use mean rate instead
     */
    bool gamma_median;

    /**
            proportion of invariable sites
     */
    double p_invar_sites;

    /** TRUE to optimize all model and rate parameters jointly by BFGS, default: FALSE */
    bool optimize_model_rate_joint;

    /**
            TRUE if you want to optimize branch lengths by Newton-Raphson method
     */
    bool optimize_by_newton;

    /**
            TRUE if you want to fix branch lengths during model optimization
     */
    bool fixed_branch_length;

    /**
            criterion to assess important quartet
     */
    IQP_ASSESS_QUARTET iqp_assess_quartet;

    /**
     *      Using IQP algorithm to do tree perturbation
     */
    bool iqp;

    /**
            the LP file is in gurobi format or not
     */
    bool gurobi_format;

    /**
            number of threads for gurobi call
     */
    bool gurobi_threads;

    /**
            TRUE if doing bootstrap on the input trees (good, bad, ugly)
     */
    int num_bootstrap_samples;

    /** bootstrap specification of the form "l1:b1,l2:b2,...,lk:bk"
        to randomly draw b1 sites from the first l1 sites, etc. Note that l1+l2+...+lk
        must equal m, where m is the alignment length. Otherwise, an error will occur.
        The default bootstrap_spec == NULL, a standard procedure is applied, i.e., randomly draw m sites.
    */
    char *bootstrap_spec;

    /**
            1 if output all intermediate trees from every IQPNNI iteration
            2 if output all intermediate trees + 1-NNI-away trees
     */
    int write_intermediate_trees;

    /**
     *  Write out all candidate trees (the locally optimal trees)
     */
    int write_local_optimal_trees;

    /**
        TRUE to avoid duplicated trees while writing intermediate trees
     */
    bool avoid_duplicated_trees;

    /**
            Robinson-Foulds distance computation mode: RF_ADJACENT PAIR, RF_ALL_PAIR
     */
    int rf_dist_mode;

    /**
            compute the site-specific rates by Meyer & von Haeseler method
     */
    bool mvh_site_rate;

    /**
            FALSE to use MH Model, FALSE for using tree-likelihood
     */
    bool rate_mh_type;

    /**
            TRUE to discard saturated for Meyer & von Haeseler (2003) model
     */
    bool discard_saturated_site;

    /**
            rates will be normalized to this mean value
     */
    double mean_rate;

    /**
            Percentage threshold to accept a branch of the approximate likelihood ratio test
            (aLRT) with SH-like interpretation. See Guindon et al. (2010) Syst. Biol. for details.
            Default: 90%.
     */
    int aLRT_threshold;

    /**
            number of replicates, default: 1000
     */
    int aLRT_replicates;

    /**
            number of replicates for local bootstrap probabilities method of Adachi & Hasegawa (1996) in MOLPHY
     */
    int localbp_replicates;

    /**
            SSE Option
     */
    LikelihoodKernel SSE;
    /**
     	 	0: do not print anything
            1: print site log-likelihood
            2: print site log-likelihood per Gamma category
     */
    int print_site_lh;

    /**
     	 	0: do not print anything
            1: print site parsimony
            2: print site parsimony per Gamma category // not yet
     */
    int print_site_pars; // 2020-08-17 Diep added

    /**
     	 	true: print site parsimony for user input tree
     */
    bool print_site_pars_user_tree; // 2020-08-17 Diep added


    /** TRUE to print site-specific rates, default: FALSE */
    bool print_site_rate;

    /**
            TRUE to print tree log-likelihood
     */
    bool print_tree_lh;

    bool print_branch_lengths;

    /****** adaptive NNI search heuristic ******/

    /**
     *  Output log-likelihood
     */
    bool nni_lh;

    /**
     *  The number of iqp iteration before the heuristics is applied
     */
    int speedUpFromIter;

    /**
     *  Lambda in PhyML algorithm
     */
    double lambda;

    /**
     * Confidence level for the speed up heuristics
     */
    double speed_conf;

    bool new_heuristic;

    /***** WH-test (Weiss & von Haeseler 2003) *****/

    /**
            Results of Weiss & Haeseler test of model homogeneity
     */
    double whtest_simulations;
    double whtest_delta;
    double whtest_delta_quantile;
    double whtest_p_value;


    /**
            bit-wise type including MCAT_LOG, MCAT_MEAN
     */
    int mcat_type;

    /**
            initial rate file in format:
            Site Rate
            1  f_1
            2  f_2
            ...
     */
    char *rate_file;

    /***** NGS stuffs   ********/

    /**
            next-generation sequencing input file for Fritz project
     */
    char *ngs_file;

    /**
            next-generation sequencing input file containing mapped reads to the reference genome
     */
    char *ngs_mapped_reads;

    bool ngs_ignore_gaps;

    bool do_pars_multistate;

    /**
            File containing p-values of the genes, for GSS project with Roland
     */
    char *gene_pvalue_file;

    /**
            scaling factor for the p-values
     */
    double gene_scale_factor;

    /**
            transforming pvalues to logarithms
     */
    bool gene_pvalue_loga;

    /***** variables for reading NCBI taxonomy tree *******/

    /**
            NCBI taxonomy ID, for processing nodes.dmp file
     */
    int ncbi_taxid;

    /**
            NCBI taxon rank, restricting the tree to that rank
     */
    const char *ncbi_taxon_level;

    /**
            rank to ingore, e.g., "no rank", branch length to such node will be set to zero
     */
    const char *ncbi_ignore_level;

    /**
            typically names.dmp from NCBI
     */
    const char *ncbi_names_file;

    /**********************************************/
    /******* variables for ECOpd analysis *********/

	/**
		eco_dag_file - contains the food web in matrix form (n species, nxn matrix), 0 for no connection, 1 for predation of j predator on i prey
	*/
	char *eco_dag_file;

    /**
		eco_detail_file - contains IDs of species present in the final set and/or species absent in the TREE or SPLIT system, but present in the food web
	*/
	const char *eco_detail_file;

	/*
	 * the type of the phylo input - tree or network
	 */
	const char *eco_type;

	/*
		k% - percent of species to be conserved
	 */
	int k_percent;

    /*
		diet - percent of species diet to be preserved for species survival
	*/
	int diet_min;
	int diet_max;
	int diet_step;

    /*
		eco_run - run number, used when random branch length is assigned to the edges of an input tree
	*/
	int eco_run;

    /*
		eco_weighted - indicates whether to treat the food web as weighted or not weighted
	*/
	bool eco_weighted;

    /**********************************************/
    /****** variables for upper bound tests *******/
	bool upper_bound;



    /**********************************************/
    /**** variables for ultra-fast bootstrap ******/

    /**
            number of replicates for guided bootstrap
     */
    int gbo_replicates;

	/* interval (l-epsilon,l+epsilon) indicates tie for bootstrap tree
	 * in this case, one tree is picked up at random
	 */
	double ufboot_epsilon;

    /**
            TRUE to check with different max_candidate_trees
     */
    int check_gbo_sample_size;

    /**
            TRUE to use RELL method of Simodaira Hasegawa, FALSE otherwise
     */
    bool use_rell_method;

    /**
            TRUE to use ELW method of Strimmer & Rambaut for new bootstrap, FALSE otherwise
     */
    bool use_elw_method;

    /**
            TRUE to weight each bootstrap sample by its probability, FALSE otherwise
     */
    bool use_weighted_bootstrap;

    /**
            TRUE to use the single ML tree per bootstrap, FALSE to include several sup-optima
     */
    bool use_max_tree_per_bootstrap;

    /** maximum number of candidate trees to consider for new bootstrap */
    int max_candidate_trees;

    /** TRUE if user_file contains topologically distinct trees */
    bool distinct_trees;

    /** NEW: TRUE to update bootstrap trees during the search (do not store treels_ptnlh).
            FALSE to call runGuidedBootstrap() at the end */
    bool online_bootstrap;

    /** minimal correlation coefficient for bootstrap stopping rule */
    double min_correlation;

    /** number of iterations between bootstrap stopping rule check */
    int step_iterations;

    /** TRUE to store all candidate trees in memory */
    bool store_candidate_trees;

	/** true to print all UFBoot trees to a file */
	bool print_ufboot_trees;

    /****** variables for NNI cutoff heuristics ******/

    /**
            TRUE to empirically estimate nni_cutoff
     */
    bool estimate_nni_cutoff;

    /**
            logl difference with zero-branch tree, to cutoff before evaluating NNI
     */
    double nni_cutoff;

    /**
            sort the NNI before evaluating
     */
    bool nni_sort;

    /**
            Obsolete: TRUE to optimize 5 branches around NNI
     */
    //bool nni_opt_5branches;

    /** print some output info for NNI */
    bool testNNI;

    /** TRUE to do approximate NNIs with approximate branch lengths before a normal NNI */
    bool approximate_nni;


    /** TRUE to compress big file using zlib */
    bool do_compression;

    /**
            number of bootstrap samples for AvH curiosity
     */
    int avh_test;

    /**
            number of bootstrap samples for Arndt's bootstrap plot
     */
    int bootlh_test;

    /**
            partition definition for Arndt's bootstrap plot
     */
    char* bootlh_partitions;

    /** precision when printing out for floating-point number */
    int numeric_precision;

    /** file containing state-frequencies per site for site-specific state frequency model
     * each line has n+1 entries (n=number of states):
     * site_ID state1_freq state2_freq ... staten_freq
     * where site_ID is from 1 to m (m=number of sites)
     */
    char *site_freq_file;

    /*
     * Diep:
     * to search using maximum parsimony search criterion
     */
    bool maximum_parsimony;

	bool multiple_hits; // to store multiple MP trees per boot sample

	int store_top_boot_trees; // number of top trees to store per boot sample (e.g. top 10, top 5)
	int ratchet_iter; // # of iterations between two ratchet iterations
	int ratchet_wgt; // weight to add to selected site during ratchet
	int ratchet_percent; // 0-100; percent of informative sites to upweight during ratchet

    bool compute_parsimony; // to compute parsimony of a user tree
    bool newick_to_tnt; // to convert a tree in Newick format into TNT syntax
    bool newick_to_nexus; // to convert a tree in Newick format into TNT syntax

    /** TRUE to condense parsimony equivalent sites, default: false */
    bool condense_parsimony_equiv_sites;

    /*
     * Diep:
     * Name of file storing Sankoff cost matrix
     */
    char * sankoff_cost_file;

    /** TRUE if using unsigned short for pattern parsimony score */
    bool sankoff_short_int;

    /*
     * Diep: to optimize the parsimony tree using SPR implemented in PLL
     */
    bool spr_parsimony;

    /*
     * Diep: SPR's radius as required by PLL
     */
    int spr_mintrav;
    int spr_maxtrav;

    /*
     * Diep: option for comparing PLL site parsimony and IQTree
     */
    bool test_site_pars;

    /*
     * Diep: option for computing parsimony RELL with compiler's automatic vectorization
     */
    bool auto_vectorize;

    /*
     * Diep: option for sorting alignment by decreasing (score * freq) of pattern
     */
    bool sort_alignment;

	/*
	 * Diep: for determining logl_cutoff in master_pars_simple
	 */
	int cutoff_percent;

	/*
	 * Diep: to turn on the option for the 1st hill-climbing in ratchet iteration by NNI
	 */
	bool hclimb1_nni;

	/*
	 * Diep: to turn off saveCurrentTree() during the 1st hill-climbing in ratchet iteration
	 */
	bool no_hclimb1_bb;

	/*
	 * Diep: to allow optimize bootstrap trees at the end of doTreeSearch()
	 */
	bool optimize_boot_trees;

	/*
	 * Diep: to not collect any tree during bootstrap, applicable only if used with optimize_boot_trees
	 * -opt_btree -save_trees_off means not running any saveCurrentTree() but running optimizeBootTrees()
	 */
	bool save_trees_off;

	/*
	 * Diep: pickily collecting candidate trees
	 * Iteration 1 is treated special:
	 * 		- no bootstrap tree update, i.e. no storing trees into treels
	 * 		- use 100 best trees found in the iteration to calculate logl_cutoff
	 */
	bool minimize_iter1_candidates;

	/*
	 * Diep: to compute logl_cutoff by MINH's new method
	 * logl_cutoff = min {score of bootstrap trees on original alignment}
	 */
	bool cutoff_from_btrees;

	/*
	 * Diep: to only use best tree of each iteration for saveCurrentTree()
	 * Vinh's request to test with NNI
	 */
	bool ibest_as_cand;


	/*
	 * Vinh's idea: estimate RELL by just first 10% columns in the alignment
	 */
	bool do_first_rell;

	/*
	 * Diep
	 * To use SPR for refinement step
	 */
	int opt_btree_spr;

	/*
	 * Diep
	 * To pick top k trees for each boot aln to input to refinement
	 * k trees must come from k distinct iteration
	 */
	int distinct_iter_top_boot; // top k boot trees picked from distinct iterations

	/*
	 * Diep:
	 * if -distinct_iter_top_boot -top_boot_concensus
	 * then in the refinement step, for each boot aln
	 * + compute concensus of all top boots
	 * + hill-climb in this one tree only
	 */
	bool top_boot_concensus;

	/*
	 * Diep
	 * To use NNI (i.e. SPR 1) for refinement step
	 */
	bool opt_btree_nni;

	/*
	 * Diep:
	 * To call test function
	 */
	bool test_mode;

        /**
         * number of search trees saved
         */
        int savek;

	/*
	 * Diep:
	 * Use with -test_mode
     * Remove duplicate sequence from the input alignment
     * (this is to create input aln to compare with mpboot1 because
     * mpboot1 does not properly support -keep_ident)
	 */
	bool remove_dup_seq;

#ifdef _OPENMP
    int num_threads;
#endif

    /** either MTC_AIC, MTC_AICc, MTC_BIC */
    ModelTestCriterion model_test_criterion;

    /** sample size for AICc and BIC */
    int model_test_sample_size;

    /** root state, for Tina's zoombie domain */
    char *root_state;

	/**
	 * TRUE to print bootstrap alignments, default: false
	 */
	bool print_bootaln;

	/** true to print sub alignments of super alignment, default: false */
	bool print_subaln;

	/** print partition information */
	bool print_partition_info;

	/** TRUE to print concatenated alignment, default: false */
	bool print_conaln;

	/** true to count all distinct trees visited during tree search */
	bool count_trees;

	/* -1 (auto-detect): will be set to 0 if there is enough memory, 1 otherwise
	 * 0: store all partial likelihood vectors
	 * 1: only store 1 partial likelihood vector per node */
	LhMemSave lh_mem_save;

	/* TRUE to print .splits file in star-dot format */
	bool print_splits_file;
    
    /** TRUE (default) to ignore identical sequences and add them back at the end */
    bool ignore_identical_seqs;

    /** TRUE to write initial tree to a file (default: false) */
    bool write_init_tree;
};

/**
        related measures for PD
 */
struct PDRelatedMeasures {
    /**
            names of areas
     */
    vector<string> setName;

    /**
            pd scores of areas
     */
    DoubleVector PDScore;

    /**
            exclusive PD scores
     */
    DoubleVector exclusivePD;

    /**
            endemic pd scores of an area given all other areas
     */
    DoubleVector PDEndemism;

    /**
            pd-complementarity scores of an area given some provided area
     */
    DoubleVector PDComplementarity;

};



/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
        @return TRUE of ch is a control character (ascii <= 32)
 */
inline bool controlchar(char ch) {
    return ch <= 32;
}

inline bool is_newick_token(char ch) {
    return ch == ':' || ch == ';' || ch == ',' || ch == ')' || ch == '(' || ch == '[' || ch == ']';
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
        print error message then exit program
 */
//void outError(char *error);

/**
        print error message then exit program
 */
void outError(const char *error);

/**
        print error message then exit program
 */
void outError(string error);


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
        print double error messages then exit program
 */
void outError(const char *error, const char *msg);

/**
        print double error messages then exit program
 */
void outError(const char *error, string msg);

/**
        Output a warning message to screen
        @param error warning message
 */
void outWarning(const char *warn);
void outWarning(string warn);


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/


/**
        generate a random branch length under an exponential distribution
        with mean params.mean_len. Also make sure that the resulting
        length is in the range (params.min_len, params.max_len)
        @return the random branch length
 */
double randomLen(Params &params);

/**
        convert string to int, with error checking
        @param str original string
        @return the integer
 */
/**
        Compute the logarithm of the factorial of an integer number
        @param num: the number
        @return logarithm of (num! = 1*2*...*num)
 */
double logFac(const int num);

/**
 * Function to randomly select an element in a C++ container
 *
 * @param begin
 * @param end
 * @return
 */
template <typename I>
I random_element(I begin, I end);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*
        Error messages
 */
const char ERR_NO_TAXON[] = "Find no taxon with name ";
const char ERR_NO_AREA[] = "Find no area with name ";
const char ERR_NO_ROOT[] = "Root taxon not found: ";
const char ERR_ROOT_NET[] = "-root option is not available for network";
const char ERR_CONFLICT_ROOT[] = "Tree is already rooted, -o <taxon> is not allowed.";
const char ERR_DUPLICATED_TAXA[] = "Duplicated taxa name in the tree.";
const char ERR_FEW_TAXA[] = "Number of taxa must be greater than 2.";
const char ERR_NO_SPLITS[] = "No splits found!";
const char ERR_FEW_SPLITS[] = "Number of splits must be at least equal to the number of taxa";
const char ERR_NEG_BRANCH[] = "Negative branch length not allowed.";
const char ERR_NO_MEMORY[] = "Not enough memory!";

const char ERR_READ_INPUT[] = "File not found or incorrect input, pls check it again.";
const char ERR_UNEXPECTED_EOF[] = "Unexpected end of file.";
const char ERR_READ_ANY[] = "Unidentified error while reading file, pls check it carefully again.";
const char ERR_WRITE_OUTPUT[] = "Cannot write to file ";

const char ERR_NO_K[] = "You must specify the number of taxa in the PD set.";
const char ERR_TOO_SMALL_K[] = "Size of PD-set must be at least the size of initial set.";
const char ERR_NO_BUDGET[] = "Total budget is not specified or less than zero.";
const char ERR_TOO_SMALL_BUDGET[] = "Not enough budget to conserve the inital set of taxa.";

const char ERR_INTERNAL[] = "Internal error, pls contact authors!";

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
 * convert int to string
 * @param int
 * @return string
 */
string convertIntToString(int number);

string convertDoubleToString(double number);

/**
 *
 * @param SRC
 * @param DEST
 * @return bool
 */
bool copyFile(const char SRC[], const char DEST[]);

/**
 * Check if the file exists
 * @param strFilename
 * @return
 */
bool fileExists(string strFilename);

/**
        convert string to int, with error checking
        @param str original string
        @return the number
 */
int convert_int(const char *str) throw (string);

/**
        convert string to int, with error checking
        @param str original string
        @param end_pos end position
        @return the number
 */
int convert_int(const char *str, int &end_pos) throw (string);

/**
        convert comma-separated string to integer vector, with error checking
        @param str original string with integers separated by comma
        @param vec (OUT) integer vector
 */
void convert_int_vec(const char *str, IntVector &vec) throw (string);

/**
        convert string to double, with error checking
        @param str original string
        @return the double
 */
double convert_double(const char *str) throw (string);

/**
        convert string to double, with error checking
        @param str original string
        @param end_pos end position
        @return the double
 */
double convert_double(const char *str, int &end_pos) throw (string);

/**
 * Convert seconds to hour, minute, second
 * @param sec
 * @return string represent hour, minute, second
 */
string convert_time(const double sec);


/**
        convert a string to to range lower:upper:step_size with error checking
        @param str original string
        @param lower (OUT) lower bound of the range
        @param upper (OUT) upper bound of the range
        @param step_size (OUT) step size of the range
 */
void convert_range(const char *str, int &lower, int &upper, int &step_size) throw (string);

/**
        convert a string to to range lower:upper:step_size with error checking
        @param str original string
        @param lower (OUT) lower bound of the range
        @param upper (OUT) upper bound of the range
        @param step_size (OUT) step size of the range
 */
void convert_range(const char *str, double &lower, double &upper, double &step_size) throw (string);


/**
        read the file containing branch/split scaling factor and taxa weights
        @param params program parameters
        @param ntaxa total number of taxa
        @param scale (OUT) scaling factor
        @param tax_name (OUT) vector of taxa names
        @param tax_weight (OUT) vector of corresponding taxa weights
 */
void readWeightFile(Params &params, int ntaxa, double &scale, StrVector &tax_name, DoubleVector &tax_weight);

/**
        read the initial taxa set from the file
        @param params program parameters
        @param ntaxa number of taxa
        @param tax_name (OUT) vector of taxa names
 */
void readInitTaxaFile(Params &params, int ntaxa, StrVector &tax_name);

/**
        read the initial area set from the file
        @param params program parameters
        @param nareas number of areas
        @param area_name (OUT) vector of area names
 */
void readInitAreaFile(Params &params, int nareas, StrVector &area_name);


/**
        read a list of taxa set from a file, not in nexus format but as follows:
        n1
        tax-name-1
        ...
        tax-name-n1

        n2
        tax-name-1
        ...
        tax-name-n2
        ....

        @param filename file name
        @param sets (OUT) the returned sets of taxa
 */
void readTaxaSets(char *filename, MSetsBlock *sets);

/**
        read areas shared boundary file, in form of a standard distance matrix
        @param file_name file name
        @param areas the read sets block
        @param areas_shared_boundary (OUT) shared boundary length between areas.
                Diagonal elements represent the boundary length of single areas
 */
void readAreasBoundary(char *file_name, MSetsBlock *areas, double *areas_shared_boundary);

/**
        parse program argument into params
        @param argc number of arguments
        @param argv list of arguments
        @param params (OUT) program parameters
 */
void parseArg(int argc, char *argv[], Params &params);

/**
        detect the format of input file
        @param input_file file name
        @return
                IN_NEWICK if file in newick format,
                IN_NEXUS if in nexus format,
                IN_FASTA if in fasta format,
                IN_PHYLIP if in phylip format,
                IN_OTHER if file format unknown.
 */
InputType detectInputFile(char *input_file);

/**
        if file exists, ask user to overwrite it or not
        @param filename file name
        @return TRUE if agree to overwrite an existing file, or simply file does not exist
 */
bool overwriteFile(char *filename);

/**
        print usage information
        @param argv program arguments list
        @param full_command TRUE to print all available commands, FALSE to print normal usage dialog
 */
void usage(char* argv[], bool full_command);

/**
 *   Print a string into a file
 */
void printString2File(string myString, string filename);

/**
 * print usage for iq-tree
 * @param program arguments list
 * @param full_command TRUE to print all available commands, FALSE to print normal usage dialog
 */
void usage_iqtree(char* argv[], bool full_command);


/**
 * Diep:
 * print usage for MPBoot
 * @param program arguments list
 * @param full_command TRUE to print all available commands, FALSE to print normal usage dialog
 */
void usage_mpboot(char* argv[], bool full_command);

/**
        parse area name string, where names are separated by commas
        @param area_names a string of name
        @param areas (OUT) a set of name string
 */
void parseAreaName(char *area_names, set<string> &areas);

/**
 * generate 2 different random integer numbers smaller than a specific integer threshold
 * @param size integer threshold
 * @param &first first random integer number
 * @param @second second random integer number
 */
void get2RandNumb(const int size, int &first, int &second);

/*
inline double getCPUTime(clock_t startTime) {
        return double(clock() - startTime) / CLOCKS_PER_SEC;
}*/


/**
 *  Fills the range [first, last) with sequentially increasing values,
 *  starting with value and repetitively evaluating ++value.
 *  Introduced in C++11 --> this is a reimplementation
 */
template<class ForwardIterator, class T>
void iota( ForwardIterator first, ForwardIterator last, T value );

/**
        compute p-value for a chi-square value
        @param chi_square chi-square value
        @param df degree of freedom
        @return p-value
 */
double computePValueChiSquare(double x, int df);

/*--------------------------------------------------------------*/
/* random number generator */
/*--------------------------------------------------------------*/

/**
 * initialize the random number generator
 * @param seed seed for generator
 */
int init_random(int seed);

/**
 * finalize random number generator (e.g. free memory
 */
int finish_random();

/**
 * returns a random integer in the range [0; n - 1]
 * @param n upper-bound of random number
 */
int random_int(int n);

/**
 *  return a random integer in the range [a,b]
 */
//int randint(int a, int b);

/**
 * returns a random integer in the range [0; RAND_MAX - 1]
 * = random_int(RAND_MAX)
 */
int random_int();

/**
 * returns a random floating-point nuber in the range [0; 1)
 */
double random_double();

template <class T>
void my_random_shuffle (T first, T last)
{
	int n = last - first;
	for (int i=n-1; i>0; --i) {
		swap (first[i],first[random_int(i+1)]);
	}
}

/**
 * generic function for sorting by index
 */
template <class T>
void quicksort_index(T* arr, int* index, int left, int right) {
    int i = left, j = right, tmp2;
    T tmp, pivot = arr[(left + right) / 2];

    /* partition */
    while (i <= j) {
        while (arr[i] < pivot)
            i++;
        while (pivot < arr[j])
            j--;
        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            tmp2 = index[i];
            index[i] = index[j];
            index[j] = tmp2;
            i++;
            j--;
        }
    };

    /* recursion */
    if (left < j)
        quicksort_index(arr, index, left, j);
    if (i < right)
        quicksort_index(arr, index, i, right);
}

/**
 * generic function for sorting by index preseving entries in [first,last)
 * @param first first element
 * @param last last element
 * @param index (OUT) ordered index of elements from smallest to largest
 */
template <class T>
void sort_index(T* first, T* last, int *index) {
    T* x;
    int i;
    T* arr = new T[last - first];
    for (x = first, i = 0; x != last; x++, i++) {
        index[i] = i;
        arr[i] = *x;
    }
    assert(last - first == i);
    quicksort_index(arr, index, 0, (last - first) - 1);
    delete [] arr;
}

/**
 * print the header of summary file
 */
void summarizeHeader(ostream &out, Params &params, bool budget_constraint, InputType analysis_type);

/**
 * print footer of summary file
 */
void summarizeFooter(ostream &out, Params &params);

int calculateSequenceHash(string &seq); 

#endif
