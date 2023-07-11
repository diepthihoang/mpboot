//
// C++ Interface: phylotree
//
// Description:
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef PHYLOTREE_H
#define PHYLOTREE_H
// #define NDEBUG
//  comented out this for Mac

// PLEASE DONT TOUCH THESE VARIABLES ANYMORE!
#define EIGEN_NO_AUTOMATIC_RESIZING
// #define EIGEN_CACHEFRIENDLY_PRODUCT_THRESHOLD 32 // PLEASE DONT TOUCH THESE VARIABLES ANYMORE!
// #define EIGEN_UNROLLING_LIMIT 1000 // PLEASE DONT TOUCH THESE VARIABLES ANYMORE!

// #define EIGEN_TUNE_FOR_CPU_CACHE_SIZE (512*256)
// #define EIGEN_TUNE_FOR_CPU_CACHE_SIZE (8*512*512)
#include "Eigen/Core"
#include "mtree.h"
#include "alignment.h"
#include "model/modelsubst.h"
#include "model/modelfactory.h"
#include "phylonode.h"
#include "optimization.h"
#include "model/rateheterogeneity.h"
#include "phyloanalysis.h"
#include "mutation.h"

const double MIN_BRANCH_LEN = 0.000001; // NEVER TOUCH THIS CONSTANT AGAIN PLEASE!
const double MAX_BRANCH_LEN = 100.0;
const double TOL_BRANCH_LEN = 0.000001;       // NEVER TOUCH THIS CONSTANT AGAIN PLEASE!
const double TOL_LIKELIHOOD = 0.001;          // NEVER TOUCH THIS CONSTANT AGAIN PLEASE!
const double TOL_LIKELIHOOD_PARAMOPT = 0.001; // BQM: newly introduced for ModelFactory::optimizeParameters
// const static double SCALING_THRESHOLD = sqrt(DBL_MIN);
// const static double SCALING_THRESHOLD = 1e-100;
// const static double SCALING_THRESHOLD_INVER = 1 / SCALING_THRESHOLD;
// const static double LOG_SCALING_THRESHOLD = log(SCALING_THRESHOLD);
// #include "pll/pll.h"
//  2^256
#define SCALING_THRESHOLD_INVER 115792089237316195423570985008687907853269984665640564039457584007913129639936.0
#define SCALING_THRESHOLD (1.0 / SCALING_THRESHOLD_INVER)
#define LOG_SCALING_THRESHOLD log(SCALING_THRESHOLD)

const int SPR_DEPTH = 2;

using namespace Eigen;

#ifdef __AVX
#define MEM_ALIGNMENT 32
inline size_t get_safe_upper_limit(size_t cur_limit)
{
        return ((cur_limit + 3) / 4) * 4;
}

inline size_t get_safe_upper_limit_float(size_t cur_limit)
{
        return ((cur_limit + 7) / 8) * 8;
}

#else
#define MEM_ALIGNMENT 16
inline size_t get_safe_upper_limit(size_t cur_limit)
{
        return (cur_limit % 2 == 0) ? cur_limit : cur_limit + 1;
}
inline size_t get_safe_upper_limit_float(size_t cur_limit)
{
        return ((cur_limit + 3) / 4) * 4;
}

#endif

// #include "pll/mem_alloc.h"

inline double *aligned_alloc_double(size_t size)
{
#if defined WIN32 || defined _WIN32 || defined __WIN32__
        return (double *)_aligned_malloc(size * sizeof(double), MEM_ALIGNMENT);
#else
        void *res;
        posix_memalign(&res, MEM_ALIGNMENT, size * sizeof(double));
        return (double *)res;
#endif
}

template <class T>
inline T *aligned_alloc(size_t size)
{
#if defined WIN32 || defined _WIN32 || defined __WIN32__
        return (T *)_aligned_malloc(size * sizeof(T), MEM_ALIGNMENT);
#else
        void *res;
        posix_memalign(&res, MEM_ALIGNMENT, size * sizeof(T));
        return (T *)res;
#endif
}

inline void aligned_free(void *mem)
{
#if defined WIN32 || defined _WIN32 || defined __WIN32__
        _aligned_free(mem);
#else
        free(mem);
#endif
}

#ifdef __AVX
#define VectorClassMaster Vec4d
#define VectorClassFloat Vec8f
#define VectorClassInt Vec8i
#define VectorClassUShort Vec16us
#define VCSIZE_MASTER 4
#define VCSIZE_FLOAT 8
#define VCSIZE_INT 8
#define VCSIZE_USHORT 16

// #pragma message "Using AVX instructions"
#else
#define VectorClassMaster Vec2d
#define VectorClassFloat Vec4f
#define VectorClassInt Vec4i
#define VectorClassUShort Vec8us
#define VCSIZE_MASTER 2
#define VCSIZE_FLOAT 4
#define VCSIZE_INT 4
#define VCSIZE_USHORT 8

// #pragma message "Using SS3 instructions"
#endif

/**
 *  Row Major Array For Eigen
 */
typedef Array<double, Dynamic, Dynamic, RowMajor> RowMajorArrayXXd;

typedef std::map<string, double> StringDoubleMap;
typedef std::map<int, PhyloNode *> IntPhyloNodeMap;

#define MappedMat(NSTATES) Map<Matrix<double, NSTATES, NSTATES>>
#define MappedArr2D(NSTATES) Map<Array<double, NSTATES, NSTATES>>
#define MappedRowVec(NSTATES) Map<Matrix<double, 1, NSTATES>>
#define MappedVec(NSTATES) Map<Matrix<double, NSTATES, 1>>
#define Matrix(NSTATES) Matrix<double, NSTATES, NSTATES>
#define RowVector(NSTATES) Matrix<double, 1, NSTATES>
#define MappedRowArr2DDyn Map<Array<double, Dynamic, Dynamic, RowMajor>>
#define MappedArrDyn Map<Array<double, Dynamic, 1>>
#define MappedVecDyn(NSTATES) Map<Matrix<double, Dynamic, NSTATES>>

const int MAX_SPR_MOVES = 20;

/**
        an SPR move.
 */
struct SPRMove
{
        PhyloNode *prune_dad;
        PhyloNode *prune_node;
        PhyloNode *regraft_dad;
        PhyloNode *regraft_node;
        double score;
};

struct SPR_compare
{

        bool operator()(SPRMove s1, SPRMove s2) const
        {
                return s1.score > s2.score;
        }
};

class SPRMoves : public set<SPRMove, SPR_compare>
{
public:
        void add(PhyloNode *prune_node, PhyloNode *prune_dad,
                 PhyloNode *regraft_node, PhyloNode *regraft_dad, double score);
};

/*
left_node-----------dad-----------right_node
                     |
                     |
                     |inline
                    node
 */
struct PruningInfo
{
        NeighborVec::iterator dad_it_left, dad_it_right, left_it, right_it;
        Neighbor *dad_nei_left, *dad_nei_right, *left_nei, *right_nei;
        Node *node, *dad, *left_node, *right_node;
        double left_len, right_len;
        double *dad_lh_left, *dad_lh_right;
};

struct SwapNNIParam
{
        double nni1_score;
        double nni1_brlen;
        double nni2_score;
        double nni2_brlen;
        Neighbor *node1_nei;
        Neighbor *node2_nei;
        double *nni1_ptnlh;
        double *nni2_ptnlh;
};

struct NNIMove
{
        // Two nodes representing the central branch
        PhyloNode *node1, *node2;
        // Roots of the two subtree that are swapped
        NeighborVec::iterator node1Nei_it, node2Nei_it;

        // log-likelihood of the tree after applying the NNI
        double newloglh;

        int swap_id;

        // old branch lengths of 5 branches before doing NNI
        // double oldLen[5];

        // new branch lengths of 5 branches corresponding to the NNI
        double newLen[5];

        // pattern likelihoods
        double *ptnlh;

        bool operator<(const NNIMove &rhs) const
        {
                return newloglh > rhs.newloglh;
                // return delta > rhs.delta;
        }
};

struct LeafFreq
{
        int leaf_id;

        int freq;

        bool operator<(const LeafFreq &rhs) const
        {
                return (freq < rhs.freq);
        }
};

/**
 * PLEASE CALL THIS AT THE BEGINNING OF PROGRAM:
 * Precompute some information for fast Fitch parsimony algorithm
 */
void precomputeFitchInfo();

/**
Phylogenetic Tree class

        @author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
 */
class PhyloTree : public MTree, public Optimization
{

        friend class PhyloSuperTree;
        friend class PhyloSuperTreePlen;
        friend class RateGamma;
        friend class RateKategory;

public:
        /**
           default constructor ( everything is initialized to NULL)
         */
        PhyloTree();

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /**
         * Constructor with given alignment
         * @param alignment
         */
        PhyloTree(Alignment *aln);

        void init();

        /**
            destructor
         */

        virtual ~PhyloTree();

        /**
                copy the phylogenetic tree structure into this tree, override to take sequence names
                in the alignment into account
                @param tree the tree to copy
         */

        bool add_row;
        UINT *save_branch_states_dad, *save_branch_fitch_result;

        // init mutation for calculate
        void initMutation(vector<int> &permCol);

        // compute mutation for a branch
        void computeMutationBranch(vector<int> &permCol, PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst = NULL);

        void computePartialMutation(UINT *states_dad, vector<int> &permCol, PhyloNeighbor *dad_branch, PhyloNode *dad);

        // compute parsimony score using mutation
        int computeParsimonyScoreMutation();
        
        int computeParsimonyBranchMutation(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst = NULL);

        int computePartialParsimonyMutation(PhyloNeighbor *dad_branch, PhyloNode *dad);

        std::vector<std::pair<PhyloNode*, PhyloNeighbor*> > breadth_first_expansion();
        
        // calculate new mutations
        void calculatePlacementMutation(CandidateNode &input, bool compute_parsimony_scores = false, bool compute_vecs = false);

        void addNewSample(PhyloNode* best_node, PhyloNeighbor *best_node_branch, std::vector<Mutation> node_excess_mutations, int index, std::string name);

        void checkMutation(vector<int> &pos);

        void checkMutationBranch(vector<int> &pos, PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst = NULL);

        string checkPartialMutation(vector<int> &pos, PhyloNeighbor *dad_branch, PhyloNode *dad);

        virtual void copyTree(MTree *tree);
        /**
                copy the sub-tree structure into this tree
                @param tree the tree to copy
                @param taxa_set 0-1 string of length leafNum (1 to keep the leaf)
         */
        virtual void copyTree(MTree *tree, string &taxa_set);

        /**
                copy the phylogenetic tree structure into this tree, designed specifically for PhyloTree.
                So there is some distinction with copyTree.
                @param tree the tree to copy
         */
        void copyPhyloTree(PhyloTree *tree);

        /**
                set the alignment, important to compute parsimony or likelihood score
                @param alignment associated alignment
         */
        void setAlignment(Alignment *alignment);

        /** set the root by name */
        void setRootNode(char *my_root);

        /**
                set the substitution model, important to compute the likelihood
                @param amodel associated substitution model
         */
        void setModel(ModelSubst *amodel);

        /**
                set the model factory
                @param model_fac model factory
         */
        void setModelFactory(ModelFactory *model_fac);

        /**
                set rate heterogeneity, important to compute the likelihood
                @param rate associated rate heterogeneity class
         */
        void setRate(RateHeterogeneity *rate);

        /**
                get rate heterogeneity
                @return associated rate heterogeneity class
         */
        RateHeterogeneity *getRate();

        void discardSaturatedSite(bool val);

        /**
                get the name of the model
         */
        virtual string getModelName();

        /**
         * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}+I{pinvar}+G{alpha}
         */
        virtual string getModelNameParams();

        ModelSubst *getModel()
        {
                return model;
        }

        ModelFactory *getModelFactory()
        {
                return model_factory;
        }

        virtual bool isSuperTree()
        {
                return false;
        }

        virtual bool isParsimonyTree()
        {
                return false;
        }

        /**
                allocate a new node. Override this if you have an inherited Node class.
                @param node_id node ID
                @param node_name node name
                @return a new node
         */
        virtual Node *newNode(int node_id = -1, const char *node_name = NULL);

        /**
                allocate a new node. Override this if you have an inherited Node class.
                @param node_id node ID
                @param node_name node name issued by an interger
                @return a new node
         */
        virtual Node *newNode(int node_id, int node_name);

        /**
         *		@return number of alignment patterns
         */
        virtual int getAlnNPattern()
        {
                return aln->getNPattern();
        }

        /**
         *		@return number of alignment sites
         */
        virtual int getAlnNSite()
        {
                return aln->getNSite();
        }

        /**
                this function return the parsimony or likelihood score of the tree. Default is
                to compute the parsimony score. Override this function if you define a new
                score function.
                @return the tree score
         */
        // virtual double computeScore() { return -computeLikelihood(); }
        // virtual double computeScore() { return (double)computeParsimonyScore(); }

        /****************************************************************************
                Parsimony function
         ****************************************************************************/

        /**
         * 		Return the approximated branch length estimation using corrected parsimony branch length
         * 		This is usually used as the starting point before using Newton-Raphson
         */
        double computeCorrectedParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad);

        /**
                initialize partial_pars vector of all PhyloNeighbors, allocating central_partial_pars
         */
        virtual void initializeAllPartialPars();

        /**
                initialize partial_pars vector of all PhyloNeighbors, allocating central_partial_pars
                @param node the current node
                @param dad dad of the node, used to direct the search
                @param index the index
         */
        virtual void initializeAllPartialPars(int &index, PhyloNode *node = NULL, PhyloNode *dad = NULL);

        /**
                compute the tree parsimony score
                @return parsimony score of the tree
         */
        virtual int computeParsimony();

        /**
                TODO: Compute partial parsimony score of the subtree rooted at dad
                @param dad_branch the branch leading to the subtree
                @param dad its dad, used to direct the tranversal
         */
        virtual void computePartialParsimony(PhyloNeighbor *dad_branch, PhyloNode *dad);

        /**
                compute tree parsimony score on a branch
                @param dad_branch the branch leading to the subtree
                @param dad its dad, used to direct the tranversal
                @param branch_subst (OUT) if not NULL, the number of substitutions on this branch
                @return parsimony score of the tree
         */
        virtual int computeParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst = NULL);

        void printParsimonyStates(PhyloNeighbor *dad_branch = NULL, PhyloNode *dad = NULL);

        /**
                SLOW VERSION: compute the parsimony score of the tree, given the alignment
                @return the parsimony score
         */
        int computeParsimonyScore();

        /**
                SLOW VERSION: compute the parsimony score of the tree, given the alignment
                @return the parsimony score
                @param node the current node
                @param dad dad of the node, used to direct the search
                @param ptn pattern ID
                @param states set of admissible states at the current node (in binary code)
         */
        int computeParsimonyScore(int ptn, int &states, PhyloNode *node = NULL, PhyloNode *dad = NULL);

        /****************************************************************************
                likelihood function
         ****************************************************************************/

        /**
                initialize partial_lh vector of all PhyloNeighbors, allocating central_partial_lh
         */
        virtual void initializeAllPartialLh();

        /**
                de-allocate central_partial_lh
         */
        virtual void deleteAllPartialLh();

        /**
                initialize partial_lh vector of all PhyloNeighbors, allocating central_partial_lh
                @param node the current node
                @param dad dad of the node, used to direct the search
                @param index the index
         */
        virtual void initializeAllPartialLh(int &index, int &indexlh, PhyloNode *node = NULL, PhyloNode *dad = NULL);

        /**
                clear all partial likelihood for a clean computation again
         */
        virtual void clearAllPartialLH();

        /**
         * compute all partial likelihoods if not computed before
         */
        void computeAllPartialLh(PhyloNode *node = NULL, PhyloNode *dad = NULL);

        /**
                allocate memory for a partial likelihood vector
         */
        double *newPartialLh();

        /**
                allocate memory for a partial parsimony vector
         */
        virtual UINT *newPartialPars();

        /** get the number of bytes occupied by partial_lh */
        int getPartialLhBytes();

        /**
                allocate memory for a scale num vector
         */
        UBYTE *newScaleNum();

        /** get the number of bytes occupied by scale_num */
        int getScaleNumBytes();

        /**
         * this stores partial_lh for each state at the leaves of the tree because they are the same between leaves
         * e.g. (1,0,0,0) for A,  (0,0,0,1) for T
         */
        double *tip_partial_lh;
        bool tip_partial_lh_computed;

        void computeTipPartialLikelihood();
        void computePtnInvar();

        /**
                compute the partial likelihood at a subtree
                @param dad_branch the branch leading to the subtree
                @param dad its dad, used to direct the tranversal
         */
        virtual void computePartialLikelihood(PhyloNeighbor *dad_branch, PhyloNode *dad = NULL);

        void computePartialLikelihoodNaive(PhyloNeighbor *dad_branch, PhyloNode *dad = NULL);

        /**
         * this implements the SSE version using Eigen library
         */
        template <int NSTATES>
        inline void computePartialLikelihoodSSE(PhyloNeighbor *dad_branch, PhyloNode *dad = NULL);

        template <const int nstates>
        void computePartialLikelihoodEigen(PhyloNeighbor *dad_branch, PhyloNode *dad = NULL);

        template <class VectorClass, const int VCSIZE, const int nstates>
        void computePartialLikelihoodEigenTipSSE(PhyloNeighbor *dad_branch, PhyloNode *dad = NULL);

        /**
                compute tree likelihood on a branch. used to optimize branch length
                @param dad_branch the branch leading to the subtree
                @param dad its dad, used to direct the tranversal
                @param pattern_lh (OUT) if not NULL, the function will assign pattern log-likelihoods to this vector
                                assuming pattern_lh has the size of the number of patterns
                @return tree likelihood
         */
        virtual double computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh = NULL);

        /**
         * this implements the SSE version using Eigen library
         */
        template <int NSTATES>
        inline double computeLikelihoodBranchSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh = NULL);

        /**
         * MINH: this implements the fast alternative strategy for reversible model (March 2013)
         * where partial likelihoods at nodes store real partial likelihoods times eigenvectors
         */
        template <int NSTATES>
        inline double computeLikelihoodBranchFast(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh = NULL);

        template <const int nstates>
        double computeLikelihoodBranchEigen(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh = NULL);

        template <class VectorClass, const int VCSIZE, const int nstates>
        double computeLikelihoodBranchEigenTipSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh = NULL);

        //    double computeLikelihoodBranchNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh = NULL, double *pattern_rate = NULL);
        double computeLikelihoodBranchNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh = NULL);

        /**
                compute tree likelihood when a branch length collapses to zero
                @param dad_branch the branch leading to the subtree
                @param dad its dad, used to direct the tranversal
                @return tree likelihood
         */
        virtual double computeLikelihoodZeroBranch(PhyloNeighbor *dad_branch, PhyloNode *dad);

        /**
            compute likelihood of rooted tree with virtual root (FOR TINA)
            @param dad_branch the branch leading to the subtree
            @param dad its dad, used to direct the tranversal
            @return tree likelihood
         */
        virtual double computeLikelihoodRooted(PhyloNeighbor *dad_branch, PhyloNode *dad);

        /**
                compute the tree likelihood
                @param pattern_lh (OUT) if not NULL, the function will assign pattern log-likelihoods to this vector
                                assuming pattern_lh has the size of the number of patterns
                @return tree likelihood
         */
        virtual double computeLikelihood(double *pattern_lh = NULL);

        /**
                compute pattern likelihoods only if the accumulated scaling factor is non-zero.
                Otherwise, copy the pattern_lh attribute
                @param pattern_lh (OUT) pattern log-likelihoods,
                                assuming pattern_lh has the size of the number of patterns
                @param cur_logl current log-likelihood (for sanity check)
                @param pattern_lh_cat (OUT) if not NULL, store all pattern-likelihood per category
         */
        virtual void computePatternLikelihood(double *pattern_lh, double *cur_logl = NULL,
                                              double *pattern_lh_cat = NULL);

        /**
         * compute the array of pattern parsimony
         * to be serving the parsimony UFBoot
         * @param ptn_npars: array to store negative pattern parsimony
         * @cur_npars: pointer to the variable storing tree's negative parsimony score
         */
        virtual void computePatternParsimony(double *ptn_npars, double *cur_npars = NULL);

        /**
                Compute the variance in tree log-likelihood
                (Kishino & Hasegawa 1989, JME 29:170-179)
                @param pattern_lh pattern log-likelihoods, will be computed if NULL
                @param tree_lh tree log-likelihood, will be computed if ZERO
         */
        double computeLogLVariance(double *pattern_lh = NULL, double tree_lh = 0.0);

        /**
                Compute the variance in log-likelihood difference
                between the current tree and another tree.
                (Kishino & Hasegawa 1989, JME 29:170-179)
                @param pattern_lh_other pattern log-likelihoods of the other tree
                @param pattern_lh pattern log-likelihoods of current tree, will be computed if NULL
         */
        double computeLogLDiffVariance(double *pattern_lh_other, double *pattern_lh = NULL);

        /**
         *  \brief Estimate the observed branch length between \a dad_branch and \a dad analytically.
         *	The ancestral states of the 2 nodes are first computed (Yang, 2006).
         *	Branch length are then computed using analytical formula.
         *	@param[in] dad_branch must be an internal node
         *	@param[in] dad must be an internal node
         *	@return estimated branch length or -1.0 if one of the 2 nodes is leaf
         */
        double computeBayesianBranchLength(PhyloNeighbor *dad_branch, PhyloNode *dad);

        /**
         * \brief Approximate the branch legnth between \a dad_branch and \a dad using Least Square instead
         * of Newton Raphson
         * @param[in] dad_branch
         * @param[in] dad
         * @return approximated branch length
         */
        double computeLeastSquareBranLen(PhyloNeighbor *dad_branch, PhyloNode *dad);

        /**
         * Update all subtree distances that are affect by doing an NNI on branch (node1-node2)
         * @param nni NNI move that is carried out
         */
        void updateSubtreeDists(NNIMove &nni);

        /**
         * Compute all pairwise distance of subtree rooted at \a source and other subtrees
         */
        void computeSubtreeDists();

        void getUnmarkedNodes(PhyloNodeVector &unmarkedNodes, PhyloNode *node = NULL, PhyloNode *dad = NULL);

        void computeAllSubtreeDistForOneNode(PhyloNode *source, PhyloNode *nei1, PhyloNode *nei2, PhyloNode *node, PhyloNode *dad);

        double correctBranchLengthF81(double observedBran, double alpha = -1.0);

        double computeCorrectedBayesianBranchLength(PhyloNeighbor *dad_branch, PhyloNode *dad);

        /**
                Compute the variance in log-likelihood difference
                between the current tree and another tree.
                (Kishino & Hasegawa 1989, JME 29:170-179)
                @param other_tree the other tree to compare
                @param pattern_lh pattern log-likelihoods of current tree, will be computed if NULL
         */
        double computeLogLDiffVariance(PhyloTree *other_tree, double *pattern_lh = NULL);

        /**
                Roll back the tree saved with only Taxon IDs and branch lengths.
                For this function to work, one must printTree before with WT_TAXON_ID + WT_BR_LEN
                @param best_tree_string input stream to read from
         */
        void rollBack(istream &best_tree_string);

        /**
                Read the tree saved with Taxon Names and branch lengths.
                @param tree_string tree string to read from
         */
        void readTreeString(const string &tree_string);

        /**
                Read the tree saved with Taxon Names and branch lengths.
                @param tree_string tree string to read from
         */
        void readTreeFile(const string &file_name);

        /**
         * Return the tree string contining taxon names and branch lengths
         * @return
         */
        string getTreeString();

        /**
         *  Return the sorted topology without branch length, used to compare tree topology
         */
        string getTopology();

        bool checkEqualScalingFactor(double &sum_scaling, PhyloNode *node = NULL, PhyloNode *dad = NULL);

        /****************************************************************************
                computing derivatives of likelihood function
         ****************************************************************************/

        double computeLikelihoodDervNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

        /**
         * this implements the SSE version using Eigen library
         */
        template <int NSTATES>
        inline double computeLikelihoodDervSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

        template <const int nstates>
        double computeLikelihoodDervEigen(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

        template <class VectorClass, const int VCSIZE, const int nstates>
        double computeLikelihoodDervEigenTipSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

        /**
         * MINH: this implements the fast alternative strategy for reversible model (March 2013)
         * where partial likelihoods at nodes store real partial likelihoods times eigenvectors
         */
        template <int NSTATES>
        inline double computeLikelihoodDervFast(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

        /**
                compute tree likelihood and derivatives on a branch. used to optimize branch length
                @param dad_branch the branch leading to the subtree
                @param dad its dad, used to direct the tranversal
                @param df (OUT) first derivative
                @param ddf (OUT) second derivative
                @return tree likelihood
         */
        double computeLikelihoodDerv(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

        template <int NSTATES>
        double computeLikelihoodDervSSE_Test(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

        template <int NSTATES>
        double computeLikelihoodDervSSE_Test2(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

        template <int NSTATES>
        inline void computeLikelihoodDervSSE_PTN(int ptn, double *partial_lh_site, double *partial_lh_child, double *trans_state, double *derv1_state, double *derv2_state);
        /****************************************************************************
                Stepwise addition (greedy) by maximum parsimony
         ****************************************************************************/

        /**
                FAST VERSION: used internally by computeParsimonyTree() to find the best target branch to add into the tree
                @param added_node node to add
                @param target_node (OUT) one end of the best branch found
                @param target_dad (OUT) the other end of the best branch found
                @param node the current node
                @param dad dad of the node, used to direct the search
                @return the parsimony score of the tree
         */
        int addTaxonMPFast(Node *added_node, Node *&target_node, Node *&target_dad, Node *node, Node *dad);

        /**
         * FAST VERSION: compute parsimony tree by step-wise addition
         * @param out_prefix prefix for .parstree file
         * @param alignment input alignment
         */
        void computeParsimonyTree(const char *out_prefix, Alignment *alignment);

        /**
                SLOW VERSION: grow the tree by step-wise addition
                @param alignment input alignment
         */
        void growTreeMP(Alignment *alignment);

        /**
                used internally by growTreeMP() to find the best target branch to add into the tree
                @param added_node node to add
                @param target_node (OUT) one end of the best branch found
                @param target_dad (OUT) the other end of the best branch found
                @param node the current node
                @param dad dad of the node, used to direct the search
                @return the parsimony score of the tree
         */
        int addTaxonMP(Node *added_node, Node *&target_node, Node *&target_dad, Node *node, Node *dad);

        /****************************************************************************
                Nearest Neighbor Interchange with parsimony
         ****************************************************************************/
        /**
                search by a nearest neigbor interchange with parsimony
         */
        void searchNNI();

        /**
                search by a nearest neigbor interchange with parsimony
                @param node the current node
                @param dad dad of the node, used to direct the search
                @param cur_score current score
                @return best score
         */
        double searchNNI(double cur_score, PhyloNode *node = NULL, PhyloNode *dad = NULL);

        /**
                try to swap the tree with nearest neigbor interchange at the branch connecting node1-node2.
                If a swap shows better score, return the swapped tree and the score.
                @param cur_score current score
                @param node1 1st end node of the branch
                @param node2 2nd end node of the branch
                @return best score
         */
        double swapNNI(double cur_score, PhyloNode *node1, PhyloNode *node2);

        /****************************************************************************
                Branch length optimization by maximum likelihood
         ****************************************************************************/

        /**
                optimize one branch length by ML
                @param node1 1st end node of the branch
                @param node2 2nd end node of the branch
                @param clearLH true to clear the partial likelihood, otherwise false
                @param maxNRStep maximum number of Newton-Raphson steps
                @return likelihood score
         */
        virtual double optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH = true, int maxNRStep = 100);

        /**
                optimize all branch lengths of the children of node
                @param node the current node
                @param dad dad of the node, used to direct the search
                @return the likelihood of the tree
         */
        double optimizeChildBranches(PhyloNode *node, PhyloNode *dad = NULL);

        /**
                optimize all branch lengths at the subtree rooted at node step-by-step.
                @param node the current node
                @param dad dad of the node, used to direct the search
                @return the likelihood of the tree
         */
        virtual double optimizeAllBranches(PhyloNode *node, PhyloNode *dad = NULL, int maxNRStep = 100);

        /**
         * optimize all branch lengths at the subtree rooted at node step-by-step.
         * Using Least Squares instead of Newton Raphson.
         * @param node the current node
         * @param dad dad of the node, used to direct the search
         */
        void optimizeAllBranchesLS(PhyloNode *node = NULL, PhyloNode *dad = NULL);

        /**
                optimize all branch lengths of the tree
                @param iterations number of iterations to loop through all branches
                @return the likelihood of the tree
         */
        virtual double optimizeAllBranches(int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100);

        /**
                inherited from Optimization class, to return to likelihood of the tree
                when the current branch length is set to value
                @param value current branch length
                @return negative of likelihood (for minimization)
         */
        virtual double computeFunction(double value);

        /**
                Inherited from Optimization class.
                This function calculate f(value), first derivative f'(value) and 2nd derivative f''(value).
                used by Newton raphson method to minimize the function.
                @param value current branch length
                @param df (OUT) first derivative
                @param ddf (OUT) second derivative
                @return negative of likelihood (for minimization)
         */
        virtual double computeFuncDerv(double value, double &df, double &ddf);

        /****************************************************************************
               Branch length optimization by Least Squares
        ****************************************************************************/

        /**
         * Estimate the current branch using least squares
         * @param node1 first node of the branch
         * @param node2 second node of the branch
         * @return
         */
        double optimizeOneBranchLS(PhyloNode *node1, PhyloNode *node2);

        /****************************************************************************
                Auxilary functions and varialbes for speeding up branch length optimization (RAxML Trick)
         ****************************************************************************/

        bool theta_computed;

        /**
         *	NSTATES x NUMCAT x (number of patterns) array
         *	Used to store precomputed values when optimizing branch length
         *	See Tung's report on 07.05.2012 for more information
         */
        double *theta_all;

        /**
         * frequencies of alignment patterns, used as buffer for likelihood computation
         */
        double *ptn_freq;

        /**
         * used as buffer for faster likelihood computation
         * for const pattern: it stores product of p_invar and state frequency
         * for other pattern: zero
         */
        double *ptn_invar;

        /****************************************************************************
                Nearest Neighbor Interchange by maximum likelihood
         ****************************************************************************/

        /**
                search by a nearest neigbor interchange, then optimize branch lengths. Do it
                until tree does not improve
                @return the likelihood of the tree
         */
        double optimizeNNIBranches();

        /**
                search by a nearest neigbor interchange
                @return the likelihood of the tree
         */
        double optimizeNNI();

        /**
                search by a nearest neigbor interchange
                @param cur_score current likelihood score
                @param node the current node
                @param dad dad of the node, used to direct the search
                @return the likelihood of the tree
         */
        double optimizeNNI(double cur_score, PhyloNode *node = NULL, PhyloNode *dad = NULL
                           /*,ostream *out = NULL, int brtype = 0, ostream *out_lh = NULL, ostream *site_lh = NULL,
                   StringIntMap *treels = NULL, vector<double*> *treels_ptnlh = NULL, DoubleVector *treels_logl = NULL,
                   int *max_trees = NULL, double *logl_cutoff = NULL*/
        );

        /**
           search for the best NNI move corresponding to this branch
           @return NNIMove the best NNI, this NNI could be worse than the current tree
           according to the evaluation scheme in use
           @param node1 1 of the 2 nodes on the branch
           @param node2 1 of the 2 nodes on the branch
           @param nniMoves (IN/OUT) detailed information of the 2 NNIs, set .ptnlh to compute pattern likelihoods
         */
        virtual NNIMove getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, NNIMove *nniMoves = NULL);

        /**
                Do an NNI
                @param move reference to an NNI move object containing information about the move
                @param clearLH decides whether or not the partial likelihood should be cleared
         */
        virtual void doNNI(NNIMove &move, bool clearLH = true);

        /**
         * Randomly choose perform an NNI, out of the two defined by branch node1-node2.
         * This function also clear the corresponding partial likelihood vectors
         * @param node1 one node of the branch
         * @param node2 one node of the branch
         */
        void doOneRandomNNI(Node *node1, Node *node2);

        /**
         *   Apply 5 new branch lengths stored in the NNI move
         *   @param nnimove the NNI move currently in consideration
         */
        virtual void changeNNIBrans(NNIMove nnimove);

        /****************************************************************************
                Stepwise addition (greedy) by maximum likelihood
         ****************************************************************************/

        /**
                grow the tree by step-wise addition
                @param alignment input alignment
         */
        void growTreeML(Alignment *alignment);

        /**
                used internally by growTreeML() to find the best target branch to add into the tree
                @param added_node node to add
                @param target_node (OUT) one end of the best branch found
                @param target_dad (OUT) the other end of the best branch found
                @param node the current node
                @param dad dad of the node, used to direct the search
                @return the likelihood of the tree
         */
        double addTaxonML(Node *added_node, Node *&target_node, Node *&target_dad, Node *node, Node *dad);

        /****************************************************************************
                Distance function
         ****************************************************************************/

        /**
                compute the distance between 2 sequences.
                @param seq1 index of sequence 1
                @param seq2 index of sequence 2
                @param initial_dist initial distance
                @param (OUT) variance of distance between seq1 and seq2
                @return distance between seq1 and seq2
         */

        virtual double computeDist(int seq1, int seq2, double initial_dist, double &var);

        virtual double computeDist(int seq1, int seq2, double initial_dist);

        /**
                compute distance and variance matrix, assume dist_mat and var_mat are allocated by memory of size num_seqs * num_seqs.
                @param dist_mat (OUT) distance matrix between all pairs of sequences in the alignment
                @param var_mat (OUT) variance matrix for distance matrix
                @return the longest distance
         */
        double computeDist(double *dist_mat, double *var_mat);

        /**
                compute observed distance matrix, assume dist_mat is allocated by memory of size num_seqs * num_seqs.
                @param dist_mat (OUT) distance matrix between all pairs of sequences in the alignment
                @return the longest distance
         */
        double computeObsDist(double *dist_mat);

        /**
                compute distance matrix, allocating memory if necessary
                @param params program parameters
                @param alignment input alignment
                @param dist_mat (OUT) distance matrix between all pairs of sequences in the alignment
                @param dist_file (OUT) name of the distance file
                @return the longest distance
         */
        double computeDist(Params &params, Alignment *alignment, double *&dist_mat, double *&var_mat, string &dist_file);

        /**
                compute observed distance matrix, allocating memory if necessary
                @param params program parameters
                @param alignment input alignment
                @param dist_mat (OUT) distance matrix between all pairs of sequences in the alignment
                @param dist_file (OUT) name of the distance file
                @return the longest distance
         */
        double computeObsDist(Params &params, Alignment *alignment, double *&dist_mat, string &dist_file);

        /**
                correct the distances to follow metric property of triangle inequalities.
                Using the Floyd alogrithm.
                @param dist_mat (IN/OUT) the shortest path between all pairs of taxa
        @return the longest distance
         */
        double correctDist(double *dist_mat);

        /****************************************************************************
                compute BioNJ tree, a more accurate extension of Neighbor-Joining
         ****************************************************************************/

        /**
                compute BioNJ tree
                @param params program parameters
                @param alignment input alignment
                @param dist_file distance matrix file
         */
        void computeBioNJ(Params &params, Alignment *alignment, string &dist_file);

        /**
                Neighbor-joining/parsimony tree might contain negative branch length. This
                function will fix this.
                @param fixed_length fixed branch length to set to negative branch lengths
                @param node the current node
                @param dad dad of the node, used to direct the search
                @return The number of branches that have no/negative length
         */
        virtual int fixNegativeBranch(bool force = false, Node *node = NULL, Node *dad = NULL);

        // OBSOLETE: assignRandomBranchLengths no longer needed, use fixNegativeBranch instead!
        //    int assignRandomBranchLengths(bool force = false, Node *node = NULL, Node *dad = NULL);

        /* compute Bayesian branch lengths based on ancestral sequence reconstruction */
        void computeAllBayesianBranchLengths(Node *node = NULL, Node *dad = NULL);

        /****************************************************************************
                Subtree Pruning and Regrafting by maximum likelihood
                NOTE: NOT DONE YET
         ****************************************************************************/

        /**
                search by Subtree pruning and regrafting
                @return the likelihood of the tree
         */
        double optimizeSPR();

        /**
                search by Subtree pruning and regrafting, then optimize branch lengths. Iterative until
                no tree improvement found.
                @return the likelihood of the tree
         */
        double optimizeSPRBranches();

        /**
                search by Subtree pruning and regrafting at a current subtree
                @param cur_score current likelihood score
                @param node the current node
                @param dad dad of the node, used to direct the search
                @return the likelihood of the tree
         */
        double optimizeSPR(double cur_score, PhyloNode *node = NULL, PhyloNode *dad = NULL);

        /**
         *  original implementation by Minh
         */
        double optimizeSPR_old(double cur_score, PhyloNode *node = NULL, PhyloNode *dad = NULL);

        /**
         *  original implementation by Minh
         */
        double swapSPR_old(double cur_score, int cur_depth, PhyloNode *node1, PhyloNode *dad1,
                           PhyloNode *orig_node1, PhyloNode *orig_node2,
                           PhyloNode *node2, PhyloNode *dad2, vector<PhyloNeighbor *> &spr_path);

        /**
                move the subtree (dad1-node1) to the branch (dad2-node2)
         */
        double swapSPR(double cur_score, int cur_depth, PhyloNode *node1, PhyloNode *dad1,
                       PhyloNode *orig_node1, PhyloNode *orig_node2,
                       PhyloNode *node2, PhyloNode *dad2, vector<PhyloNeighbor *> &spr_path);

        double assessSPRMove(double cur_score, const SPRMove &spr);

        void pruneSubtree(PhyloNode *node, PhyloNode *dad, PruningInfo &info);

        void regraftSubtree(PruningInfo &info,
                            PhyloNode *in_node, PhyloNode *in_dad);

        /****************************************************************************
                Approximate Likelihood Ratio Test with SH-like interpretation
         ****************************************************************************/

        void computeNNIPatternLh(double cur_lh,
                                 double &lh2, double *pattern_lh2,
                                 double &lh3, double *pattern_lh3,
                                 PhyloNode *node1, PhyloNode *node2);

        /**
                Resampling estimated log-likelihood (RELL)
         */
        void resampleLh(double **pat_lh, double *lh_new);

        /**
                Test one branch of the tree with aLRT SH-like interpretation
         */
        double testOneBranch(
            double best_score, double *pattern_lh, int reps, int lbp_reps,
            PhyloNode *node1, PhyloNode *node2, double &lbp_support);

        /**
                Test all branches of the tree with aLRT SH-like interpretation
         */
        int testAllBranches(int threshold,
                            double best_score, double *pattern_lh, int reps, int lbp_reps,
                            PhyloNode *node = NULL, PhyloNode *dad = NULL);

        /****************************************************************************
                Collapse stable (highly supported) clades by one representative
         ****************************************************************************/

        /**
                delete a leaf from the tree, assume tree is birfucating
                @param leaf the leaf node to remove
         */
        void deleteLeaf(Node *leaf);

        /**
                reinsert one leaf back into the tree
                @param leaf the leaf to reinsert
                @param adjacent_node the node adjacent to the leaf, returned by deleteLeaves() function
                @param node one end node of the reinsertion branch in the existing tree
                @param dad the other node of the reinsertion branch in the existing tree
         */
        void reinsertLeaf(Node *leaf, Node *node, Node *dad);

        bool isSupportedNode(PhyloNode *node, int min_support);

        /**
                Collapse stable (highly supported) clades by one representative
                @return the number of taxa prunned
         */
        int collapseStableClade(int min_support, NodeVector &pruned_taxa, StrVector &linked_name, double *&dist_mat);

        int restoreStableClade(Alignment *original_aln, NodeVector &pruned_taxa, StrVector &linked_name);

        /**
                randomize the neighbor orders of all nodes
         */
        void randomizeNeighbors(Node *node = NULL, Node *dad = NULL);

        virtual void changeLikelihoodKernel(LikelihoodKernel lk);

        /****************************************************************************
                Public variables
         ****************************************************************************/

        /**
                associated alignment
         */
        Alignment *aln;

        /**
         * Distance matrix
         */
        double *dist_matrix;

        /**
         * Variance matrix
         */
        double *var_matrix;

        /**
                TRUE if you want to optimize branch lengths by Newton-Raphson method
         */
        bool optimize_by_newton;

        /**
         *      TRUE if the loglikelihood is computed using SSE
         */
        LikelihoodKernel sse;

        /**
         * Current score of the tree;
         */
        double curScore;

        /**
         * for UpperBounds: Initial tree log-likelihood
         */
        double mlInitial;

        /**
         * for UpperBounds: Log-likelihood after optimization of model parameters in the beginning of tree search
         */
        double mlFirstOpt;

        /*
         * for UpperBounds: mlCheck = 1, if previous two values were already saved.
         * Needed, because parameter optimization is done twice before and after tree search
         */

        int mlCheck;

        /*
         * 		Store the all the parameters for the program
         */
        Params *params;

        /**
                assign the leaf names with the alignment sequence names, using the leaf ID for assignment.
                @param node the starting node, NULL to start from the root
                @param dad dad of the node, used to direct the search
         */
        void assignLeafNames(Node *node = NULL, Node *dad = NULL);

        /**
         * initialize partition information for super tree
         */
        virtual void initPartitionInfo()
        {
        }

        /**
         * print transition matrix for all branches
         *
         */
        void printTransMatrices(Node *node = NULL, Node *dad = NULL);

        /**
         * compute the memory size required for storing partial likelihood vectors
         * @return memory size required in bytes
         */
        virtual uint64_t getMemoryRequired();

        /****** following variables are for ultra-fast bootstrap *******/
        /** 2 to save all trees, 1 to save intermediate trees */
        int save_all_trees;

        set<int> computeNodeBranchDists(Node *node = NULL, Node *dad = NULL);

        /*
         * Manuel's approach for analytic approximation of branch length given initial guess
            b0: initial guess for the maximum
            @return approximted branch length
        */
        double approxOneBranch(PhyloNode *node, PhyloNode *dad, double b0);

        void approxAllBranches(PhyloNode *node = NULL, PhyloNode *dad = NULL);

        //     void addRemainRow(const vector<string>& remainRowName, const vector<string>& remainRow, const vector<int>& perm, const vector<int> &permCol);

protected:
        /**
         *  is the subtree distance matrix need to be computed or updated
         */
        bool subTreeDistComputed;

        /**
         * Map data structure to store distance between subtree.
         * The key is a string which is constructed by concatenating IDs of
         * the 2 nodes, e.g. 15-16
         */
        StringDoubleMap subTreeDists;

        StringDoubleMap subTreeWeights;

        /** distance (# of branches) between 2 nodes */
        int *nodeBranchDists;

        /**
         * A list containing all the marked list. This is used in the dynamic programming
         * algorithm for compute inter subtree distances
         */
        IntPhyloNodeMap markedNodeList;

        /** converted root state, for Tina's zoombie domain */
        char root_state;

        /**
                internal pattern log-likelihoods, always stored after calling computeLikelihood()
                or related functions. Note that scaling factors are not incorporated here.
                If you want to get real pattern log-likelihoods, please use computePatternLikelihood()
         */
        double *_pattern_lh;

        /**
         * Store array of pattern parsimony computed in computeParsimonyBranch()
         */
        BootValTypePars *_pattern_pars;

        /**
                internal pattern likelihoods per category, always stored after calling computeLikelihood()
                or related functions. Note that scaling factors are not incorporated here.
                If you want to get real pattern likelihoods, please use computePatternLikelihood()
         */
        double *_pattern_lh_cat;

        /**
                associated substitution model
         */
        ModelSubst *model;

        /**
                Model factory includes SubstModel and RateHeterogeneity
                stores transition matrices computed before for efficiency purpose, eps. AA or CODON model.
         */
        ModelFactory *model_factory;

        /**
                among-site rates
         */
        RateHeterogeneity *site_rate;

        /**
                current branch iterator, used by computeFunction() to optimize branch lengths
                and by computePatternLikelihood() to compute all pattern likelihoods
         */
        PhyloNeighbor *current_it;
        /**
                current branch iterator of the other end, used by computeFunction() to optimize branch lengths
                and by computePatternLikelihood() to compute all pattern likelihoods
         */
        PhyloNeighbor *current_it_back;

        /**
                spr moves
         */
        SPRMoves spr_moves;

        /**
                SPR radius
         */
        int spr_radius;

        /**
                the main memory storing all partial likelihoods for all neighbors of the tree.
                The variable partial_lh in PhyloNeighbor will be assigned to a region inside this variable.
         */
        double *central_partial_lh;

        /**
                the main memory storing all scaling event numbers for all neighbors of the tree.
                The variable scale_num in PhyloNeighbor will be assigned to a region inside this variable.
         */
        UBYTE *central_scale_num;

        /**
                the main memory storing all partial parsimony states for all neighbors of the tree.
                The variable partial_pars in PhyloNeighbor will be assigned to a region inside this variable.
         */
        UINT *central_partial_pars;

        /**
                TRUE to discard saturated for Meyer & von Haeseler (2003) model
         */
        bool discard_saturated_site;

        /**
         * Temporary partial likelihood array: used when swapping branch and recalculate the
         * likelihood --> avoid calling malloc everytime
         */
        double *tmp_partial_lh1;
        double *tmp_partial_lh2;

        /**
         *  Temporary array containing anscentral states.
         *  Used to avoid calling malloc
         */

        double *tmp_anscentral_state_prob1;
        double *tmp_anscentral_state_prob2;
        /** pattern-specific rates */
        // double *tmp_ptn_rates;

        /**
         * Temporary scale num array: used when swapping branch and recalculate the
         * likelihood --> avoid calling malloc
         */
        UBYTE *tmp_scale_num1;
        UBYTE *tmp_scale_num2;

        /****************************************************************************
                Vector of bit blocks, used for parsimony function
         ****************************************************************************/

        /**
                @return size of the bits block vector for one node
         */
        size_t getBitsBlockSize();

        /**
                allocate new memory for a bit block vector
                @return the allocated memory
         */
        virtual UINT *newBitsBlock();

        /**
                @return size of the bits entry (for storing num_states bits)
         */
        int getBitsEntrySize();

        /**
                @param bits_entry
                @return TRUE if bits_entry contains all 0s, FALSE otherwise
         */
        bool isEmptyBitsEntry(UINT *bits_entry);

        /**
                @param bits_entry1
                @param bits_entry1
                @param bits_union (OUT) union of bits_entry1 and bits_entry2
         */
        void unionBitsEntry(UINT *bits_entry1, UINT *bits_entry2, UINT *&bits_union);

        /**
                set a single bit to 1
                @param bits_entry
                @param id index of the bit in the entry to set to 1
         */
        void setBitsEntry(UINT *&bits_entry, int id);

        /**
                get a single bit content
                @param bits_entry
                @param id index of the bit in the entry
                @return TRUE if bit ID is 1, FALSE otherwise
         */
        bool getBitsEntry(UINT *&bits_entry, int id);

        /**
                get bit blocks, each block span num_state bits
                @param bit_vec bit block vector
                @param index block index
                @param bits_entry (OUT) content of the block at index
         */
        void getBitsBlock(UINT *bit_vec, int index, UINT *&bits_entry);

        /**
                set bit blocks, each block span num_state bits
                @param bit_vec (OUT) bit block vector
                @param index block index
                @param bits_entry the content of the block at index
         */
        void setBitsBlock(UINT *&bit_vec, int index, UINT *bits_entry);

        virtual void saveCurrentTree(double logl)
        {
        } // save current tree
};

#endif