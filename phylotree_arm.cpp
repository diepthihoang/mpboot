#include "phylotree.h"

void PhyloTree::changeLikelihoodKernel(LikelihoodKernel lk) {
    if (sse == lk) return;
    if ((sse == LK_EIGEN || sse == LK_EIGEN_SSE) && (lk == LK_NORMAL || lk == LK_SSE)) {
        sse = lk;
        deleteAllPartialLh();
        initializeAllPartialLh();
        clearAllPartialLH();
    } else {
        sse = lk;
    }
}

void PhyloTree::computePtnInvar() {
    size_t nptn = aln->getNPattern(), ptn;
    size_t maxptn = get_safe_upper_limit(nptn + model_factory->unobserved_ptns.size());
    int nstates = aln->num_states;

    double *state_freq = aligned_alloc_double(nstates);
    model->getStateFrequency(state_freq);
    memset(ptn_invar, 0, maxptn * sizeof(double));
    double p_invar = site_rate->getPInvar();
    if (p_invar != 0.0) {
        for (ptn = 0; ptn < nptn; ptn++) {
            if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
                ptn_invar[ptn] = p_invar * state_freq[(int)(*aln)[ptn][0]];
            }
        }
        for (ptn = 0; ptn < model_factory->unobserved_ptns.size(); ptn++) {
            ptn_invar[nptn + ptn] = p_invar * state_freq[(int)model_factory->unobserved_ptns[ptn]];
        }
    }
    aligned_free(state_freq);
}

void PhyloTree::computePartialLikelihood(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    computePartialLikelihoodNaive(dad_branch, dad);
}

double PhyloTree::computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
    return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
}

double PhyloTree::computeLikelihoodDerv(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);
}

template <const int nstates>
double PhyloTree::computeLikelihoodBranchEigen(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
    return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
}

template double PhyloTree::computeLikelihoodBranchEigen<2>(PhyloNeighbor*, PhyloNode*, double*);
template double PhyloTree::computeLikelihoodBranchEigen<4>(PhyloNeighbor*, PhyloNode*, double*);
template double PhyloTree::computeLikelihoodBranchEigen<20>(PhyloNeighbor*, PhyloNode*, double*);
