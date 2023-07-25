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

#include "tools.h"
#include "timeutil.h"

VerboseMode verbose_mode;

/*
        WIN32 does not define gettimeofday() function.
        Here declare it extra for WIN32 only.
 */
//#if defined(WIN32) && !defined(HAVE_GETTIMEOFDAY)
#if defined(WIN32)
#include <sstream>
#endif
//
//struct timezone {
//};
//
//void gettimeofday(struct timeval* t, void* timezone) {
//    struct _timeb timebuffer;
//    _ftime(&timebuffer);
//    t->tv_sec = timebuffer.time;
//    t->tv_usec = 1000 * timebuffer.millitm;
//}
//#else
//#include <sys/time.h>
//#endif


/********************************************************
        Defining DoubleMatrix methods
 ********************************************************/

/*DoubleMatrix::DoubleMatrix(int arows, int acols) {
        rows = arows;
        cols = acols;
        size =  rows * cols;
        value = new double[size];
}

void DoubleMatrix::setZero() {
        memset(value, 0, size * sizeof(double));
}


DoubleMatrix::~DoubleMatrix() {
        if (value) delete value;
        value = NULL;
}
 */

/********************************************************
        Miscellaneous
 ********************************************************/

/**
        Output an error to screen, then exit program
        @param error error message
 */
/*
void outError(char *error)
{
        cerr << "ERROR: " << error << endl;
        exit(2);
}
 */

/**
        Output an error to screen, then exit program
        @param error error message
 */
void outError(const char *error) {
    cerr << "ERROR: " << error << endl;
    exit(2);
}

/**
        Output an error to screen, then exit program
        @param error error message
 */
void outError(string error) {
    outError(error.c_str());
}

void outError(const char *error, const char *msg) {
    string str = error;
    str += msg;
    outError(str);
}

void outError(const char *error, string msg) {
    string str = error;
    str += msg;
    outError(str);
}

/**
        Output a warning message to screen
        @param error warning message
 */
void outWarning(const char *warn) {
    cerr << "WARNING: " << warn << endl;
}

void outWarning(string warn) {
    outWarning(warn.c_str());
}

double randomLen(Params &params) {
    double ran = static_cast<double> (random_int(999) + 1) / 1000;
    double len = -params.mean_len * log(ran);

    if (len < params.min_len) {
        int fac = random_int(1000);
        double delta = static_cast<double> (fac) / 1000.0; //delta < 1.0
        len = params.min_len + delta / 1000.0;
    }

    if (len > params.max_len) {
        int fac = random_int(1000);
        double delta = static_cast<double> (fac) / 1000.0; //delta < 1.0
        len = params.max_len - delta / 1000.0;
    }
    return len;
}

//From Tung

string convertIntToString(int number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

string convertDoubleToString(double number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

//From Tung

bool copyFile(const char SRC[], const char DEST[]) {
    std::ifstream src; // the source file
    std::ofstream dest; // the destination file

    src.open(SRC, std::ios::binary); // open in binary to prevent jargon at the end of the buffer
    dest.open(DEST, std::ios::binary); // same again, binary
    if (!src.is_open() || !dest.is_open())
        return false; // could not be copied

    dest << src.rdbuf(); // copy the content
    dest.close(); // close destination file
    src.close(); // close source file

    return true; // file copied successfully
}

bool fileExists(string strFilename) {
    struct stat stFileInfo;
    bool blnReturn;
    int intStat;

    // Attempt to get the file attributes
    intStat = stat(strFilename.c_str(), &stFileInfo);
    if (intStat == 0) {
        // We were able to get the file attributes
        // so the file obviously exists.
        blnReturn = true;
    } else {
        // We were not able to get the file attributes.
        // This may mean that we don't have permission to
        // access the folder which contains this file. If you
        // need to do that level of checking, lookup the
        // return values of stat which will give you
        // more details on why stat failed.
        blnReturn = false;
    }
    return (blnReturn);
}

int convert_int(const char *str) throw (string) {
    char *endptr;
    int i = strtol(str, &endptr, 10);

    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    return i;
}

int convert_int(const char *str, int &end_pos) throw (string) {
	char *endptr;
	int i = strtol(str, &endptr, 10);

	if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL) {
		string err = "Expecting integer, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}
	end_pos = endptr - str;
	return i;
}

void convert_int_vec(const char *str, IntVector &vec) throw (string) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    do {
		int i = strtol(beginptr, &endptr, 10);

		if ((i == 0 && endptr == beginptr) || abs(i) == HUGE_VALL) {
			string err = "Expecting integer, but found \"";
			err += beginptr;
			err += "\" instead";
			throw err;
		}
		vec.push_back(i);
		if (*endptr == ',') endptr++;
		beginptr = endptr;
    } while (*endptr != 0);
}


double convert_double(const char *str) throw (string) {
    char *endptr;
    double d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || *endptr != 0) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    return d;
}

double convert_double(const char *str, int &end_pos) throw (string) {
	char *endptr;
	double d = strtod(str, &endptr);
	if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF) {
		string err = "Expecting floating-point number, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}
	end_pos = endptr - str;
	return d;
}

string convert_time(const double sec) {
    int sec_int = (int) floor(sec);
    int secs = sec_int % 60;
    int mins = (sec_int % 3600) / 60;
    int hours = sec_int / 3600;
    stringstream ss;
    ss << hours << "h:" << mins << "m:" << secs << "s";
    return ss.str();
}

void convert_range(const char *str, int &lower, int &upper, int &step_size) throw (string) {
    char *endptr;
    char *beginptr = (char*) str;

    // parse the lower bound of the range
    int d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    //lower = d;
    int d_save = d;
    upper = d;
    if (*endptr == 0) return;


    // parse the upper bound of the range
    str = endptr + 1;
    d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    lower = d_save;
    upper = d;
    if (*endptr == 0) return;

    // parse the step size of the range
    str = endptr + 1;
    d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    step_size = d;
    str = beginptr;

}

void convert_range(const char *str, double &lower, double &upper, double &step_size) throw (string) {
    char *endptr;
    char *beginptr = (char*) str;

    // parse the lower bound of the range
    double d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    //lower = d;
    double d_save = d;
    upper = d;
    if (*endptr == 0) return;


    // parse the upper bound of the range
    str = endptr + 1;
    d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    lower = d_save;
    upper = d;
    if (*endptr == 0) return;

    // parse the step size of the range
    str = endptr + 1;
    d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || *endptr != 0) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    step_size = d;
    str = beginptr;

}

void readWeightFile(Params &params, int ntaxa, double &scale, StrVector &tax_name, DoubleVector &tax_weight) {
    cout << "Reading scale factor and taxa weights file " << params.param_file << " ..." << endl;
    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(params.param_file);
        string name, tmp;

        in >> tmp;
        scale = convert_double(tmp.c_str());

        for (; !in.eof() && ntaxa > 0; ntaxa--) {
            // remove the failbit
            in.exceptions(ios::badbit);
            if (!(in >> name)) break;
            // set the failbit again
            in.exceptions(ios::failbit | ios::badbit);

            tax_name.push_back(name);
            // read the sequence weight
            in >> tmp;
            tax_weight.push_back(convert_double(tmp.c_str()));
        }
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
    } catch (ios::failure) {
        outError(ERR_READ_INPUT);
    } catch (string str) {
        outError(str);
    }
}

void readStringFile(const char* filename, int max_num, StrVector &strv) {
    try {
        ifstream in;
        // set the failbit and badbit
        in.exceptions(ios::failbit | ios::badbit);
        in.open(filename);
        string name;

        // remove the failbit
        in.exceptions(ios::badbit);
        for (; !in.eof() && max_num > 0; max_num--) {
            if (!(in >> name)) break;
            strv.push_back(name);
        }
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
    } catch (ios::failure) {
        outError(ERR_READ_INPUT);
    }
}

void readInitTaxaFile(Params &params, int ntaxa, StrVector &tax_name) {
    cout << "Reading initial taxa set file " << params.initial_file << " ..." << endl;
    readStringFile(params.initial_file, ntaxa, tax_name);
}

void printString2File(string myString, string filename) {
    ofstream myfile(filename.c_str());
    if (myfile.is_open()) {
        myfile << myString;
        myfile.close();
    } else {
        cout << "Unable to open file " << filename << endl;
    }
}

void readInitAreaFile(Params &params, int nareas, StrVector &area_name) {
    cout << "Reading initial area file " << params.initial_area_file << " ..." << endl;
    readStringFile(params.initial_area_file, nareas, area_name);
}

void readAreasBoundary(char *file_name, MSetsBlock *areas, double *areas_boundary) {

    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(file_name);

        int nset;
        in >> nset;
        if (nset != areas->getNSets())
            throw "File has different number of areas";
        int pos = 0, seq1, seq2;
        for (seq1 = 0; seq1 < nset; seq1++) {
            string seq_name;
            in >> seq_name;
            if (seq_name != areas->getSet(seq1)->name)
                throw "Area name " + seq_name + " is different from " + areas->getSet(seq1)->name;
            for (seq2 = 0; seq2 < nset; seq2++) {
                in >> areas_boundary[pos++];
            }
        }
        // check for symmetric matrix
        for (seq1 = 0; seq1 < nset - 1; seq1++) {
            if (areas_boundary[seq1 * nset + seq1] <= 1e-6)
                throw "Diagonal elements of distance matrix should represent the boundary of single areas";
            for (seq2 = seq1 + 1; seq2 < nset; seq2++)
                if (areas_boundary[seq1 * nset + seq2] != areas_boundary[seq2 * nset + seq1])
                    throw "Shared boundary between " + areas->getSet(seq1)->name + " and " + areas->getSet(seq2)->name + " is not symmetric";
        }


        in.close();
        cout << "Areas relation matrix was read from " << file_name << endl;
    } catch (const char *str) {
        outError(str);
    } catch (string str) {
        outError(str);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, file_name);
    }

}

void readTaxaSets(char *filename, MSetsBlock *sets) {
    TaxaSetNameVector *allsets = sets->getSets();
    try {
        int count = 0;
        ifstream in;
        // set the failbit and badbit
        in.exceptions(ios::failbit | ios::badbit);
        in.open(filename);
        string name;

        // remove the failbit
        in.exceptions(ios::badbit);
        while (!in.eof()) {
            int ntaxa = 0;
            string str;
            if (!(in >> str)) break;
            ntaxa = convert_int(str.c_str());
            if (ntaxa <= 0) throw "Number of taxa must be > 0";
            count++;
            //allsets->resize(allsets->size()+1);
            TaxaSetName *myset = new TaxaSetName;
            allsets->push_back(myset);
            myset->name = "";
            myset->name += count;
            for (; ntaxa > 0; ntaxa--) {
                string str;
                if (!(in >> str)) throw "Cannot read in taxon name";
                if ((ntaxa > 1) && in.eof()) throw "Unexpected end of file while reading taxon names";
                myset->taxlist.push_back(str);
            }
        }
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
        if (count == 0) throw "No set found, you must specify at least 1 set";
    } catch (ios::failure) {
        outError(ERR_READ_INPUT);
    } catch (const char *str) {
        outError(str);
    } catch (string str) {
        outError(str);
    }
}

void get2RandNumb(const int size, int &first, int &second) {
    // pick a random element
    first = random_int(size);
    // pick a random element from what's left (there is one fewer to choose from)...
    second = random_int(size - 1);
    // ...and adjust second choice to take into account the first choice
    if (second >= first) {
        ++second;
    }
}

void parseArg(int argc, char *argv[], Params &params) {
    int cnt;
    verbose_mode = VB_MIN;
    params.tree_gen = NONE;
    params.user_file = NULL;
    params.out_prefix = NULL;
    params.out_file = NULL;
    params.sub_size = 0;
    params.pd_proportion = 0.0;
    params.min_proportion = 0.0;
    params.step_proportion = 0.01;
    params.min_size = 0;
    params.step_size = 1;
    params.find_all = false;
    params.run_mode = DETECTED;
    params.detected_mode = DETECTED;
    params.param_file = NULL;
    params.initial_file = NULL;
    params.initial_area_file = NULL;
    params.pdtaxa_file = NULL;
    params.areas_boundary_file = NULL;
    params.boundary_modifier = 1.0;
    params.dist_file = NULL;
    params.compute_obs_dist = false;
    params.compute_jc_dist = true;
    params.compute_ml_dist = true;
    params.compute_ml_tree = true;
    params.budget_file = NULL;
    params.overlap = 0;
    params.is_rooted = false;
    params.sample_size = -1;
    params.repeated_time = 1;
    //params.nr_output = 10000;
    params.nr_output = 0;
    //params.smode = EXHAUSTIVE;
    params.intype = IN_OTHER;
    params.budget = -1;
    params.min_budget = -1;
    params.step_budget = 1;
    params.root = NULL;
    params.num_splits = 0;
    params.min_len = 0.001;
    params.mean_len = 0.1;
    params.max_len = 0.999;
    params.num_zero_len = 0;
    params.pd_limit = 100;
    params.calc_pdgain = false;
    params.multi_tree = false;
    params.second_tree = NULL;
    params.tree_weight_file = NULL;
    params.consensus_type = CT_NONE;
    params.find_pd_min = false;
    params.branch_cluster = 0;
    params.taxa_order_file = NULL;
    params.endemic_pd = false;
    params.exclusive_pd = false;
    params.complement_area = NULL;
    params.scaling_factor = -1;
    params.numeric_precision = -1;
    params.binary_programming = false;
    params.quad_programming = false;
    params.test_input = TEST_NONE;
    params.tree_burnin = 0;
    params.tree_max_count = 1000000;
    params.split_threshold = 0.0;
    params.split_weight_threshold = -1000;
    params.split_weight_summary = SW_SUM;
    params.gurobi_format = true;
    params.gurobi_threads = 1;
    params.num_bootstrap_samples = 0;
    params.bootstrap_spec = NULL;

    params.aln_file = NULL;
    params.treeset_file = NULL;
    params.topotest_replicates = 0;
    params.do_weighted_test = false;
    params.do_au_test = false;
    params.siteLL_file = NULL; //added by MA
    params.partition_file = NULL;
    params.partition_type = 0;
    params.remove_empty_seq = true;
    params.terrace_aware = true;
    params.sequence_type = NULL;
    params.aln_output = NULL;
    params.aln_site_list = NULL;
    params.aln_output_format = ALN_PHYLIP;
    params.gap_masked_aln = NULL;
    params.concatenate_aln = NULL;
    params.aln_nogaps = false;
    params.aln_no_const_sites = false;
//    params.parsimony = false;
//    params.parsimony_tree = false;
    params.tree_spr = false;
    params.nexus_output = false;
    params.k_representative = 4;
    params.loglh_epsilon = 0.001;
    params.numSmoothTree = 1;
    params.nni5 = false; // Diep: Revert for UFBoot-MP release
    params.leastSquareBranch = false;
    params.pars_branch_length = false;
    params.bayes_branch_length = false;
    params.manuel_analytic_approx = false;
    params.leastSquareNNI = false;
    params.ls_var_type = OLS;
    params.maxCandidates = 100; // Less RAS === less runtime
    params.popSize = 5; // Diep: reset for weighted
    params.p_delete = -1;
    params.min_iterations = -1;
    params.max_iterations = 1;
    params.stop_condition = SC_UNSUCCESS_ITERATION;
    params.stop_confidence = 0.95;
    params.model_name = "";
    params.model_set = NULL;
    params.store_trans_matrix = false;
    //params.freq_type = FREQ_EMPIRICAL;
    params.freq_type = FREQ_UNKNOWN;
    params.num_rate_cats = 4;
    params.gamma_shape = -1.0;
    params.gamma_median = false;
    params.p_invar_sites = -1.0;
    params.optimize_model_rate_joint = false;
    params.optimize_by_newton = true;
    params.fixed_branch_length = false;
    params.iqp_assess_quartet = IQP_DISTANCE;
    params.iqp = false;
    params.write_intermediate_trees = 0;
    params.avoid_duplicated_trees = false;
    params.rf_dist_mode = 0;
    params.mvh_site_rate = false;
    params.rate_mh_type = true;
    params.discard_saturated_site = false;
    params.mean_rate = 1.0;
    params.aLRT_threshold = 101;
    params.aLRT_replicates = 0;
    params.localbp_replicates = 0;
    params.SSE = LK_EIGEN_SSE;
    params.print_site_lh = 0;
    params.print_site_rate = false;
    params.print_tree_lh = false;
    params.print_site_pars = 0;
    params.print_site_pars_user_tree = false;
    params.lambda = 1;
    params.speed_conf = 1.0;
    params.whtest_simulations = 1000;
    params.mcat_type = MCAT_LOG + MCAT_PATTERN;
    params.rate_file = NULL;
    params.ngs_file = NULL;
    params.ngs_mapped_reads = NULL;
    params.ngs_ignore_gaps = true;
    params.do_pars_multistate = false;
    params.gene_pvalue_file = NULL;
    params.gene_scale_factor = -1;
    params.gene_pvalue_loga = false;
    params.second_align = NULL;
    params.ncbi_taxid = 0;
    params.ncbi_taxon_level = NULL;
    params.ncbi_names_file = NULL;
    params.ncbi_ignore_level = NULL;

	params.eco_dag_file  = NULL;
	params.eco_type = NULL;
	params.eco_detail_file = NULL;
	params.k_percent = 0;
	params.diet_min = 0;
	params.diet_max = 0;
	params.diet_step = 0;
	params.eco_weighted = false;
	params.eco_run = 0;

	params.upper_bound = false;

    params.gbo_replicates = 0;
	params.ufboot_epsilon = 0.5;
    params.check_gbo_sample_size = 0;
    params.use_rell_method = true;
    params.use_elw_method = false;
    params.use_weighted_bootstrap = false;
    params.use_max_tree_per_bootstrap = true;
    params.max_candidate_trees = 0;
    params.distinct_trees = false;
    params.online_bootstrap = true;
    params.min_correlation = 0.99;
    params.step_iterations = 100;
    params.store_candidate_trees = false;
	params.print_ufboot_trees = false;
    //const double INF_NNI_CUTOFF = -1000000.0;
    params.nni_cutoff = -1000000.0;
    params.estimate_nni_cutoff = false;
    params.nni_sort = false;
    //params.nni_opt_5branches = false;
    params.testNNI = false;
    params.approximate_nni = false;
    params.do_compression = false;

    params.new_heuristic = true;
    params.write_best_trees = false;
    params.iteration_multiple = 1;
    params.initPerStrength = 0.5;
#ifdef USING_PLL
    params.pll = true;
#else
    params.pll = false;
#endif
    params.modeps = 0.001;
    params.parbran = false;
    params.binary_aln_file = NULL;
    params.maxtime = 1000000;
    params.reinsert_par = false;
    params.bestStart = true;
    params.snni = true; // turn on sNNI default now
//    params.autostop = true; // turn on auto stopping rule by default now
    params.unsuccess_iteration = -1;
    params.speednni = true; // turn on reduced hill-climbing NNI by default now
    params.adaptPert = false;
    params.numParsTrees = 100;
    params.sprDist = -1;
    params.numNNITrees = 20;
    params.avh_test = 0;
    params.bootlh_test = 0;
    params.bootlh_partitions = NULL;
    params.site_freq_file = NULL;

    params.maximum_parsimony = true;// Diep: Revert for UFBoot-MP release
    params.multiple_hits = false;
    params.store_top_boot_trees = 0;
    params.ratchet_iter = 1;// Diep: Revert for UFBoot-MP release
    params.ratchet_wgt = 1; // default if just specify -ratchet
    params.ratchet_percent = 50; // default if just specify -ratchet
    params.compute_parsimony = false;
    params.newick_to_tnt = false;
    params.newick_to_nexus = false;
    params.sankoff_cost_file = NULL;
    params.sankoff_short_int = true; // Diep: Revert for MPBoot release
    params.condense_parsimony_equiv_sites = false;
    params.spr_parsimony = true;// Diep: Revert for UFBoot-MP release
    params.spr_mintrav = 1; // same as PLL
    params.spr_maxtrav = 6; // PLL default is 20
    params.test_site_pars = false;
    params.auto_vectorize = false;
    params.sort_alignment = true;
    params.cutoff_percent = 10;// Diep: Default = 10 for UFBoot-MP release; Hidden test value: >100 for UFBoot logl-cutoff
    params.hclimb1_nni = false;
    params.no_hclimb1_bb = false;
    params.optimize_boot_trees = true;// Diep: Revert for UFBoot-MP release
    params.save_trees_off = false;
    params.minimize_iter1_candidates = false; // Diep: to go for speed
    params.cutoff_from_btrees = false;
    params.ibest_as_cand = false;
    params.opt_btree_nni = false;
    params.opt_btree_spr = 0;
    params.distinct_iter_top_boot = 0; // Diep: if not specify -distinct_iter_top_boot <int>
    params.top_boot_concensus = false;
    params.do_first_rell = false;
    params.remove_dup_seq = false;
    params.test_mode = false;
	params.savek = 0;

#ifdef _OPENMP
    params.num_threads = 0;
#endif
    params.model_test_criterion = MTC_BIC;
    params.model_test_sample_size = 0;
    params.root_state = NULL;
    params.print_bootaln = false;
	params.print_subaln = false;
	params.print_partition_info = false;
	params.print_conaln = false;
	params.count_trees = false;
	params.print_branch_lengths = false;
	params.lh_mem_save = LM_DETECT; // auto detect
	params.start_tree = STT_PLL_PARSIMONY;
	params.print_splits_file = false;
    params.ignore_identical_seqs = true;
    params.write_init_tree = false;
    params.write_local_optimal_trees = false;

	if (params.nni5) {
	    params.nni_type = NNI5;
	} else {
	    params.nni_type = NNI1;
	}

    struct timeval tv;
    struct timezone tz;
    // initialize random seed based on current time
    gettimeofday(&tv, &tz);
    //params.ran_seed = (unsigned) (tv.tv_sec+tv.tv_usec);
    params.ran_seed = (unsigned) (tv.tv_usec);

    for (cnt = 1; cnt < argc; cnt++) {
        try {

            if (strcmp(argv[cnt], "-h") == 0 || strcmp(argv[cnt], "--help") == 0) {
#ifdef IQ_TREE
//                usage_iqtree(argv, false);
				usage_mpboot(argv, false);

#else
                usage(argv, false);
#endif
                continue;
            }
			if (strcmp(argv[cnt], "-ho") == 0 || strcmp(argv[cnt], "-?") == 0) {
//				usage_iqtree(argv, false);
				usage_mpboot(argv, false);
				continue;
			}
			if (strcmp(argv[cnt], "-hh") == 0
					|| strcmp(argv[cnt], "-hhh") == 0) {
				usage(argv, true);
				continue;
			}
			if (strcmp(argv[cnt], "-v0") == 0) {
				verbose_mode = VB_QUIET;
				continue;
			}
			if (strcmp(argv[cnt], "-v") == 0 || strcmp(argv[cnt], "-v1") == 0) {
				verbose_mode = VB_MED;
				continue;
			}
			if (strcmp(argv[cnt], "-vv") == 0
					|| strcmp(argv[cnt], "-v2") == 0) {
				verbose_mode = VB_MAX;
				continue;
			}
			if (strcmp(argv[cnt], "-vvv") == 0
					|| strcmp(argv[cnt], "-v3") == 0) {
				verbose_mode = VB_DEBUG;
				continue;
			}
			if (strcmp(argv[cnt], "-k") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -k <num_taxa>";
				convert_range(argv[cnt], params.min_size, params.sub_size,
						params.step_size);
				params.k_representative = params.min_size;
				continue;
			}
			if (strcmp(argv[cnt], "-pre") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pre <output_prefix>";
				params.out_prefix = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-pp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pp <pd_proportion>";
				convert_range(argv[cnt], params.min_proportion,
						params.pd_proportion, params.step_proportion);
				if (params.pd_proportion < 0 || params.pd_proportion > 1)
					throw "PD proportion must be between 0 and 1";
				continue;
			}
			if (strcmp(argv[cnt], "-mk") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mk <min_taxa>";
				params.min_size = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bud") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bud <budget>";
				convert_range(argv[cnt], params.min_budget, params.budget,
						params.step_budget);
				continue;
			}
			if (strcmp(argv[cnt], "-mb") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mb <min_budget>";
				params.min_budget = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-o") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -o <taxon>";
				params.root = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-root") == 0) {
				params.is_rooted = true;
				continue;
			}
			if (strcmp(argv[cnt], "-all") == 0) {
				params.find_all = true;
				continue;
			}
			if (strcmp(argv[cnt], "-g") == 0
					|| strcmp(argv[cnt], "--greedy") == 0) {
				params.run_mode = GREEDY;
				continue;
			}
			if (strcmp(argv[cnt], "-pr") == 0
					|| strcmp(argv[cnt], "--pruning") == 0) {
				params.run_mode = PRUNING;
				//continue; } if (strcmp(argv[cnt],"--both") == 0) {
				//params.run_mode = BOTH_ALG;
				continue;
			}
			if (strcmp(argv[cnt], "-e") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -e <file>";
				params.param_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-if") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -if <file>";
				params.initial_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-nni_nr_step") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nni_nr_step <newton_raphson_steps>";
				NNI_MAX_NR_STEP = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-ia") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ia <file>";
				params.initial_area_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-u") == 0) {
				// file containing budget information
				cnt++;
				if (cnt >= argc)
					throw "Use -u <file>";
				params.budget_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-dd") == 0) {
				// compute distribution of PD score on random sets
				cnt++;
				if (cnt >= argc)
					throw "Use -dd <sample_size>";
				params.run_mode = PD_DISTRIBUTION;
				params.sample_size = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-ts") == 0) {
				// calculate PD score a taxa set listed in the file
				cnt++;
				//params.run_mode = PD_USER_SET;
				if (cnt >= argc)
					throw "Use -ts <taxa_file>";
				params.pdtaxa_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-bound") == 0) {
				// boundary length of areas
				cnt++;
				if (cnt >= argc)
					throw "Use -bound <file>";
				params.areas_boundary_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-blm") == 0) {
				// boundary length modifier
				cnt++;
				if (cnt >= argc)
					throw "Use -blm <boundary_modifier>";
				params.boundary_modifier = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-dist") == 0
					|| strcmp(argv[cnt], "-d") == 0) {
				// calculate distance matrix from the tree
				params.run_mode = CALC_DIST;
				cnt++;
				if (cnt >= argc)
					throw "Use -dist <distance_file>";
				params.dist_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-djc") == 0) {
				params.compute_ml_dist = false;
				continue;
			}
			if (strcmp(argv[cnt], "-dobs") == 0) {
				params.compute_obs_dist = true;
				continue;
			}
			if (strcmp(argv[cnt], "-r") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -r <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = YULE_HARDING;
				continue;
			}
			if (strcmp(argv[cnt], "-rs") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rs <alignment_file>";
				params.tree_gen = YULE_HARDING;
				params.aln_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-rstar") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rstar <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = STAR_TREE;
				continue;
			}
			if (strcmp(argv[cnt], "-ru") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ru <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = UNIFORM;
				continue;
			}
			if (strcmp(argv[cnt], "-rcat") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcat <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = CATERPILLAR;
				continue;
			}
			if (strcmp(argv[cnt], "-rbal") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rbal <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = BALANCED;
				continue;
			}
            if (strcmp(argv[cnt], "-keep_ident") == 0) {
                params.ignore_identical_seqs = false;
                continue;
            }
			if (strcmp(argv[cnt], "-rcsg") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcsg <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = CIRCULAR_SPLIT_GRAPH;
				continue;
			}
			if (strcmp(argv[cnt], "-rpam") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rpam <num_splits>";
				params.num_splits = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-rlen") == 0) {
				cnt++;
				if (cnt >= argc - 2)
					throw "Use -rlen <min_len> <mean_len> <max_len>";
				params.min_len = convert_double(argv[cnt]);
				params.mean_len = convert_double(argv[cnt + 1]);
				params.max_len = convert_double(argv[cnt + 2]);
				cnt += 2;
				continue;
			}
			if (strcmp(argv[cnt], "-rzero") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rzero <num_zero_branch>";
				params.num_zero_len = convert_int(argv[cnt]);
				if (params.num_zero_len < 0)
					throw "num_zero_len must not be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-rset") == 0) {
				cnt++;
				if (cnt >= argc - 1)
					throw "Use -rset <overlap> <outfile>";
				params.overlap = convert_int(argv[cnt]);
				cnt++;
				params.pdtaxa_file = argv[cnt];
				params.tree_gen = TAXA_SET;
				continue;
			}
			if (strcmp(argv[cnt], "-rep") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rep <repeated_times>";
				params.repeated_time = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-lim") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lim <pd_limit>";
				params.pd_limit = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-noout") == 0) {
				params.nr_output = 0;
				continue;
			}
			if (strcmp(argv[cnt], "-1out") == 0) {
				params.nr_output = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-oldout") == 0) {
				params.nr_output = 100;
				continue;
			}
			if (strcmp(argv[cnt], "-nexout") == 0) {
				params.nexus_output = true;
				continue;
			}
			if (strcmp(argv[cnt], "-exhaust") == 0) {
				params.run_mode = EXHAUSTIVE;
				continue;
			}
			if (strcmp(argv[cnt], "-seed") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -seed <random_seed>";
				params.ran_seed = (unsigned) convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-pdgain") == 0) {
				params.calc_pdgain = true;
				continue;
			}
			if (strcmp(argv[cnt], "-sup") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sup <target_tree_file>";
				params.second_tree = argv[cnt];
				params.consensus_type = CT_ASSIGN_SUPPORT;
				continue;
			}
			if (strcmp(argv[cnt], "-sup2") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sup2 <target_tree_file>";
				params.second_tree = argv[cnt];
				params.consensus_type = CT_ASSIGN_SUPPORT_EXTENDED;
				continue;
			}
			if (strcmp(argv[cnt], "-treew") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -treew <tree_weight_file>";
				params.tree_weight_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-con") == 0) {
				params.consensus_type = CT_CONSENSUS_TREE;
				continue;
			}
			if (strcmp(argv[cnt], "-net") == 0) {
				params.consensus_type = CT_CONSENSUS_NETWORK;
			} /**MINH ANH: to serve some statistics on tree*/
			else if (strcmp(argv[cnt], "-comp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -comp <treefile>";
				params.consensus_type = COMPARE;
				params.second_tree = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-stats") == 0) {
				params.run_mode = STATS;
				continue;
			}
			if (strcmp(argv[cnt], "-gbo") == 0) { //guided bootstrap
				cnt++;
				if (cnt >= argc)
					throw "Use -gbo <site likelihod file>";
				params.siteLL_file = argv[cnt];
				//params.run_mode = GBO;
			} // MA
			else if (strcmp(argv[cnt], "-mprob") == 0) { //compute multinomial distribution probability
				cnt++;
				if (cnt >= argc)
					throw "Use -mprob <ref_alignment>";
				params.second_align = argv[cnt];
				//params.run_mode = MPRO;
			} // MA
			else if (strcmp(argv[cnt], "-min") == 0) {
				params.find_pd_min = true;
				continue;
			}
			if (strcmp(argv[cnt], "-excl") == 0) {
				params.exclusive_pd = true;
				continue;
			}
			if (strcmp(argv[cnt], "-endem") == 0) {
				params.endemic_pd = true;
				continue;
			}
			if (strcmp(argv[cnt], "-compl") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -compl <area_name>";
				params.complement_area = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-cluster") == 0) {
				params.branch_cluster = 4;
				cnt++;
				if (cnt >= argc)
					throw "Use -cluster <taxa_order_file>";
				params.taxa_order_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-taxa") == 0) {
				params.run_mode = PRINT_TAXA;
				continue;
			}
			if (strcmp(argv[cnt], "-area") == 0) {
				params.run_mode = PRINT_AREA;
				continue;
			}
			if (strcmp(argv[cnt], "-scale") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -scale <scaling_factor>";
				params.scaling_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-scaleg") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -scaleg <gene_scale_factor>";
				params.gene_scale_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-scalebranch") == 0) {
				params.run_mode = SCALE_BRANCH_LEN;
				cnt++;
				if (cnt >= argc)
					throw "Use -scalebranch <scaling_factor>";
				params.scaling_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-scalenode") == 0) {
				params.run_mode = SCALE_NODE_NAME;
				cnt++;
				if (cnt >= argc)
					throw "Use -scalenode <scaling_factor>";
				params.scaling_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-prec") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -prec <numeric_precision>";
				params.numeric_precision = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-lp") == 0) {
				params.run_mode = LINEAR_PROGRAMMING;
				continue;
			}
			if (strcmp(argv[cnt], "-lpbin") == 0) {
				params.run_mode = LINEAR_PROGRAMMING;
				params.binary_programming = true;
				continue;
			}
			if (strcmp(argv[cnt], "-qp") == 0) {
				params.gurobi_format = true;
				params.quad_programming = true;
				continue;
			}
			if (strcmp(argv[cnt], "-q") == 0) {
				verbose_mode = VB_QUIET;
				continue;
			}
			if (strcmp(argv[cnt], "-mult") == 0) {
				params.multi_tree = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bi") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bi <burnin_value>";
				params.tree_burnin = convert_int(argv[cnt]);
				if (params.tree_burnin < 0)
					throw "Burnin value must not be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-tm") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -tm <tree_max_count>";
				params.tree_max_count = convert_int(argv[cnt]);
				if (params.tree_max_count < 0)
					throw "tree_max_count must not be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-t") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -t <split_threshold>";
				params.split_threshold = convert_double(argv[cnt]);
				if (params.split_threshold < 0 || params.split_threshold > 1)
					throw "Split threshold must be between 0 and 1";
				continue;
			}
			if (strcmp(argv[cnt], "-tw") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -tw <split_weight_threshold>";
				params.split_weight_threshold = convert_double(argv[cnt]);
				if (params.split_weight_threshold < 0)
					throw "Split weight threshold is negative";
				continue;
			}
			if (strcmp(argv[cnt], "-swc") == 0) {
				params.split_weight_summary = SW_COUNT;
				continue;
			}
			if (strcmp(argv[cnt], "-swa") == 0) {
				params.split_weight_summary = SW_AVG_ALL;
				continue;
			}
			if (strcmp(argv[cnt], "-swp") == 0) {
				params.split_weight_summary = SW_AVG_PRESENT;
				continue;
			}
			if (strcmp(argv[cnt], "-iwc") == 0) {
				params.test_input = TEST_WEAKLY_COMPATIBLE;
				continue;
			}
			if (strcmp(argv[cnt], "-aln") == 0
					|| strcmp(argv[cnt], "-s") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -aln, -s <alignment_file>";
				params.aln_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-z") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -aln, -z <user_trees_file>";
				params.treeset_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-zb") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -zb <#replicates>";
				params.topotest_replicates = convert_int(argv[cnt]);
				if (params.topotest_replicates < 1000)
					throw "Please specify at least 1000 replicates";
				continue;
			}
			if (strcmp(argv[cnt], "-zw") == 0) {
				params.do_weighted_test = true;
				continue;
			}
			if (strcmp(argv[cnt], "-zau") == 0) {
				params.do_au_test = true;
				continue;
			}
			if (strcmp(argv[cnt], "-sp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sp <partition_file>";
				params.partition_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-spp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -spp <type of partition model>";
				params.partition_file = argv[cnt];
				params.partition_type = 'p';
				continue;
			}
			if (strcmp(argv[cnt], "-spj") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -spj <type of partition model>";
				params.partition_file = argv[cnt];
				params.partition_type = 'j';
				continue;
			}
			if (strcmp(argv[cnt], "-keep_empty_seq") == 0) {
				params.remove_empty_seq = false;
				continue;
			}
			if (strcmp(argv[cnt], "-no_terrace") == 0) {
				params.terrace_aware = false;
				continue;
			}
			if (strcmp(argv[cnt], "-sf") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sf <ngs_file>";
				params.ngs_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-sm") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sm <ngs_mapped_read_file>";
				params.ngs_mapped_reads = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-ngs_gap") == 0) {
				params.ngs_ignore_gaps = false;
				continue;
			}
			if (strcmp(argv[cnt], "-st") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -st BIN or -st DNA or -st AA or -st CODON or -st MORPH";
				params.sequence_type = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-starttree") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -starttree BIONJ|PARS|PLLPARS";
				if (strcmp(argv[cnt], "BIONJ") == 0)
					params.start_tree = STT_BIONJ;
				else if (strcmp(argv[cnt], "PARS") == 0)
					params.start_tree = STT_PARSIMONY;
				else if (strcmp(argv[cnt], "PLLPARS") == 0)
					params.start_tree = STT_PLL_PARSIMONY;
				else
					throw "Invalid option, please use -starttree with BIONJ or PARS or PLLPARS";
				continue;
			}

			if (strcmp(argv[cnt], "-ao") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ao <alignment_file>";
				params.aln_output = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-as") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -as <aln_site_list>";
				params.aln_site_list = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-an") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -an <ref_seq_name>";
				params.ref_seq_name = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-af") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -af phy|fasta";
				if (strcmp(argv[cnt], "phy") == 0)
					params.aln_output_format = ALN_PHYLIP;
				else if (strcmp(argv[cnt], "fasta") == 0)
					params.aln_output_format = ALN_FASTA;
				else
					throw "Unknown output format";
				continue;
			}
			if (strcmp(argv[cnt], "-am") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -am <gap_masked_aln>";
				params.gap_masked_aln = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-ac") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ac <concatenate_aln>";
				params.concatenate_aln = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-nogap") == 0) {
				params.aln_nogaps = true;
				continue;
			}
			if (strcmp(argv[cnt], "-noconst") == 0) {
				params.aln_no_const_sites = true;
				continue;
			}
//			if (strcmp(argv[cnt], "-parstree") == 0) {
				// maximum parsimony
//				params.parsimony_tree = true;
//            continue; } if (strcmp(argv[cnt], "-pars") == 0) {
//                // maximum parsimony
//                params.parsimony = true;
//				continue;
//			}
			if (strcmp(argv[cnt], "-spr") == 0) {
				// subtree pruning and regrafting
				params.tree_spr = true;
				continue;
			}
			if (strcmp(argv[cnt], "-krep") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -krep <num_k>";
				params.k_representative = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-pdel") == 0
					|| strcmp(argv[cnt], "-p") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pdel <probability>";
				params.p_delete = convert_double(argv[cnt]);
				if (params.p_delete < 0.0 || params.p_delete > 1.0)
					throw "Probability of deleting a leaf must be between 0 and 1";
				continue;
			}
			if (strcmp(argv[cnt], "-pers") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pers <perturbation_strength>";
				params.initPerStrength = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-n") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -n <#iterations>";
				params.min_iterations = convert_int(argv[cnt]);
				params.stop_condition = SC_FIXED_ITERATION;
//                params.autostop = false;
				continue;
			}
			if (strcmp(argv[cnt], "-nb") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nb <#bootstrap_replicates>";
				params.min_iterations = convert_int(argv[cnt]);
				params.iqp_assess_quartet = IQP_BOOTSTRAP;
				params.avoid_duplicated_trees = true;
				continue;
			}
			if (strcmp(argv[cnt], "-mod") == 0
					|| strcmp(argv[cnt], "-m") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mod <model_name>";
				params.model_name = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-mset") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mset <model_set>";
				params.model_set = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-mh") == 0) {
				params.mvh_site_rate = true;
				params.discard_saturated_site = false;
				params.SSE = LK_NORMAL;
				continue;
			}
			if (strcmp(argv[cnt], "-mhs") == 0) {
				params.mvh_site_rate = true;
				params.discard_saturated_site = true;
				params.SSE = LK_NORMAL;
				continue;
			}
			if (strcmp(argv[cnt], "-rl") == 0) {
				params.rate_mh_type = false;
				continue;
			}
			if (strcmp(argv[cnt], "-nr") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nr <mean_rate>";
				params.mean_rate = convert_double(argv[cnt]);
				if (params.mean_rate < 0)
					throw "Wrong mean rate for MH model";
				continue;
			}
			if (strcmp(argv[cnt], "-mstore") == 0) {
				params.store_trans_matrix = true;
				continue;
			}
			if (strcmp(argv[cnt], "-nni_lh") == 0) {
				params.nni_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-lmd") == 0) {
				cnt++;
				params.lambda = convert_double(argv[cnt]);
				if (params.lambda > 1.0)
					throw "Lambda must be in (0,1]";
				continue;
			}
			if (strcmp(argv[cnt], "-nosse") == 0) {
				params.SSE = LK_NORMAL;
				continue;
			}
			if (strcmp(argv[cnt], "-slowsse") == 0) {
				params.SSE = LK_SSE;
				continue;
			}
			if (strcmp(argv[cnt], "-fastlk") == 0) {
				params.SSE = LK_EIGEN;
				continue;
			}
			if (strcmp(argv[cnt], "-fastsse") == 0
					|| strcmp(argv[cnt], "-fasttipsse") == 0) {
				params.SSE = LK_EIGEN_SSE;
				continue;
			}
			if (strcmp(argv[cnt], "-f") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -f <c | o | u | q>";
				if (strcmp(argv[cnt], "q") == 0 || strcmp(argv[cnt], "EQ") == 0)
					params.freq_type = FREQ_EQUAL;
				else if (strcmp(argv[cnt], "c") == 0
						|| strcmp(argv[cnt], "EM") == 0)
					params.freq_type = FREQ_EMPIRICAL;
				else if (strcmp(argv[cnt], "o") == 0
						|| strcmp(argv[cnt], "ES") == 0)
					params.freq_type = FREQ_ESTIMATE;
				else if (strcmp(argv[cnt], "u") == 0
						|| strcmp(argv[cnt], "UD") == 0)
					params.freq_type = FREQ_USER_DEFINED;
				else
					throw "Use -f <c | o | u | q>";
				continue;
			}
			if (strcmp(argv[cnt], "-fs") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -fs <site_freq_file>";
				params.site_freq_file = argv[cnt];
				params.SSE = LK_NORMAL;
				continue;
			}
			if (strcmp(argv[cnt], "-c") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -c <#rate_category>";
				params.num_rate_cats = convert_int(argv[cnt]);
				if (params.num_rate_cats < 1)
					throw "Wrong number of rate categories";
				continue;
			}
			if (strcmp(argv[cnt], "-a") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -a <gamma_shape>";
				params.gamma_shape = convert_double(argv[cnt]);
				if (params.gamma_shape < 0)
					throw "Wrong number of gamma shape parameter (alpha)";
				continue;
			}
			if (strcmp(argv[cnt], "-gmean") == 0) {
				params.gamma_median = false;
				continue;
			}
			if (strcmp(argv[cnt], "-gmedian") == 0) {
				params.gamma_median = true;
				continue;
			}
			if (strcmp(argv[cnt], "-i") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -i <p_invar_sites>";
				params.p_invar_sites = convert_double(argv[cnt]);
				if (params.p_invar_sites < 0)
					throw "Wrong number of proportion of invariable sites";
				continue;
			}
			if (strcmp(argv[cnt], "-brent") == 0) {
				params.optimize_by_newton = false;
				continue;
			}
			if (strcmp(argv[cnt], "-jointopt") == 0) {
				params.optimize_model_rate_joint = true;
				continue;
			}
			if (strcmp(argv[cnt], "-brent_ginvar") == 0) {
				params.optimize_model_rate_joint = false;
				continue;
			}
			if (strcmp(argv[cnt], "-fixbr") == 0) {
				params.fixed_branch_length = true;
				continue;
			}
			if (strcmp(argv[cnt], "-sr") == 0) {
				params.stop_condition = SC_WEIBULL;
				cnt++;
				if (cnt >= argc)
					throw "Use -sr <#max_iteration>";
				params.max_iterations = convert_int(argv[cnt]);
				if (params.max_iterations <= params.min_iterations)
					throw "Specified max iteration must be greater than min iteration";
				continue;
			}
			if (strcmp(argv[cnt], "-nm") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nm <#max_iteration>";
				params.max_iterations = convert_int(argv[cnt]);
				if (params.max_iterations <= params.min_iterations)
					throw "Specified max iteration must be greater than min iteration";
				continue;
			}
			if (strcmp(argv[cnt], "-sc") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sc <stop_confidence_value>";
				params.stop_confidence = convert_double(argv[cnt]);
				if (params.stop_confidence <= 0.5
						|| params.stop_confidence >= 1)
					throw "Stop confidence value must be in range (0.5,1)";
				continue;
			}
			if (strcmp(argv[cnt], "-gurobi") == 0) {
				params.gurobi_format = true;
				continue;
			}
			if (strcmp(argv[cnt], "-gthreads") == 0) {
				params.gurobi_format = true;
				cnt++;
				if (cnt >= argc)
					throw "Use -gthreads <gurobi_threads>";
				params.gurobi_threads = convert_int(argv[cnt]);
				if (params.gurobi_threads < 1)
					throw "Wrong number of threads";
				continue;
			}
			if (strcmp(argv[cnt], "-b") == 0 || strcmp(argv[cnt], "-bo") == 0) {
				params.multi_tree = true;
				if (strcmp(argv[cnt], "-bo") == 0)
					params.compute_ml_tree = false;
				if (strcmp(argv[cnt], "-b") == 0)
					params.consensus_type = CT_CONSENSUS_TREE;
				cnt++;
				if (cnt >= argc)
					throw "Use -b <num_bootstrap_samples>";
				params.num_bootstrap_samples = convert_int(argv[cnt]);
				if (params.num_bootstrap_samples < 1)
					throw "Wrong number of bootstrap samples";
				if (params.num_bootstrap_samples == 1)
					params.compute_ml_tree = false;
				if (params.num_bootstrap_samples == 1)
					params.consensus_type = CT_NONE;
				continue;
			}
			if (strcmp(argv[cnt], "-bspec") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bspec <bootstrap_specification>";
				params.bootstrap_spec = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-bc") == 0) {
				params.multi_tree = true;
				params.compute_ml_tree = false;
				cnt++;
				if (cnt >= argc)
					throw "Use -bc <num_bootstrap_samples>";
				params.num_bootstrap_samples = convert_int(argv[cnt]);
				if (params.num_bootstrap_samples < 1)
					throw "Wrong number of bootstrap samples";
				if (params.num_bootstrap_samples > 1)
					params.consensus_type = CT_CONSENSUS_TREE;
				continue;
			}
			if (strcmp(argv[cnt], "-iqppars") == 0) {
				params.iqp_assess_quartet = IQP_PARSIMONY;
				continue;
			}
			if (strcmp(argv[cnt], "-iqp") == 0) {
				params.iqp = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wlt") == 0) {
				// write all candidate trees
				params.write_local_optimal_trees = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wt") == 0) {
				params.write_intermediate_trees = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-wt2") == 0) {
				params.write_intermediate_trees = 2;
				params.avoid_duplicated_trees = true;
				params.print_tree_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wt3") == 0) {
				params.write_intermediate_trees = 3;
				params.avoid_duplicated_trees = true;
				params.print_tree_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wbl") == 0) {
				params.print_branch_lengths = true;
				continue;
			}
            if (strcmp(argv[cnt], "-wit") == 0) {
                params.write_init_tree = true;
                continue;
            }
			if (strcmp(argv[cnt], "-nodup") == 0) {
				params.avoid_duplicated_trees = true;
				continue;
			}
			if (strcmp(argv[cnt], "-rf_all") == 0) {
				params.rf_dist_mode = RF_ALL_PAIR;
				continue;
			}
			if (strcmp(argv[cnt], "-rf_adj") == 0) {
				params.rf_dist_mode = RF_ADJACENT_PAIR;
				continue;
			}
			if (strcmp(argv[cnt], "-rf") == 0) {
				params.rf_dist_mode = RF_TWO_TREE_SETS;
				cnt++;
				if (cnt >= argc)
					throw "Use -rf <second_tree>";
				params.second_tree = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-rf2") == 0) {
				params.rf_dist_mode = RF_TWO_TREE_SETS_EXTENDED;
				cnt++;
				if (cnt >= argc)
					throw "Use -rf2 <second_tree>";
				params.second_tree = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-aLRT") == 0) {
				cnt++;
				if (cnt + 1 >= argc)
					throw "Use -aLRT <threshold%> <#replicates>";
				params.aLRT_threshold = convert_int(argv[cnt]);
				if (params.aLRT_threshold < 85 || params.aLRT_threshold > 101)
					throw "aLRT threshold must be between 85 and 100";
				cnt++;
				params.aLRT_replicates = convert_int(argv[cnt]);
				if (params.aLRT_replicates < 1000
						&& params.aLRT_replicates != 0)
					throw "aLRT replicates must be at least 1000";
				continue;
			}
			if (strcmp(argv[cnt], "-alrt") == 0) {
				cnt++;
				params.aLRT_replicates = convert_int(argv[cnt]);
				if (params.aLRT_replicates < 1000
						&& params.aLRT_replicates != 0)
					throw "aLRT replicates must be at least 1000";
				continue;
			}
			if (strcmp(argv[cnt], "-lbp") == 0) {
				cnt++;
				params.localbp_replicates = convert_int(argv[cnt]);
				if (params.localbp_replicates < 1000
						&& params.localbp_replicates != 0)
					throw "Local bootstrap (LBP) replicates must be at least 1000";
				continue;
			}
			if (strcmp(argv[cnt], "-wsl") == 0) {
				params.print_site_lh = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-wspars") == 0) {
				params.print_site_pars = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-wslg") == 0) {
				params.print_site_lh = 2;
				continue;
			}
			if (strcmp(argv[cnt], "-wsr") == 0) {
				params.print_site_rate = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wba") == 0) {
				params.print_bootaln = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wsa") == 0) {
				params.print_subaln = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wtl") == 0) {
				params.print_tree_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wpi") == 0) {
				params.print_partition_info = true;
				params.print_conaln = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wca") == 0) {
				params.print_conaln = true;
				continue;
			}

			if (strcmp(argv[cnt], "-wsplits") == 0) {
				params.print_splits_file = true;
				continue;
			}
			if (strcmp(argv[cnt], "-ns") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ns <num_simulations>";
				params.whtest_simulations = convert_int(argv[cnt]);
				if (params.whtest_simulations < 1)
					throw "Wrong number of simulations for WH-test";
				continue;
			}
			if (strcmp(argv[cnt], "-mr") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mr <rate_file>";
				params.rate_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-cat_mean") == 0) {
				params.mcat_type |= MCAT_MEAN;
				continue;
			}
			if (strcmp(argv[cnt], "-cat_nolog") == 0) {
				params.mcat_type &= (127 - MCAT_LOG);
				continue;
			}
			if (strcmp(argv[cnt], "-cat_site") == 0) {
				params.mcat_type &= (127 - MCAT_PATTERN);
				continue;
			}
			if (strcmp(argv[cnt], "-tina") == 0) {
				params.do_pars_multistate = true;
				continue;
			}
			if (strcmp(argv[cnt], "-pval") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pval <gene_pvalue_file>";
				params.gene_pvalue_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-nnitest") == 0) {
				params.testNNI = true;
				continue;
			}
			if (strcmp(argv[cnt], "-anni") == 0) {
				params.approximate_nni = true;
				continue;
			}
			if (strcmp(argv[cnt], "-nnicut") == 0) {
				params.estimate_nni_cutoff = true;
				//nni_cutoff = -5.41/2;
				continue;
			}
			if (strcmp(argv[cnt], "-nnichi2") == 0) {
				params.nni_cutoff = -5.41 / 2;
				continue;
			}
			if (strcmp(argv[cnt], "-nnicutval") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nnicutval <log_diff_value>";
				params.nni_cutoff = convert_double(argv[cnt]);
				if (params.nni_cutoff >= 0)
					throw "cutoff value for -nnicutval must be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-nnisort") == 0) {
				params.nni_sort = true;
				continue;
			}
			if (strcmp(argv[cnt], "-plog") == 0) {
				params.gene_pvalue_loga = true;
				continue;
			}
			if (strcmp(argv[cnt], "-dmp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmp <ncbi_taxid>";
				params.ncbi_taxid = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-dmplevel") == 0
					|| strcmp(argv[cnt], "-dmprank") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmprank <ncbi_taxon_rank>";
				params.ncbi_taxon_level = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-dmpignore") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmpignore <ncbi_ignore_level>";
				params.ncbi_ignore_level = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-dmpname") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmpname <ncbi_names_file>";
				params.ncbi_names_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-eco") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -eco <eco_dag_file>";
				params.eco_dag_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-k%") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -k% <k in %>";
				//convert_range(argv[cnt], params.k_percent, params.sub_size, params.step_size);
				params.k_percent = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-diet") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -diet <d in %>";
				convert_range(argv[cnt], params.diet_min, params.diet_max,
						params.diet_step);
				//params.diet = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-up") == 0) {
				params.upper_bound = true;
				continue;
			}
			if (strcmp(argv[cnt], "-ecoR") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ecoR <run number>";
				params.eco_run = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bb") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bb <#replicates>";
				params.gbo_replicates = convert_int(argv[cnt]);
				params.avoid_duplicated_trees = true;
				if (params.gbo_replicates < 1000)
					throw "#replicates must be >= 1000";
				params.consensus_type = CT_CONSENSUS_TREE;
//				params.stop_condition = SC_BOOTSTRAP_CORRELATION;
				params.stop_condition = SC_UNSUCCESS_ITERATION; // Diep: because MP already has refinement
				//params.nni5Branches = true;
				continue;
			}
			if (strcmp(argv[cnt], "-beps") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -beps <epsilon>";
				params.ufboot_epsilon = convert_double(argv[cnt]);
				if (params.ufboot_epsilon <= 0.0)
					throw "Epsilon must be positive";
				continue;
			}
			if (strcmp(argv[cnt], "-wbt") == 0) {
				params.print_ufboot_trees = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bs") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bs <begin_sampling_size>";
				params.check_gbo_sample_size = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bmax") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bmax <max_candidate_trees>";
				params.max_candidate_trees = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bcor") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bcor <min_correlation>";
				params.min_correlation = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-nstep") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nstep <step_iterations>";
				params.step_iterations = convert_int(argv[cnt]);
				if (params.step_iterations < 10
						|| params.step_iterations % 2 == 1)
					throw "At least step size of 10 and even number please";
				params.min_iterations = params.step_iterations;
				continue;
			}
			if (strcmp(argv[cnt], "-boff") == 0) {
				params.online_bootstrap = false;
				continue;
			}
			if (strcmp(argv[cnt], "-nostore") == 0
					|| strcmp(argv[cnt], "-memsave") == 0) {
				params.store_candidate_trees = false;
				continue;
			}
			if (strcmp(argv[cnt], "-lhmemsave") == 0) {
				params.lh_mem_save = LM_PER_NODE;
				continue;
			}
			if (strcmp(argv[cnt], "-nolhmemsave") == 0) {
				params.lh_mem_save = LM_ALL_BRANCH;
				continue;
			}
			if (strcmp(argv[cnt], "-storetrees") == 0) {
				params.store_candidate_trees = true;
				continue;
			}
			if (strcmp(argv[cnt], "-nodiff") == 0) {
				params.distinct_trees = false;
				continue;
			}
			if (strcmp(argv[cnt], "-treediff") == 0) {
				params.distinct_trees = true;
				continue;
			}
			if (strcmp(argv[cnt], "-norell") == 0) {
				params.use_rell_method = false;
				continue;
			}
			if (strcmp(argv[cnt], "-elw") == 0) {
				params.use_elw_method = true;
				continue;
			}
			if (strcmp(argv[cnt], "-noweight") == 0) {
				params.use_weighted_bootstrap = false;
				continue;
			}
			if (strcmp(argv[cnt], "-nomore") == 0) {
				params.use_max_tree_per_bootstrap = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bweight") == 0) {
				params.use_weighted_bootstrap = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bmore") == 0) {
				params.use_max_tree_per_bootstrap = false;
				continue;
			}
			if (strcmp(argv[cnt], "-gz") == 0) {
				params.do_compression = true;
				continue;
			}
			if (strcmp(argv[cnt], "-newheu") == 0) {
				params.new_heuristic = true;
				// Enable RAxML kernel
				continue;
			}
			if (strcmp(argv[cnt], "-maxtime") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -maxtime <time_in_minutes>";
				params.maxtime = convert_double(argv[cnt]);
				params.min_iterations = 1000000;
				params.stop_condition = SC_REAL_TIME;
				continue;
			}
			if (strcmp(argv[cnt], "-numpars") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -numpars <number_of_parsimony_trees>";
				params.numParsTrees = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-toppars") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -toppars <number_of_top_parsimony_trees>";
				params.numNNITrees = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-poplim") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -poplim <max_pop_size>";
				params.maxCandidates = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-popsize") == 0
					|| strcmp(argv[cnt], "-numcand") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -numcand <number_of_candidate_trees>";
				params.popSize = convert_int(argv[cnt]);
				assert(params.popSize < params.numParsTrees);
				continue;
			}
			if (strcmp(argv[cnt], "-savek") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -savek <number_of_saving_trees>";
				params.savek = convert_int(argv[cnt]);
				params.numParsTrees = max(params.savek, params.numParsTrees);
				assert(params.savek > 0);
				continue;
			}
			if (strcmp(argv[cnt], "-beststart") == 0) {
				params.bestStart = true;
				cnt++;
				if (cnt >= argc)
					throw "Use -best_start <binary_alignment_file>";
				params.binary_aln_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-pll") == 0) {
				params.pll = true;
				continue;
			}


//			if(strcmp(argv[cnt], "-mpars") == 0){
//            	params.maximum_parsimony = true;
//            	params.nni5 = false;
//            	params.nni_type = NNI1;
//            	continue;
//            }
			if(strcmp(argv[cnt], "-mpars_off") == 0){
            	params.maximum_parsimony = false;
            	params.nni5 = true;
            	params.nni_type = NNI5;
            	continue;
            }

            if(strcmp(argv[cnt], "-mulhits") == 0){
            	params.multiple_hits = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-topboot") == 0){
            	cnt++;
                if (cnt >= argc)
                    throw "Use -topboot <number of top trees stored for each boot sample>";
            	params.store_top_boot_trees = convert_int(argv[cnt]);
            	continue;
            }
            if(strcmp(argv[cnt], "-ratchet") == 0){
            	if(params.ratchet_iter < 0) params.ratchet_iter = 1;
            	continue;
            }
            if(strcmp(argv[cnt], "-ratchet_off") == 0){
            	params.ratchet_iter = -1;
            	continue;
            }
            if(strcmp(argv[cnt], "-ratchet_iter") == 0){
            	cnt++;
                if (cnt >= argc)
                    throw "Use -ratchet_iter <number iterations between 2 ratchets>";
            	params.ratchet_iter = convert_int(argv[cnt]);
            	continue;
            }
            if(strcmp(argv[cnt], "-ratchet_wgt") == 0){
            	cnt++;
                if (cnt >= argc)
                    throw "Use -ratchet_wgt <weight to add to each resampled character>";
            	params.ratchet_wgt = convert_int(argv[cnt]);
            	continue;
            }
            if(strcmp(argv[cnt], "-ratchet_percent") == 0){
            	cnt++;
                if (cnt >= argc)
                    throw "Use -ratchet_percent <percentage of characters in subset picked for ratchet>";
            	params.ratchet_percent = convert_int(argv[cnt]);
            	continue;
            }
			if(strcmp(argv[cnt], "-comppars") == 0){
            	params.compute_parsimony = true;
            	params.nni5 = false;
            	params.nni_type = NNI1;
            	continue;
            }
			if(strcmp(argv[cnt], "-wspars-user-tree") == 0){
            	params.print_site_pars_user_tree = true;
            	params.nni5 = false;
            	params.nni_type = NNI1;
            	continue;
            }
			if(strcmp(argv[cnt], "-totnt") == 0){
            	params.newick_to_tnt = true;
            	continue;
            }
			if(strcmp(argv[cnt], "-tonex") == 0){
            	params.newick_to_nexus = true;
            	continue;
            }
			if(strcmp(argv[cnt], "-cost") == 0){
            	cnt++;
                if (cnt >= argc)
                    throw "Use -cost <sankoff cost file name>";
            	if(params.sankoff_cost_file == NULL)
            		params.sankoff_cost_file = argv[cnt];
            	continue;
            }
			if(strcmp(argv[cnt], "-short") == 0) {
                params.sankoff_short_int = true;
                continue;
            }
			if(strcmp(argv[cnt], "-short_off") == 0) {
                params.sankoff_short_int = false; // To turn off using short for pattern parsimony
                continue;
            }
			if(strcmp(argv[cnt], "-mpcondense") == 0){
            	params.condense_parsimony_equiv_sites = true;
            	continue;
            }
			if(strcmp(argv[cnt], "-spr_pars") == 0){
            	params.spr_parsimony = true;
            	params.maximum_parsimony = true;
            	continue;
            }
			if(strcmp(argv[cnt], "-nni_pars") == 0){
            	params.spr_parsimony = false;
            	params.maximum_parsimony = true;
            	continue;
            }

			if(strcmp(argv[cnt], "-spr_mintrav") == 0){
            	cnt++;
                if (cnt >= argc)
                    throw "Use -spr_mintrav <minimal SPR radius>";
            	params.spr_mintrav = convert_int(argv[cnt]);
            	continue;
            }
			if(strcmp(argv[cnt], "-spr_maxtrav") == 0 || strcmp(argv[cnt], "-spr_rad") == 0){
            	cnt++;
                if (cnt >= argc)
                    throw "Use " + string(argv[cnt]) + " <maximal SPR radius>";
            	params.spr_maxtrav = convert_int(argv[cnt]);
            	params.sprDist = params.spr_maxtrav; // Diep: hopefully this speed the pllMakeParsimonyTreeFast...
            	continue;
            }
			if(strcmp(argv[cnt], "-sitepars") == 0){
            	params.test_site_pars = true;
            	continue;
            }
			if(strcmp(argv[cnt], "-autovec") == 0){
            	params.auto_vectorize = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-keep_aln") == 0){
            	params.sort_alignment = false;
            	continue;
            }
			if(strcmp(argv[cnt], "-cand_cutoff") == 0){
            	cnt++;
                if (cnt >= argc)
                    throw "Use -cand_cutoff <integer #s from 0 to 100 for selecting top #s percentile as bootstrap candidates>";
            	params.cutoff_percent = convert_int(argv[cnt]);
            	continue;
            }
            if(strcmp(argv[cnt], "-hclimb1_nni") == 0){
            	params.hclimb1_nni = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-no_hclimb1_bb") == 0){
            	params.no_hclimb1_bb = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-opt_btree") == 0){
            	params.optimize_boot_trees = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-opt_btree_off") == 0){
            	params.optimize_boot_trees = false;
            	continue;
            }
            if(strcmp(argv[cnt], "-save_trees_off") == 0){
            	params.save_trees_off = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-min_iter1_cand") == 0){
            	params.minimize_iter1_candidates = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-cutoff_from_btrees") == 0){
            	params.cutoff_from_btrees = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-ibest_as_cand") == 0){
            	params.ibest_as_cand = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-opt_btree_nni") == 0){
            	params.opt_btree_nni = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-do_first_rell") == 0){
            	params.do_first_rell = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-test_mode") == 0){
//            	if(argc != 7)
//            		throw "These options are to compute weighted parsimony score.\nUse -s [alignment_file] -test_mode [tree_file] -cost fitch/[cost_file]";
            	// if(argc != 5)
            	// 	throw "These options are to convert newick tree string of taxa to one of id.\nUse -s [alignment_file] -test_mode [tree_file]";

            	params.test_mode = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-remove_dup_seq") == 0){
            	params.remove_dup_seq = true;
            	continue;
            }
            if(strcmp(argv[cnt], "-opt_btree_spr") == 0){
            	cnt++;
                if (cnt >= argc)
                    throw "Use -opt_btree_spr <spr_radius>";
            	params.opt_btree_spr = convert_int(argv[cnt]);
            	continue;
            }
            if(strcmp(argv[cnt], "-distinct_iter_top_boot") == 0){
            	cnt++;
                if (cnt >= argc)
                    throw "Use -distinct_iter_top_boot <# of boot trees (before refinement step) for each boot aln>";
            	params.distinct_iter_top_boot = convert_int(argv[cnt]);
            	continue;
            }
            if(strcmp(argv[cnt], "-top_boot_concensus") == 0){
            	params.top_boot_concensus = true;
            	continue;
            }
			if (strcmp(argv[cnt], "-me") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -me <model_epsilon>";
				params.modeps = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-pllmod") == 0) {
				params.pll = true;
				continue;
			}
			if (strcmp(argv[cnt], "-pars_ins") == 0) {
				params.reinsert_par = true;
				continue;
			}
			if (strcmp(argv[cnt], "-nospeednni") == 0) {
				params.speednni = false;
				continue;
			}
			if (strcmp(argv[cnt], "-adapt") == 0) {
				params.adaptPert = true;
				continue;
			}
			if (strcmp(argv[cnt], "-snni") == 0) {
				params.snni = true;
				// dont need to turn this on here
				//params.autostop = true;
				//params.speednni = true;
				// Minh: why do you turn this on? it doubles curPerStrength at some point
				//params.adaptPert = true;
				continue;
			}
			if (strcmp(argv[cnt], "-iqpnni") == 0) {
				params.snni = false;
				params.start_tree = STT_BIONJ;
//            continue; } if (strcmp(argv[cnt], "-auto") == 0) {
//            	params.autostop = true;
				continue;
			}
			if (strcmp(argv[cnt], "-stop_cond") == 0 || strcmp(argv[cnt], "-numstop") == 0) {
				if (params.stop_condition != SC_BOOTSTRAP_CORRELATION)
					params.stop_condition = SC_UNSUCCESS_ITERATION;
				cnt++;
				params.unsuccess_iteration = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-lsbran") == 0) {
				params.leastSquareBranch = true;
				continue;
			}
			if (strcmp(argv[cnt], "-manuel") == 0) {
				params.manuel_analytic_approx = true;
				continue;
			}
			if (strcmp(argv[cnt], "-parsbran") == 0) {
				params.pars_branch_length = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bayesbran") == 0) {
				params.bayes_branch_length = true;
				continue;
			}
			if (strcmp(argv[cnt], "-fivebran") == 0
					|| strcmp(argv[cnt], "-nni5") == 0) {
				params.nni5 = true;
				params.nni_type = NNI5;
				continue;
			}
			if (strcmp(argv[cnt], "-onebran") == 0
					|| strcmp(argv[cnt], "-nni1") == 0) {
				params.nni_type = NNI1;
				params.nni5 = false;
				continue;
			}
			if (strcmp(argv[cnt], "-smooth") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -smooth <num_iterations>";
				params.numSmoothTree = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-lsnni") == 0) {
				params.leastSquareNNI = true;
				continue;
			}
			if (strcmp(argv[cnt], "-lsvar") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lsvar <o|ft|fm|st|p>";
				if (strcmp(argv[cnt], "o") == 0
						|| strcmp(argv[cnt], "ols") == 0) {
					params.ls_var_type = OLS;
					continue;
				}
				if (strcmp(argv[cnt], "ft") == 0
						|| strcmp(argv[cnt], "first_taylor") == 0) {
					params.ls_var_type = WLS_FIRST_TAYLOR;
					continue;
				}
				if (strcmp(argv[cnt], "fm") == 0
						|| strcmp(argv[cnt], "fitch_margoliash") == 0) {
					params.ls_var_type = WLS_FITCH_MARGOLIASH;
					continue;
				}
				if (strcmp(argv[cnt], "st") == 0
						|| strcmp(argv[cnt], "second_taylor") == 0) {
					params.ls_var_type = WLS_SECOND_TAYLOR;
					continue;
				}
				if (strcmp(argv[cnt], "p") == 0
						|| strcmp(argv[cnt], "pauplin") == 0) {
					params.ls_var_type = WLS_PAUPLIN;
				} else {
					throw "Use -lsvar <o|ft|fm|st|p>";
				}
				continue;
			}
			if (strcmp(argv[cnt], "-eps") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -eps <log-likelihood epsilon>";
				params.loglh_epsilon = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-pb") == 0) { // Enable parsimony branch length estimation
				params.parbran = true;
				continue;
			}
			if (strcmp(argv[cnt], "-write_best_trees") == 0) {
				params.write_best_trees = true;
				continue;
			}
			if (strcmp(argv[cnt], "-x") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -x <iteration_multiple>";
				params.iteration_multiple = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-sp_iter") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sp_iter <number_iteration>";
				params.speedup_iter = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-avh") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -avh <arndt_#bootstrap>";
				params.avh_test = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bootlh") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bootlh <#replicates>";
				params.bootlh_test = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bootpart") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bootpart <part1_length,part2_length,...>";
				params.bootlh_partitions = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-AIC") == 0) {
				params.model_test_criterion = MTC_AIC;
				continue;
			}
			if (strcmp(argv[cnt], "-AICc") == 0) {
				params.model_test_criterion = MTC_AICC;
				continue;
			}
			if (strcmp(argv[cnt], "-ms") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ms <model_test_sample_size>";
				params.model_test_sample_size = convert_int(argv[cnt]);
#ifdef _OPENMP
				continue;}if (strcmp(argv[cnt], "-omp") == 0) {
				cnt++;
				if (cnt >= argc)
				throw "Use -omp <num_threads>";
				params.num_threads = convert_int(argv[cnt]);
				if (params.num_threads < 1)
				throw "At least 1 thread please";
#endif
				continue;
			}
			if (strcmp(argv[cnt], "-rootstate") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -rootstate <rootstate>";
                params.root_state = argv[cnt];
                params.SSE = LK_NORMAL;
                continue;
			}
			if (strcmp(argv[cnt], "-ct") == 0) {
            	params.count_trees = true;
            	continue;
			}
			if (strcmp(argv[cnt], "-sprdist") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sprdist <SPR distance used in parsimony search>";
				params.sprDist = convert_int(argv[cnt]);
				continue;
			}
			if (argv[cnt][0] == '-') {
                string err = "Invalid \"";
                err += argv[cnt];
                err += "\" option.";
                throw err;
            } else {
                if (params.user_file == NULL)
                    params.user_file = argv[cnt];
                else
                    params.out_file = argv[cnt];
            }
            if (params.root != NULL && params.is_rooted)
                throw "Not allowed to specify both -o <taxon> and -root";

        }
        // try
        catch (const char *str) {
            outError(str);
            //} catch (char *str) {
            //outError(str);
        } catch (string str) {
            outError(str);
        } catch (...) {
            string err = "Unknown argument \"";
            err += argv[cnt];
            err += "\"";
            outError(err);
        }

    } // for
    if (!params.user_file && !params.aln_file && !params.ngs_file && !params.ngs_mapped_reads && !params.partition_file)
#ifdef IQ_TREE
//        usage_iqtree(argv, false);
//		usage_mpboot(argv, false);
		usage(argv, false);
#else
        usage(argv, false);
#endif
    if (!params.out_prefix) {
    	if (params.eco_dag_file)
    		params.out_prefix = params.eco_dag_file;
    	else if (params.partition_file)
            params.out_prefix = params.partition_file;
        else if (params.aln_file)
            params.out_prefix = params.aln_file;
        else if (params.ngs_file)
            params.out_prefix = params.ngs_file;
        else if (params.ngs_mapped_reads)
            params.out_prefix = params.ngs_mapped_reads;
        else
            params.out_prefix = params.user_file;
    }

    // Diep: for old mpars.spr (must work with snni)
    if(params.spr_parsimony && !params.snni){
    	outError("-spr_pars must work with -snni");
    }

	// Diep: for storing multiple best trees for each boot aln
	if(params.multiple_hits && !params.maximum_parsimony){
		outError("-mulhits must work with -mpars");
	}

    // Diep: NNI is equivalent to SPR radius = 1
    if(params.sprDist == -1){
		if(params.maximum_parsimony){
			if(!params.spr_parsimony) params.sprDist = 1;
			else params.sprDist = params.spr_maxtrav;
		}else
			params.sprDist = 20;
    }

    if(params.optimize_boot_trees == false && params.save_trees_off == true){
    	outError("-save_trees_off must work with -opt_btree");
    }else if(params.optimize_boot_trees == true && params.save_trees_off == true){
    	params.stop_condition = SC_UNSUCCESS_ITERATION;
    }
}

extern void printCopyright(ostream &out);
extern void printCopyrightMP(ostream &out);

void usage(char* argv[], bool full_command) {
//    printCopyright(cout);
	printCopyrightMP(cout); // to print UFBoot-MP info

	cout << "Minimal command-line examples (replace 'mpboot ...' with actual path to executable):" << endl << endl;

 	cout << "1. Reconstruct maximum parsimony tree from a sequence alignment" << endl
 		<< "   (example.phy):" << endl
		<< "     mpboot -s example.phy" << endl << endl;

	cout << "2. Reconstruct MP tree and assess branch supports with the MPBoot method" << endl
		<< "   (1000 replicates):" << endl
		<< "     mpboot -s example.phy -bb 1000" << endl << endl;

	cout << "To show all available options: run 'mpboot -h'" << endl << endl;

	cout << "Have a look at the tutorial and manual for more information:" << endl
		<< "     http://www.cibiv.at/software/mpboot" << endl << endl;

//    cout << "Usage: " << argv[0] << " [OPTIONS] <file_name> [<output_file>]" << endl;
//    cout << "GENERAL OPTIONS:" << endl;
//    cout << "  -hh               Print this help dialog" << endl;
//    cout << "  -h                Print help options for phylogenetic inference" << endl;
//    cout << "  <file_name>       User tree in NEWICK format or split network in NEXUS format" << endl;
//    cout << "  <output_file>     Output file to store results, default is '<file_name>.pda'" << endl;
//    cout << "  -k <num_taxa>     Find optimal set of size <num_taxa>" << endl;
//    cout << "  -k <min>:<max>    Find optimal sets of size from <min> to <max>" << endl;
//    cout << "  -k <min>:<max>:<step>" << endl;
//    cout << "                    Find optimal sets of size min, min+step, min+2*step,..." << endl;
//    cout << "  -o <taxon>        Root name to compute rooted PD (default: unrooted)" << endl;
//    cout << "  -if <file>        File containing taxa to be included into optimal sets" << endl;
//    cout << "  -e <file>         File containing branch/split scale and taxa weights" << endl;
//    cout << "  -all              Identify all multiple optimal sets" << endl;
//    cout << "  -lim <max_limit>  The maximum number of optimal sets for each k if -a is specified" << endl;
//    cout << "  -min              Compute minimal sets (default: maximal)" << endl;
//    cout << "  -1out             Print taxa sets and scores to separate files" << endl;
//    cout << "  -oldout           Print output compatible with version 0.3" << endl;
//    cout << "  -v                Verbose mode" << endl;
//    cout << endl;
//    cout << "OPTIONS FOR PHYLOGENETIC DIVERSITY (PD):" << endl;
//    cout << "  -root             Make the tree ROOTED, default is unrooted" << endl;
//    cout << "    NOTE: this option and -o <taxon> cannot be both specified" << endl;
//    cout << "  -g                Run greedy algorithm only (default: auto)" << endl;
//    cout << "  -pr               Run pruning algorithm only (default: auto)" << endl;
//    cout << endl;
//    /*
//    cout << "OPTIONS FOR SPLIT DIVERSITY:" << endl;
//    cout << "  -exhaust          Force to use exhaustive search" << endl;
//    cout << "    NOTE: by default, the program applies dynamic programming algorithm" << endl;
//    cout << "          on circular networks and exhaustive search on general networks" << endl;
//    cout << endl;*/
//    cout << "OPTIONS FOR BUDGET CONSTRAINTS:" << endl;
//    cout << "  -u <file>         File containing total budget and taxa preservation costs" << endl;
//    cout << "  -b <budget>       Total budget to conserve taxa" << endl;
//    cout << "  -b <min>:<max>    Find all sets with budget from <min> to <max>" << endl;
//    cout << "  -b <min>:<max>:<step>" << endl;
//    cout << "                    Find optimal sets with budget min, min+step, min+2*step,..." << endl;
//    cout << endl;
//    cout << "OPTIONS FOR AREA ANALYSIS:" << endl;
//    cout << "  -ts <taxa_file>   Compute/maximize PD/SD of areas (combine with -k to maximize)" << endl;
//    cout << "  -excl             Compute exclusive PD/SD" << endl;
//    cout << "  -endem            Compute endemic PD/SD" << endl;
//    cout << "  -compl <areas>    Compute complementary PD/SD given the listed <areas>" << endl;
//    cout << endl;
//
//    cout << "OPTIONS FOR VIABILITY CONSTRAINTS:" << endl;
//    cout << "  -eco <food_web>   File containing food web matrix" << endl;
//    cout << "  -k% <n>           Find optimal set of size relative the total number of taxa" << endl;
//    cout << "  -diet <min_diet>  Minimum diet portion (%) to be preserved for each predator" << endl;
//    cout << endl;
//    //if (!full_command) exit(0);
//
//    cout << "MISCELLANEOUS:" << endl;
//    cout << "  -dd <sample_size> Compute PD distribution of random sets of size k" << endl;
//    /*
//    cout << "  -gbo <sitelh_file> Compute and output the alignment of (normalized)" << endl;
//    cout << "                    expected frequencies given in site_ll_file" << endl;
//	*/
//
//    //	cout << "  -rep <times>        Repeat algorithm a number of times." << endl;
//    //	cout << "  -noout              Print no output file." << endl;
//    cout << endl;
//    //cout << "HIDDEN OPTIONS: see the source code file pda.cpp::parseArg()" << endl;

    exit(0);
}

void usage_iqtree(char* argv[], bool full_command) {
    printCopyright(cout);
    cout << "Usage: " << argv[0] << " -s <alignment> [OPTIONS] [<treefile>] " << endl << endl;
    cout << "GENERAL OPTIONS:" << endl
            << "  -?                   Printing this help dialog" << endl
            << "  -s <alignment>       Input alignment (REQUIRED) in PHYLIP/FASTA/NEXUS format" << endl
            << "  -st <data_type>      BIN, DNA, AA, CODON, or MORPH (default: auto-detect)" << endl
            << "  -sp <partition_file> Partition model specification in NEXUS format." << endl
            << "                       For single model use the -m option (see below)" << endl
            << "  -z <trees_file>      Compute log-likelihoods for all trees in the given file" << endl
            << "  <treefile>           Initial tree for tree reconstruction (default: MP)" << endl
            << "  -o <outgroup_taxon>  Outgroup taxon name for writing .treefile" << endl
            << "  -pre <PREFIX>        Using <PREFIX> for output files (default: alignment name)" << endl
#ifdef _OPENMP
            << "  -omp <#cpu_cores>    Number of cores/threads to use (default: all cores)" << endl
#endif
            << "  -seed <number>       Random seed number, normally used for debugging purpose" << endl
            << "  -v, -vv, -vvv        Verbose mode, printing more messages to screen" << endl
            << endl << "NEW STOCHASTIC TREE SEARCH ALGORITHM:" << endl
            << "  -pll                 Use phylogenetic likelihood library (PLL) (default: off)" << endl
            << "  -numpars <number>    Number of initial parsimony trees (default: 100)" << endl
            << "  -toppars <number>    Number of best parsimony trees (default: 20)" << endl
            << "  -numcand <number>    Size of the candidate tree set (defaut: 5)" << endl
            << "  -pers <perturbation> Perturbation strength for randomized NNI (default: 0.5)" << endl
            << "  -numstop <number>    Number of unsuccessful iterations to stop (default: 100)" << endl
            << "  -n <#iterations>     Fix number of iterations to <#iterations> (default: auto)" << endl
            << "  -iqp                 Use the IQP tree perturbation (default: randomized NNI)" << endl
            << "  -iqpnni              Switch back to the old IQPNNI tree search algorithm" << endl
            << endl << "ULTRAFAST BOOTSTRAP:" << endl
            << "  -bb <#replicates>    Ultrafast bootstrap (>=1000)" << endl
//            << "  -n <#iterations>     Minimum number of iterations (default: 100)" << endl
            << "  -nm <#iterations>    Maximum number of iterations (default: 1000)" << endl
			<< "  -nstep <#iterations> #Iterations for UFBoot stopping rule (default: 100)" << endl
            << "  -bcor <min_corr>     Minimum correlation coefficient (default: 0.99)" << endl
			<< "  -beps <epsilon>      RELL epsilon to break tie (default: 0.5)" << endl
            << endl << "STANDARD NON-PARAMETRIC BOOTSTRAP:" << endl
            << "  -b <#replicates>     Bootstrap + ML tree + consensus tree (>=100)" << endl
            << "  -bc <#replicates>    Bootstrap + consensus tree" << endl
            << "  -bo <#replicates>    Bootstrap only" << endl
            << "  -t <threshold>       Minimum bootstrap support [0...1) for consensus tree" << endl
            << endl << "SINGLE BRANCH TEST:" << endl
            << "  -alrt <#replicates>  SH-like approximate likelihood ratio test (SH-aLRT)" << endl
            << "  -lbp <#replicates>   Fast local bootstrap probabilities" << endl
            << endl << "SUBSTITUTION MODEL:" << endl
            << "  -m <model_name>" << endl
            << "                  DNA: HKY (default), JC, F81, K2P, K3P, K81uf, TN/TrN, TNef," << endl
            << "                       TIM, TIMef, TVM, TVMef, SYM, GTR, or 6-digit model" << endl
            << "                       specification (e.g., 010010 = HKY)" << endl
            << "              Protein: WAG (default), Poisson, cpREV, mtREV, Dayhoff, mtMAM," << endl
            << "                       JTT, LG, mtART, mtZOA, VT, rtREV, DCMut, PMB, HIVb," << endl
            << "                       HIVw, JTTDCMut, FLU, Blosum62" << endl
            << "               Binary: JC2 (default), GTR2" << endl
            << "                Codon: GY (default), ECM, MG" << endl
            << "       Morphology/SNP: MK (default), ORDERED" << endl
            << "      Model selection: TEST or TESTONLY to auto-select the best-fit model." << endl
            << "                       TESTONLY will stop the run after model selection" << endl
            << "            Otherwise: Name of file containing user-model parameters" << endl
            << "                       (rate parameters and state frequencies)" << endl
            << "  -m <model_name>+F or +FO or +FU or +FQ (default: auto)" << endl
            << "                       counted, optimized, user-defined, equal state frequency" << endl
            << "  -m <model_name>+F1x4 or +F3x4 or +F3x4C" << endl
            << "                       Codon frequencies" << endl
            << "  -m <model_name>+ASC  Ascertainment bias correction for morphological/SNP data" << endl
            << endl << "RATE HETEROGENEITY:" << endl
            << "  -m <model_name>+I or +G[n] or +I+G[n]" << endl
            << "                       Invar, Gamma, or Invar+Gamma rates. 'n' is number of" << endl
            << "                       categories for Gamma rates (default: n=4)" << endl
            << "  -a <Gamma_shape>     Gamma shape parameter for site rates (default: estimate)" << endl
            << "  -gmedian             Computing mean for Gamma rate category (default: mean)" << endl
            << "  -i <p_invar>         Proportion of invariable sites (default: estimate)" << endl
            << "  -mh                  Computing site-specific rates to .mhrate file using" << endl
            << "                       Meyer & von Haeseler (2003) method" << endl
            //<< "  -c <#categories>     Number of Gamma rate categories (default: 4)" << endl
            << endl << "TEST OF MODEL HOMOGENEITY:" << endl
            << "  -m WHTEST            Testing model (GTR+G) homogeneity assumption using" << endl
            << "                       Weiss & von Haeseler (2003) method" << endl
            << "  -ns <#simulations>   #Simulations to obtain null-distribution (default: 1000)" << endl
//            << endl << "TREE INFERENCE:" << endl
//            << "  -p <probability>     IQP: Probability of deleting leaves (default: auto)" << endl
//            << "  -k <#representative> IQP: Size of representative leaf set (default: 4)" << endl
//            << "  -n <#iterations>     Number of iterations  (default: auto)" << endl
//            << "  -sr <#iterations>    Stopping rule with max. #iterations (default: off)" << endl
//            << "  -sc <confidence>     Confidence value for stopping rule (default: 0.95)" << endl
//            << "  -spc <level>         Confidence level for NNI adaptive search (default 0.95)" << endl
//            << "  -sp_iter <number>    #iterations before NNI adaptive heuristic is started" << endl
//            << "  -lmd <lambda>        lambda parameter for the PhyML search (default 0.75)" << endl
//            << "  -nosse               Disable SSE instructions" << endl
//            << "  -wt                  Writing all intermediate trees into .treels file" << endl
//            << "  -d <file>            Reading genetic distances from file (default: JC)" << endl
//            << "  -fixbr               Fix branch lengths of <treefile>" << endl
//            << "  -seed <number>       Random seed number, normally used for debugging purpose" << endl
//            << "  -v, -vv, -vvv        Verbose mode, printing more messages to screen" << endl
            << endl << "CONSENSUS RECONSTRUCTION:" << endl
            << "  <tree_file>          Set of input trees for consensus reconstruction" << endl
            << "  -t <threshold>       Min split support in range [0,1]. 0.5 for majority-rule" << endl
            << "                       consensus (default: 0, i.e. extended consensus)" << endl
            << "  -bi <burnin>         Discarding <burnin> trees at beginning of <treefile>" << endl
            << "  -con                 Computing consensus tree to .contree file" << endl
            << "  -net                 Computing consensus network to .nex file" << endl
            << "  -sup <target_tree>   Assigning support values for <target_tree> to .suptree" << endl
            << endl << "ROBINSON-FOULDS DISTANCE:" << endl
            << "  -rf_all              Computing all-to-all RF distances of trees in <treefile>" << endl
            << "  -rf <treefile2>      Computing all RF distances between two sets of trees" << endl
            << "                       stored in <treefile> and <treefile2>" << endl
            << "  -rf_adj              Computing RF distances of adjacent trees in <treefile>" << endl
            << endl << "TREE TOPOLOGY TEST:" << endl
            << "  -zb <#replicates>    BP,KH,SH,ELW tests with RELL for trees passed via -z" << endl
            << "  -zw                  Also performing weighted-KH and weighted-SH tests" << endl
            << endl;

			cout << "GENERATING RANDOM TREES:" << endl;
			cout << "  -r <num_taxa>        Create a random tree under Yule-Harding model." << endl;
			cout << "  -ru <num_taxa>       Create a random tree under Uniform model." << endl;
			cout << "  -rcat <num_taxa>     Create a random caterpillar tree." << endl;
			cout << "  -rbal <num_taxa>     Create a random balanced tree." << endl;
			cout << "  -rcsg <num_taxa>     Create a random circular split network." << endl;
			cout << "  -rlen <min_len> <mean_len> <max_len>  " << endl;
			cout << "                       min, mean, and max branch lengths of random trees." << endl;

			cout << endl << "MISCELLANEOUS:" << endl
		    << "  -wt                  Write locally optimal trees into .treels file" << endl
			<< "  -fixbr               Fix branch lengths of <treefile>." << endl
            << "                       Used with -n 0 to compute log-likelihood of <treefile>" << endl
			<< "  -wsl                 Writing site log-likelihoods to .sitelh file" << endl
            << "  -wslg                Writing site log-likelihoods per Gamma category" << endl;
//            << "  -d <file>            Reading genetic distances from file (default: JC)" << endl
//			<< "  -d <outfile>         Calculate the distance matrix inferred from tree" << endl
//			<< "  -stats <outfile>     Output some statistics about branch lengths" << endl
//			<< "  -comp <treefile>     Compare tree with each in the input trees" << endl;


			cout << endl;

    if (full_command) {
        //TODO Print other options here (to be added)
    }

    exit(0);
}

void usage_mpboot(char* argv[], bool full_command) {
	printCopyrightMP(cout);

    cout << "Usage: " << argv[0] << " -s <alignment> [OPTIONS] [<treefile>] " << endl << endl;
    cout << "GENERAL OPTIONS:" << endl
            << "  -?                   Printing this help dialog" << endl
            << "  -s <alignment>       Input alignment (REQUIRED) in PHYLIP/FASTA/NEXUS format" << endl
            << "  -st <data_type>      BIN, DNA, AA, CODON, or MORPH (default: auto-detect)" << endl
            << "  <treefile>           Initial tree for tree reconstruction (default: MP)" << endl
            << "  -pre <PREFIX>        Using <PREFIX> for output files (default: alignment name)" << endl
            << "  -seed <number>       Random seed number, normally used for debugging purpose" << endl
            << "  -v, -vv, -vvv        Verbose mode, printing more messages to screen" << endl

			<< endl << "MPBOOT - MAXIMUM PARSIMONY BOOTSTRAP APPROXIMATION:" << endl
			<< "  -mulhits                  Store multiple equally parsimonious trees per bootstrap replicate" << endl
			<< "  -ratchet_iter <number>    Number of non-ratchet iterations before each ratchet iteration (default: 1)" << endl
			<< "  -ratchet_wgt <number>     Weight to add to each site selected for perturbation during ratchet (default: 1)" << endl
			<< "  -ratchet_percent <number> Percentage of informative sites selected for perturbation during ratchet (default: 50)" << endl
			<< "  -ratchet_off              Turn of ratchet, i.e. Only use tree perturbation" << endl
			<< "  -spr_rad <number>         Maximum radius of SPR (default: 3)" << endl
			<< "  -cand_cutoff <#s>         Use top #s percentile as cutoff for selecting bootstrap candidates (default: 10)" << endl
			<< "  -opt_btree_off            Turn off refinement step on the final bootstrap tree set" << endl
			<< "  -nni_pars                 Hill-climb by NNI instead of SPR" << endl
			<< "  -cost <file>              Read <file> for the matrix of transition cost between character states" << endl
			<< "                            Replace <file> by letter e for uniform cost." << endl

            << endl << "NEW STOCHASTIC TREE SEARCH ALGORITHM:" << endl
            << "  -numpars <number>    Number of initial parsimony trees (default: 100)" << endl
            << "  -toppars <number>    Number of best parsimony trees (default: 20)" << endl
            << "  -numcand <number>    Size of the candidate tree set (defaut: 5)" << endl
            << "  -pers <perturbation> Perturbation strength for randomized NNI (default: 0.5)" << endl
            << "  -numstop <number>    Number of unsuccessful iterations to stop (default: 100)" << endl
            << "  -n <#iterations>     Fix number of iterations to <#iterations> (default: auto)" << endl
            << "  -iqp                 Use the IQP tree perturbation (default: randomized NNI)" << endl
            << "  -iqpnni              Switch back to the old IQPNNI tree search algorithm" << endl
            << endl << "ULTRAFAST BOOTSTRAP:" << endl
            << "  -bb <#replicates>    Ultrafast bootstrap (>=1000)" << endl
//            << "  -n <#iterations>     Minimum number of iterations (default: 100)" << endl
            << "  -nm <#iterations>    Maximum number of iterations (default: 1000)" << endl
			<< "  -nstep <#iterations> #Iterations for UFBoot stopping rule (default: 100)" << endl
            << "  -bcor <min_corr>     Minimum correlation coefficient (default: 0.99)" << endl
			<< "  -beps <epsilon>      RELL epsilon to break tie (default: 0.5)" << endl
            << endl << "CONSENSUS RECONSTRUCTION:" << endl
            << "  <tree_file>          Set of input trees for consensus reconstruction" << endl
            << "  -t <threshold>       Min split support in range [0,1]. 0.5 for majority-rule" << endl
            << "                       consensus (default: 0, i.e. extended consensus)" << endl
            << "  -bi <burnin>         Discarding <burnin> trees at beginning of <treefile>" << endl
            << "  -con                 Computing consensus tree to .contree file" << endl
            << "  -net                 Computing consensus network to .nex file" << endl
            << "  -sup <target_tree>   Assigning support values for <target_tree> to .suptree" << endl
            << endl << "ROBINSON-FOULDS DISTANCE:" << endl
            << "  -rf_all              Computing all-to-all RF distances of trees in <treefile>" << endl
            << "  -rf <treefile2>      Computing all RF distances between two sets of trees" << endl
            << "                       stored in <treefile> and <treefile2>" << endl
            << "  -rf_adj              Computing RF distances of adjacent trees in <treefile>" << endl
            << endl;

			cout << "GENERATING RANDOM TREES:" << endl;
			cout << "  -r <num_taxa>        Create a random tree under Yule-Harding model." << endl;
			cout << "  -ru <num_taxa>       Create a random tree under Uniform model." << endl;
			cout << "  -rcat <num_taxa>     Create a random caterpillar tree." << endl;
			cout << "  -rbal <num_taxa>     Create a random balanced tree." << endl;
			cout << "  -rcsg <num_taxa>     Create a random circular split network." << endl;
			cout << "  -rlen <min_len> <mean_len> <max_len>  " << endl;
			cout << "                       min, mean, and max branch lengths of random trees." << endl
				<< endl;

			cout << "PRINTING SITE PARSIMONY SCORES:" << endl;
			cout << "  -wspars              When using together with parsimony tree inference, print site parsimony scores of the best tree found." << endl;
            cout << "  -wspars-user-tree <treefile> Print site parsimony scores of the user tree in <treefile>" << endl
            << endl;

			cout << endl;

    if (full_command) {
        //TODO Print other options here (to be added)
    }

    exit(0);
}

InputType detectInputFile(char *input_file) {

    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(input_file);

        unsigned char ch;
        int count = 0;
        do {
            in >> ch;
        } while (ch <= 32 && !in.eof() && count++ < 20);
        in.close();
        switch (ch) {
            case '#': return IN_NEXUS;
            case '(': return IN_NEWICK;
            case '[': return IN_NEWICK;
            case '>': return IN_FASTA;
            default:
                if (isdigit(ch)) return IN_PHYLIP;
                return IN_OTHER;
        }
    } catch (ios::failure) {
        outError("Cannot read file ", input_file);
    }
    return IN_OTHER;
}

bool overwriteFile(char *filename) {
    ifstream infile(filename);
    if (infile.is_open()) {
        cout << "Overwrite " << filename << " (y/n)? ";
        char ch;
        cin >> ch;
        if (ch != 'Y' && ch != 'y') {
            infile.close();
            return false;
        }
    }
    infile.close();
    return true;
}

void parseAreaName(char *area_names, set<string> &areas) {
    string all = area_names;
    int pos;
    while (!all.empty()) {
        pos = all.find(',');
        if (pos < 0) pos = all.length();
        areas.insert(all.substr(0, pos));
        if (pos >= all.length())
            all = "";
        else
            all = all.substr(pos + 1);
    }
}

double logFac(const int num) {
    if (num < 0) return -1.0;
    if (num == 0) return 0.0;
    double ret = 0;
    for (int i = 1; i <= num; i++)
        ret += log((double) i);
    return ret;
}

template <typename I>
I random_element(I begin, I end)
{
    const unsigned long n = std::distance(begin, end);
    const unsigned long divisor = (RAND_MAX + 1) / n;

    unsigned long k;
    do { k = std::rand() / divisor; } while (k >= n);

    return std::advance(begin, k);
}

template <class T>
inline T quantile(const vector<T>& v, const double q) {
    unsigned int size = v.size();
    if (q <= 0) return *std::min_element(v.begin(), v.end());
    if (q >= 1) return *std::max_element(v.begin(), v.end());
    //double pos = (size - 1) * q;
    //unsigned int ind = (unsigned int)(pos);
    //double delta = pos - ind;
    vector<T> w(size);
    std::copy(v, v.begin() + size, w.begin());
}

#define RAN_STANDARD 1
#define RAN_SPRNG    2
#define RAN_RAND4    3

#define RAN_TYPE 2

#if RAN_TYPE == RAN_STANDARD

int init_random(int seed) {
    srand(seed);
    cout << "(Using rand() - Standard Random Number Generator)" << endl;
    return seed;
}

int finish_random() {
	return 0;
}


#elif RAN_TYPE == RAN_RAND4
/******************************************************************************/
/* random numbers generator  (Numerical recipes)                              */
/******************************************************************************/

/* variable */
long _idum;

/* definitions */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double randomunitintervall()
/* Long period (> 2e18) random number generator. Returns a uniform random
   deviate between 0.0 and 1.0 (exclusive of endpoint values).

   Source:
   Press et al., "Numerical recipes in C", Cambridge University Press, 1992
   (chapter 7 "Random numbers", ran2 random number generator) */ {
    int j;
    long k;
    static long _idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (_idum <= 0) {
        if (-(_idum) < 1)
            _idum = 1;
        else
            _idum = -(_idum);
        _idum2 = (_idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (_idum) / IQ1;
            _idum = IA1 * (_idum - k * IQ1) - k*IR1;
            if (_idum < 0)
                _idum += IM1;
            if (j < NTAB)
                iv[j] = _idum;
        }
        iy = iv[0];
    }
    k = (_idum) / IQ1;
    _idum = IA1 * (_idum - k * IQ1) - k*IR1;
    if (_idum < 0)
        _idum += IM1;
    k = _idum2 / IQ2;
    _idum2 = IA2 * (_idum2 - k * IQ2) - k*IR2;
    if (_idum2 < 0)
        _idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - _idum2;
    iv[j] = _idum;
    if (iy < 1)
        iy += IMM1;
    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
} /* randomunitintervall */

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

int init_random(int seed) /* RAND4 */ {
    //    srand((unsigned) time(NULL));
    //    if (seed < 0)
    // 	seed = rand();
    _idum = -(long) seed;
#ifndef PARALLEL
    cout << "(Using RAND4 Random Number Generator)" << endl;
#else /* PARALLEL */
    {
        int n;
        if (PP_IamMaster) {
            cout << "(Using RAND4 Random Number Generator with leapfrog method)" << endl;
        }
        for (n = 0; n < PP_Myid; n++)
            (void) randomunitintervall();
        if (verbose_mode >= VB_MED) {
            cout << "(" << PP_Myid << ") !!! random seed set to " << seed << ", " << n << " drawn !!!" << endl;
        }
    }
#endif
    return (seed);
} /* initrandom */

int finish_random() {
	return 0;
}
/******************/

#else /* SPRNG */

/******************/

int *randstream;

int init_random(int seed) {
    //    srand((unsigned) time(NULL));
    if (seed < 0)
        seed = make_sprng_seed();
#ifndef PARALLEL
    cout << "(Using SPRNG - Scalable Parallel Random Number Generator)" << endl;
    randstream = init_sprng(0, 1, seed, SPRNG_DEFAULT); /*init stream*/
    if (verbose_mode >= VB_MED) {
        print_sprng(randstream);
    }
#else /* PARALLEL */
    if (PP_IamMaster) {
        cout << "(Using SPRNG - Scalable Parallel Random Number Generator)" << endl;
    }
    /* MPI_Bcast(&seed, 1, MPI_UNSIGNED, PP_MyMaster, MPI_COMM_WORLD); */
    randstream = init_sprng(PP_Myid, PP_NumProcs, seed, SPRNG_DEFAULT); /*initialize stream*/
    if (verbose_mode >= VB_MED) {
        cout << "(" << PP_Myid << ") !!! random seed set to " << seed << " !!!" << endl;
        print_sprng(randstream);
    }
#endif /* PARALLEL */
    return (seed);
} /* initrandom */

int finish_random() {
	return free_sprng(randstream);
}

#endif /* USE_SPRNG */

/******************/

/* returns a random integer in the range [0; n - 1] */
int random_int(int n) {
    return (int) floor(random_double() * n);
} /* randominteger */

//int randint(int a, int b) {
//	return a + (RAND_MAX * rand() + rand()) % (b + 1 - a);
//}
//

double random_double() {
#ifndef FIXEDINTRAND
#ifndef PARALLEL
#if RAN_TYPE == RAN_STANDARD
    return ((double) rand()) / ((double) RAND_MAX + 1);
#elif RAN_TYPE == RAN_SPRNG
    return sprng(randstream);
#else /* NO_SPRNG */
    return randomunitintervall();
#endif /* NO_SPRNG */
#else /* NOT PARALLEL */
#if RAN_TYPE == RAN_SPRNG
    return sprng(randstream);
#else /* NO_SPRNG */
    int m;
    for (m = 1; m < PP_NumProcs; m++)
        (void) randomunitintervall();
    PP_randn += (m - 1);
    PP_rand++;
    return randomunitintervall();
#endif /* NO_SPRNG */
#endif /* NOT PARALLEL */
#else /* FIXEDINTRAND */
    cerr << "!!! fixed \"random\" integers for testing purposes !!!" << endl;
    return 0.0;
#endif /* FIXEDINTRAND */

}

/* Following part is taken from ModelTest software */
#define	BIGX            20.0                                 /* max value to represent exp (x) */
#define	LOG_SQRT_PI     0.5723649429247000870717135          /* log (sqrt (pi)) */
#define	I_SQRT_PI       0.5641895835477562869480795          /* 1 / sqrt (pi) */
#define	Z_MAX           6.0                                  /* maximum meaningful z value */
#define	ex(x)           (((x) < -BIGX) ? 0.0 : exp (x))

/************** Normalz: probability of normal z value *********************/

/*
ALGORITHM:	Adapted from a polynomial approximation in:
                        Ibbetson D, Algorithm 209
                        Collected Algorithms of the CACM 1963 p. 616
                Note:
                        This routine has six digit accuracy, so it is only useful for absolute
                        z values < 6.  For z values >= to 6.0, Normalz() returns 0.0.
 */

double Normalz(double z) /*VAR returns cumulative probability from -oo to z VAR normal z value */ {
    double y, x, w;

    if (z == 0.0)
        x = 0.0;
    else {
        y = 0.5 * fabs(z);
        if (y >= (Z_MAX * 0.5))
            x = 1.0;
        else if (y < 1.0) {
            w = y*y;
            x = ((((((((0.000124818987 * w
                    - 0.001075204047) * w + 0.005198775019) * w
                    - 0.019198292004) * w + 0.059054035642) * w
                    - 0.151968751364) * w + 0.319152932694) * w
                    - 0.531923007300) * w + 0.797884560593) * y * 2.0;
        } else {
            y -= 2.0;
            x = (((((((((((((-0.000045255659 * y
                    + 0.000152529290) * y - 0.000019538132) * y
                    - 0.000676904986) * y + 0.001390604284) * y
                    - 0.000794620820) * y - 0.002034254874) * y
                    + 0.006549791214) * y - 0.010557625006) * y
                    + 0.011630447319) * y - 0.009279453341) * y
                    + 0.005353579108) * y - 0.002141268741) * y
                    + 0.000535310849) * y + 0.999936657524;
        }
    }
    return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}


/**************  ChiSquare: probability of chi square value *************/

/*ALGORITHM Compute probability of chi square value.
Adapted from: 	Hill, I. D. and Pike, M. C.  Algorithm 299.Collected Algorithms for the CACM 1967 p. 243
Updated for rounding errors based on remark inACM TOMS June 1985, page 185. Found in Perlman.lib*/

double computePValueChiSquare(double x, int df) /* x: obtained chi-square value,  df: degrees of freedom */ {
    double a, y, s;
    double e, c, z;
    int even; /* true if df is an even number */

    if (x <= 0.0 || df < 1)
        return (1.0);

    y = 1;

    a = 0.5 * x;
    even = (2 * (df / 2)) == df;
    if (df > 1)
        y = ex(-a);
    s = (even ? y : (2.0 * Normalz(-sqrt(x))));
    if (df > 2) {
        x = 0.5 * (df - 1.0);
        z = (even ? 1.0 : 0.5);
        if (a > BIGX) {
            e = (even ? 0.0 : LOG_SQRT_PI);
            c = log(a);
            while (z <= x) {
                e = log(z) + e;
                s += ex(c * z - a - e);
                z += 1.0;
            }
            return (s);
        } else {
            e = (even ? 1.0 : (I_SQRT_PI / sqrt(a)));
            c = 0.0;
            while (z <= x) {
                e = e * (a / z);
                c = c + e;
                z += 1.0;
            }
            return (c * y + s);
        }
    } else
        return (s);
}


int calculateSequenceHash(string &seq) {
	const static int modular = 1000000007;
	int hashValue = 0;
	for(char &c: seq) {
		hashValue = ((long long)hashValue * 107 + (int)c) % modular;
	}
	return hashValue;
}