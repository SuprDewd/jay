#include <bits/stdc++.h>
#include "common.hpp"
#include "coeff.hpp"
#include "polynomial.hpp"
#include "matrix.hpp"
using namespace std;

bool HOPS = false;

struct combinations {
    int diff, fdeg, deg;

    vector<vi> diffs, pows;
    vi tmp;

    combinations(int diff, int fdeg, int deg) : diff(diff), fdeg(fdeg), deg(deg) {
    }

    void gen_diffs(int at) {
        if (at == VARS) {
            diffs.push_back(tmp);
            return;
        }
        rep(i,0,diff+1) {
            tmp.push_back(i);
            gen_diffs(at+1);
            tmp.pop_back();
        }
    }

    void gen_pows(int at) {
        if (at == VARS + size(diffs)) {
            pows.push_back(tmp);
            return;
        }
        int mx = at < VARS ? deg : fdeg;
        rep(i,0,mx+1) {
            tmp.push_back(i);
            gen_pows(at+1);
            tmp.pop_back();
        }
    }

    void generate() {
        // x^0y^0F^0F'^0...F(k)^0 + ... + x^degy^degF^fdegF'^fdeg...F(diff)^fdeg
        gen_diffs(0);
        gen_pows(0);
    }

    polynomial get_sum(polynomial F) {
        polynomial res;
        int cat = 0;
        iter(it,pows) {

            polynomial here;
            vi key(VARS);
            rep(i,0,VARS) key[i] = (*it)[i];
            vi ckey(COEFF_VARS);
            ckey[cat++] = 1;
            coeff c;
            c[ckey] = 1;
            here[key] = c;

            rep(i,0,size(diffs)) {
                polynomial Fp = F;
                rep(j,0,VARS) {
                    rep(k,0,diffs[i][j]) {
                        Fp = Fp.derivative(j);
                    }
                }
                here = here * Fp.pow((*it)[VARS + i]);
            }

            res = res + here;
        }
        return res;
    }

    void print(ostream &out, const vector<mpq_class> &vals) {
        bool first2 = true;
        int cat = 0;
        iter(it,pows) {
            mpq_class var = vals[cat];

            if (var == 0) {
                cat++;
                continue;
            }

            if (!first2) {
                if (var < 0) {
                    out << "-";
                    var = -var;
                } else {
                    out << "+";
                }
            }
            first2 = false;

            bool first = true;
            if (var != 1) {
                out << var;
                first = false;
            }

            rep(i,0,VARS) {
                if ((*it)[i] == 0) continue;
                if (!first) out << "*";
                first = false;
                out << get_var(i);
                if ((*it)[i] != 1) {
                    out << "^" << (*it)[i];
                }
            }
            rep(i,0,size(diffs)) {
                if ((*it)[VARS + i] == 0) continue;
                if (!first) out << "*";
                first = false;

                bool any = false;
                rep(j,0,VARS) {
                    if (diffs[i][j] > 0) {
                        any = true;
                        break;
                    }
                }
                if (any) {
                    if (HOPS) out << "D(";
                    else out << "D[";
                }

                out << "F";

                if (!HOPS) {
                    rep(j,0,VARS) {
                        if (diffs[i][j] == 0) continue;
                        out << ",{" << get_var(j) << "," << diffs[i][j] << "}";
                    }
                }

                if (any) {
                    if (HOPS) out << ")";
                    else out << "]";
                }

                if ((*it)[VARS + i] != 1) {
                    out << "^" << (*it)[VARS + i];
                }
            }

            if (first) out << "1";

            cat++;
        }
        if (first2) out << "0";
        out << "=0";
    }
};

void usage(char* ex) {
    fprintf(stderr, "Usage: %s [OPTION]...\n", ex);
    fprintf(stderr,
            "Given initial terms of an unknown generating function F,\n"
            "conjecture a functional equation that F satisfies.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -D,  --diff        the highest derivative of F to appear (default: 0)\n");
    fprintf(stderr, "  -fd, --fdeg        the highest power of F to appear (default: 1)\n");
    fprintf(stderr, "  -d,  --deg         the highest powers of variables to appear\n");
    fprintf(stderr, "  -t,  --type        the type of the generating function: ogf, egf (default: ogf)\n");
    fprintf(stderr, "       --hops        generate output that HOPS understands\n");
    fprintf(stderr, "  -h,  --help        display this help and exit\n");
    fprintf(stderr, "  -q,  --quiet       don't log any unnecessary output\n");
    fprintf(stderr, "\n");
    fprintf(stderr,
            "The program expects the initial terms to be specified on standard\n"
            "input in the following format. On the first line an integer 'n',\n"
            "denoting the number of variables in the generating function. On\n"
            "each of the following lines, until EOF, one initial term must be\n"
            "specified in the format 'i_1 i_2 ... i_n c_{i_1, i_2, ..., i_n}',\n"
            "meaning that 'c_{i_1, i_2, ..., i_n}' is the coefficient of\n"
            "'x_1^{i_1} x_2^{i_2} ... x_n^{i_n}' in F. The list of initial\n"
            "terms must be exhaustive, meaning that if a coefficient is listed,\n"
            "then all coefficients with smaller indices must be listed as\n"
            "well.\n"
            );
    exit(1);
}

int main(int argc, char* argv[]) {
    string type = "ogf";
    int diff = 0,
        fdeg = 1,
        deg = -1;
    rep(i,1,argc) {
        if (strcmp(argv[i], "--type") == 0 || strcmp(argv[i], "-t") == 0) {
            if (++i >= argc) error("missing argument to %s", argv[i-1]);
            type = argv[i];
            if (type != "ogf" && type != "egf") error("unrecognized generating function type %s", type);
        } else if (strcmp(argv[i], "--diff") == 0 || strcmp(argv[i], "-D") == 0) {
            if (++i >= argc) error("missing argument to %s", argv[i-1]);
            diff = atoi(argv[i]);
            if (diff < 0) error("argument to %s must be nonnegative", argv[i-1]);
        } else if (strcmp(argv[i], "--fdeg") == 0 || strcmp(argv[i], "-fd") == 0) {
            if (++i >= argc) error("missing argument to %s", argv[i-1]);
            fdeg = atoi(argv[i]);
            if (fdeg < 0) error("argument to %s must be nonnegative", argv[i-1]);
        } else if (strcmp(argv[i], "--deg") == 0 || strcmp(argv[i], "-d") == 0) {
            if (++i >= argc) error("missing argument to %s", argv[i-1]);
            deg = atoi(argv[i]);
            if (deg < 0) error("argument to %s must be nonnegative", argv[i-1]);
        } else if (strcmp(argv[i], "--hops") == 0) {
            HOPS = true;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
        } else if (strcmp(argv[i], "--quiet") == 0 || strcmp(argv[i], "-q") == 0) {
            QUIET = true;
        } else {
            error("unexpected argument: %s", argv[i]);
        }
    }

    if (deg == -1) error("missing required argument: --deg\n");

    cin >> VARS;

    if (HOPS && VARS != 1) error("HOPS output only supported for one variable");

    combinations comb(diff, fdeg, deg);
    comb.generate();
    COEFF_VARS = size(comb.pows);

    log("using %d variables%s", COEFF_VARS, COEFF_VARS > 1000 ? ", this might take awhile" : "");

    polynomial F = read_polynomial();

    if (type == "egf") {
        iter(it,F) {
            mpq_class fac = 1;
            iter(jt,it->first) {
                fac *= factorial(*jt);
            }
            coeff cfac;
            cfac[vi(COEFF_VARS, 0)] = mpq_class(1) / fac;
            it->second = it->second * cfac;
        }
    }

    polynomial sum = comb.get_sum(F);
    int rows = 0;
    iter(it,sum) rows++;
    rows = max(rows, COEFF_VARS);
    matrix<mpq_class> mat(rows, COEFF_VARS+1);
    int rat = 0;
    iter(it,sum) {
        iter(jt,it->second) {
            int who = -1;
            rep(i,0,COEFF_VARS) {
                if (jt->first[i] == 0) continue;
                if (jt->first[i] != 1 || who != -1) error("system of equations is not linear (bug?)");
                who = i;
            }
            if (who == -1) {
                mat(rat,COEFF_VARS) -= jt->second;
            } else {
                mat(rat,who) += jt->second;
            }
        }
        rat++;
    }
    // TODO: Maybe get rid of det and rank?
    mpq_class det;
    int rank;
    mat = mat.rref(det, rank);

    rep(i,0,rows) {
        bool zero = true;
        rep(j,0,COEFF_VARS) {
            if (mat(i,j) != 0) {
                zero = false;
                break;
            }
        }
        if (zero && mat(i,COEFF_VARS) != 0) {
            error("no solutions found");
        }
    }

    // TODO: Try some heuristics to get simple solutions (e.g. let coefficient for F be 1, and rest be 0)
    for (int i = 0, j = 0, at = rows-1; i < rows; i++, j++) {
        while (j < COEFF_VARS && mat(i,j) == 0) {
            mat(at, j++) = 1;
            mat(at--, COEFF_VARS) = 1;
        }
    }

    mat = mat.rref(det, rank);
    vector<mpq_class> vals(COEFF_VARS);
    rep(i,0,COEFF_VARS) vals[i] = mat(i,COEFF_VARS);

    comb.print(cout, vals);
    cout << endl;

    return 0;
}

