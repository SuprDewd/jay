#ifndef JAY_POLYNOMIAL_HPP
#define JAY_POLYNOMIAL_HPP

#include "common.hpp"
#include "coeff.hpp"

int VARS;

class polynomial {
public:
    polynomial() : know(VARS, INF) {
    }
    coeff& operator[](const vi &key) {
        assert(size(key) == VARS);
        return val[key];
    }
    polynomial operator+(const polynomial &other) const {
        polynomial res;
        iter(it,val) res[it->first] = res[it->first] + it->second;
        iter(it,other.val) res[it->first] = res[it->first] + it->second;
        rep(i,0,VARS) res.know[i] = std::min(know[i], other.know[i]);
        res.clean();
        return res;
    }
    polynomial operator-(const polynomial &other) const {
        polynomial res;
        iter(it,val) res[it->first] = res[it->first] + it->second;
        iter(it,other.val) res[it->first] = res[it->first] - it->second;
        rep(i,0,VARS) res.know[i] = std::min(know[i], other.know[i]);
        res.clean();
        return res;
    }
    polynomial operator*(const polynomial &other) const {
        polynomial res;
        iter(it,val) iter(jt,other.val) {
            vi sm(VARS);
            rep(i,0,VARS) sm[i] = it->first[i] + jt->first[i];
            res[sm] = res[sm] + it->second * jt->second;
        }
        // TODO: This is too restricted when some of the lower terms are known to be zero (consider x * (1 + x + O(x^2)) which is currently x + O(x^2), but should be x + x^2 + O(x^3))
        rep(i,0,VARS) res.know[i] = std::min(know[i], other.know[i]);
        res.clean();
        return res;
    }
    polynomial derivative(int resp) {
        polynomial res;
        iter(it,val) {
            if (it->first[resp] == 0) continue;
            vi other = it->first;
            other[resp]--;
            res[other] = it->second * mpq_class(it->first[resp]);
        }
        res.know = know;
        if (res.know[resp] != INF) res.know[resp] = std::max(-1, res.know[resp]-1);
        res.clean();
        return res;
    }
    polynomial pow(int p) {
        polynomial res, base = *this;
        coeff one;
        one[vi(COEFF_VARS,0)] = 1;
        res[vi(VARS,0)] = one;
        while (p) {
            if (p & 1) {
                res = res * base;
            }
            base = base * base;
            p >>= 1;
        }
        return res;
    }
    std::map<vi,coeff>::iterator begin() { return val.begin(); }
    std::map<vi,coeff>::iterator end() { return val.end(); }
private:
    vi know;
    std::map<vi,coeff> val;
    void clean() {
        for (auto it = val.begin(); it != val.end(); ) {
            bool erase = it->second.is_zero();
            if (!erase) {
                rep(i,0,VARS) {
                    if (it->first[i] > know[i]) {
                        erase = true;
                        break;
                    }
                }
            }
            if (erase) {
                it = val.erase(it);
            } else {
                it++;
            }
        }
    }

    friend polynomial read_polynomial();
    friend std::ostream& operator <<(std::ostream &out, const polynomial &p);
};

struct validate_polynomial_complete {
    vi &deg, cur;
    std::set<vi> &seen;
    validate_polynomial_complete(vi &_deg, std::set<vi> &_seen) : deg(_deg), cur(VARS), seen(_seen) {
    }
    bool bt(int at) {
        if (at == VARS) {
            return seen.find(cur) != seen.end();
        }
        rep(i,0,deg[at]+1) {
            cur[at] = i;
            if (bt(at+1)) {
                return true;
            }
        }
        return false;
    }
    bool validate() {
        return bt(0);
    }
};

polynomial read_polynomial() {
    polynomial res;
    std::set<vi> seen;
    vi deg(VARS, -1);
    while (true) {
        vi key(VARS);
        bool ok = true;
        rep(i,0,VARS) {
            if (!(std::cin >> key[i])) {
                if (i != 0) error("badly formatted input");
                ok = false;
                break;
            }
        }
        if (!ok) break;
        if (seen.find(key) != seen.end()) error("duplicate coefficient");
        seen.insert(key);
        mpq_class val;
        if (!(std::cin >> val)) error("badly formatted input");
        coeff cval;
        cval[vi(COEFF_VARS, 0)] = val;
        res[key] = cval;
        rep(i,0,VARS) deg[i] = std::max(deg[i], key[i]);
    }
    if (!validate_polynomial_complete(deg, seen).validate()) {
        error("missing terms in input polynomial");
    }
    res.know = deg;
    res.clean();
    return res;
}

std::string get_var(int i) {
    if (VARS < 3) {
        return std::string("") + static_cast<char>('x' + i);
    } else {
        std::stringstream ss;
        ss << "x";
        ss << i;
        return ss.str();
    }
}

std::ostream& operator <<(std::ostream &out, const polynomial &p) {
    // TODO: Handle negative things
    bool first2 = true;
    iter(it,p.val) {
        if (it->second.is_zero()) continue;
        if (!first2) out << "+";
        first2 = false;

        bool first = true;

        if (!it->second.is_one()) {
            first = false;
            if (it->second.is_sum()) out << "(";
            out << it->second;
            if (it->second.is_sum()) out << ")";
        }

        rep(i,0,VARS) {
            if (it->first[i] == 0) continue;
            if (!first) out << "*";
            first = false;
            out << get_var(i);
            if (it->first[i] != 1) {
                out << "^" << it->first[i];
            }
        }

        if (first) {
            out << "1";
        }
    }
    if (first2) out << "0";
    return out;
}

#endif
