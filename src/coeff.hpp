#ifndef JAY_COEFF_HPP
#define JAY_COEFF_HPP

#include "common.hpp"

int COEFF_VARS;

class coeff {
public:
    coeff() {
    }
    mpq_class& operator[](const vi &key) {
        assert(size(key) == COEFF_VARS);
        return val[key];
    }
    coeff operator+(const coeff &other) const {
        coeff res;
        iter(it,val) res[it->first] += it->second;
        iter(it,other.val) res[it->first] += it->second;
        res.clean();
        return res;
    }
    coeff operator-(const coeff &other) const {
        coeff res;
        iter(it,val) res[it->first] += it->second;
        iter(it,other.val) res[it->first] -= it->second;
        res.clean();
        return res;
    }
    coeff operator*(const coeff &other) const {
        coeff res;
        iter(it,val) iter(jt,other.val) {
            vi sm(COEFF_VARS);
            rep(i,0,COEFF_VARS) sm[i] = it->first[i] + jt->first[i];
            res[sm] += it->second * jt->second;
        }
        res.clean();
        return res;
    }
    coeff operator*(const mpq_class &other) const {
        coeff res;
        iter(it,val) res[it->first] = it->second * other;
        res.clean();
        return res;
    }
    bool is_zero() const {
        iter(it,val) {
            if (it->second != 0) {
                return false;
            }
        }
        return true;
    }
    bool is_one() const {
        iter(it,val) {
            if (it->second == 0) continue;
            if (it->second != 1) return false;
            iter(jt,it->first) {
                if (*jt != 0) {
                    return false;
                }
            }
        }
        return true;
    }
    bool is_sum() const {
        bool found = false;
        iter(it,val) {
            if (it->second == 0) continue;
            if (found) return true;
            found = true;
        }
        return false;
    }
    std::map<vi,mpq_class>::iterator begin() { return val.begin(); }
    std::map<vi,mpq_class>::iterator end() { return val.end(); }
private:
    std::map<vi,mpq_class> val;
    void clean() {
        for (auto it = val.begin(); it != val.end(); ) {
            if (it->second == 0) {
                it = val.erase(it);
            } else {
                it++;
            }
        }
    }

    friend std::ostream& operator <<(std::ostream &out, const coeff &p);
};

std::string get_coeff_var(int i) {
    if (COEFF_VARS < 23) {
        return std::string("") + static_cast<char>('a' + i);
    } else {
        std::stringstream ss;
        ss << "a";
        ss << i;
        return ss.str();
    }
}

std::ostream& operator <<(std::ostream &out, const coeff &p) {
    // TODO: Handle negative things
    bool first2 = true;
    iter(it,p.val) {
        if (it->second == 0) continue;
        if (!first2) out << "+";
        first2 = false;

        bool first = true;

        if (it->second != 1) {
            first = false;
            out << it->second;
        }

        rep(i,0,COEFF_VARS) {
            if (it->first[i] == 0) continue;
            if (!first) out << "*";
            first = false;
            out << get_coeff_var(i);
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
