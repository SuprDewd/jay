#ifndef JAY_COMMON_HPP
#define JAY_COMMON_HPP

#include "gmp.h"
#include "gmpxx.h"
template <class T> int size(const T &x) { return x.size(); }
#define rep(i,a,b) for (__typeof(a) i=(a); i<(b); ++i)
#define iter(it,c) for (__typeof((c).begin()) it = (c).begin(); it != (c).end(); ++it)
typedef std::pair<int, int> ii;
typedef std::vector<int> vi;
typedef std::vector<ii> vii;
typedef long long ll;
const int INF = 2147483647;

void error(const char* format, ...) {
    va_list args;
    fprintf(stderr, "error: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n");
    exit(1);
}

bool QUIET = false;
void log(const char* format, ...) {
    if (QUIET) return;
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n");
}

mpq_class factorial(int n) {
    mpq_class res = 1;
    rep(i,2,n+1) {
        res *= i;
    }
    return res;
}

#endif
