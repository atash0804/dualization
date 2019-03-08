#include <map>
#include <set>
#include <ctime>
#include <cerrno>
#include <string>
#include <random>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "dualization.h"

using namespace std;

int main() {
    ofstream out;
    clock_t stop, start;
    double elapsed;
    out.open("times");
    srand(clock()/ CLOCKS_PER_SEC);

    //L1 has shape (m, n), L2 has shape (l, n)
    vector<size_t> n_vec = {10, 11, 12, 13};
    vector<size_t> m_vec = {10, 11, 12, 13};
    vector<size_t> l_vec = {10, 11, 12, 13};

    set<customset> cov1, cov2;
    map<size_t, set<size_t>> supporting_rows1, supporting_rows2;
    set<pair<set<size_t>, set<size_t>>> found_coverages;
    map<size_t, set<size_t>> supporting_rows3, supporting_rows4;

    for (size_t i = 0; i < 4; i++) {
        size_t n = n_vec[i];
        size_t m = m_vec[i];
        size_t l = l_vec[i];

        generate_matrix(m, n, "matrix1.txt", 0.5);
        generate_matrix(l, n, "matrix2.txt", 0.5);

        out << "N = " << n << " M = " << m << " L = " << l << endl;
        cout << "N = " << n << " M = " << m << " L = " << l << endl;

        out << "DUALIZATION:" << endl;
        PartialBitMatrix matr1 = PartialBitMatrix("matrix1.txt", m, n);
        PartialBitMatrix matr2 = PartialBitMatrix("matrix2.txt", l, n);
        cov1.clear();
        cov2.clear();
        supporting_rows1.clear();
        supporting_rows2.clear();

        start = clock();
        dualization(matr1, supporting_rows1, true, true, cov1);
        dualization(matr2, supporting_rows2, true, true, cov2);
        stop = clock();
        elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        out << "Dualization time: " << elapsed << endl;
        cov_count = 0;
        combine(cov1, cov2, false);
        out << "Cov total(1): " << cov_count << endl;

        stop = clock();
        elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        out << "Overall time: " << elapsed << endl;
        out << endl;

        out << "DOUBLE DUALIZATION" << endl;
        PartialBitMatrix matr3 = PartialBitMatrix("matrix1.txt", m, n);
        PartialBitMatrix matr4 = PartialBitMatrix("matrix2.txt", l, n);

        found_coverages.clear();
        supporting_rows3.clear();
        supporting_rows4.clear();

        start = clock();

        cov_count = 0;
        D1_dualization(matr3, matr4, supporting_rows3, supporting_rows4, true, true, found_coverages);
        out << "Cov total(2): " << found_coverages.size() << endl;

        stop = clock();
        elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        out << "Double time:  " << elapsed << endl;

        out << "_______________________________________________" << endl << endl;
    }


    out.close();
    return 0;
}