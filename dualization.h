#ifndef DUALIZATION_DUALIZATION_H
#define DUALIZATION_DUALIZATION_H

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

using namespace std;

typedef uint64_t ull;
extern int cov_count;

class customset {
    ull *data;
public:
    size_t sz;
    size_t chunks;

    explicit customset(size_t size);
    customset(set<size_t> source, size_t size);
    customset(const customset& source);
    customset& operator=(const customset& source);

    bool operator<(const customset &b);
    void sett(ull k);
    void clear(ull k);
    bool in(ull k);
    friend bool check_intersection(customset s1, customset s2);
};

bool check_intersection(customset s1, customset s2);

/**
 * A simple class for binary matrix stored as bit sets
 */
class BitMatrix {
    /*! The matrix where all data is stored*/
    vector<vector<ull>> matrix;
    size_t height;
    size_t width;
    size_t chunks;

public:
    BitMatrix();
    BitMatrix(istream& in, size_t n, size_t m);
    BitMatrix(const string& filename, size_t n, size_t m);
    ~BitMatrix();
    size_t getHeight();
    size_t getWidth();
    size_t getChunks();
    const vector<vector<ull>> &getMatrix();
    int at (size_t i, size_t j);
    friend ostream& operator<<(ostream& os, const BitMatrix& bm);
};

ostream& operator<<(ostream& os, const BitMatrix& bm);

class PartialBitMatrix: public BitMatrix {
    set<size_t> available_rows;
    set<size_t> available_cols;
    set<size_t> selected_cols;
    size_t cur_height;
    size_t cur_width;

    void delete_zero_columns();
    void delete_wider_rows();

public:
    PartialBitMatrix(istream& in, size_t n, size_t m);
    PartialBitMatrix(const string& filename, size_t n, size_t m);
    PartialBitMatrix();
    void delete_column(size_t col, bool outside = true);
    void delete_row(size_t row, bool outside = true);
    void update_matrix();
    size_t getCur_height() const;
    size_t getCur_width() const;
    const set<size_t> &getAvailable_rows() const;
    const set<size_t> &getAvailable_cols() const;
    const set<size_t> &getSelected_cols() const;
    pair<size_t, size_t> getLightestRow() const;
    friend ostream& operator<<(ostream& os, const PartialBitMatrix& bm);
};

ostream& operator<<(ostream& os, const PartialBitMatrix& bm);

bool check_support_rows(PartialBitMatrix& L, map<size_t, set<size_t>>& supporting_rows, size_t col);

map<size_t, set<size_t>> update_support_rows(PartialBitMatrix& L, map<size_t, set<size_t>>& supporting_rows, size_t col);

void D1_dualization (PartialBitMatrix& L1, PartialBitMatrix& L2, \
    map<size_t, set<size_t>> &supporting_rows1, map<size_t, set<size_t>> &supporting_rows2, bool weights, \
    bool save, set<pair<set<size_t >, set<size_t>>>& found_coverages);

void D2_dualization (PartialBitMatrix& L1, PartialBitMatrix& L2, \
    map<size_t, set<size_t>> &supporting_rows1, map<size_t, set<size_t>> &supporting_rows2, bool weights, \
    bool save, set<pair<set<size_t >, set<size_t>>>& found_coverages);

void print_results(set<pair<set<size_t >, set<size_t>>>& found_coverages);

void dualization (PartialBitMatrix& L1, map<size_t, set<size_t>> &supporting_rows1, bool weights, \
    bool save, set<customset>& coverages);

void combine(set <customset>& cov1, set <customset>& cov2, bool print = true);

void generate_matrix(size_t n, size_t m, const string& filename, double density = 0.5, int seed = -1);

#endif //DUALIZATION_DUALIZATION_H
