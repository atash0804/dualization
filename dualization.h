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
size_t cov_count = 0;

enum {
    CHUNK_SIZE = sizeof(ull),
};
bool MODE = true; //set MODE=true if the existence of second coverage is in question

class customset {
    ull *data;
public:
    size_t sz;
    size_t chunks;

    explicit customset(size_t size) : sz(size) {
        chunks = (size + CHUNK_SIZE - 1) / CHUNK_SIZE;
        data = (ull *) calloc(chunks, sizeof(ull));
    }

    customset(set<size_t> source, size_t size) : sz(size) {
        chunks = (size + CHUNK_SIZE - 1) / CHUNK_SIZE;
        data = (ull *) calloc(chunks, sizeof(ull));
        for (auto entry: source) {
            this->sett(entry);
        }
    }

    customset(const customset &source) : sz(source.sz), chunks(source.chunks) {
        data = source.data;
    }

    customset &operator=(const customset &source) {
        data = source.data;
        sz = source.sz;
        return *this;
    }

//    ~customset() {
//        if (!used)
//            free(data);
//    }
    bool operator<(const customset &b) const {
        for (int i = 0; i < chunks; i++) {
            if (data[i] != b.data[i]) return (data[i] < b.data[i]);
        }
        return false;
    }

    void sett(ull k) {
        data[k / CHUNK_SIZE] |= ((ull) 1) << (k % CHUNK_SIZE);
    }

    void clear(ull k) {
        data[k / CHUNK_SIZE] &= ~((ull) 1 << (k % CHUNK_SIZE));
    }

    bool in(ull k) const {
        return bool((ull) 1 & (data[k / CHUNK_SIZE] >> (k % CHUNK_SIZE)));
    }

    friend bool check_intersection(customset s1, customset s2);
};

bool check_intersection(customset s1, customset s2) {
    if (s1.sz != s2.sz) {
        cerr << "Sets should be of equal size" << endl;
    }
    size_t chunks = (s1.sz + CHUNK_SIZE - 1) / CHUNK_SIZE;
    for (int i = 0; i < chunks; i++) {
        if (s1.data[i] & s2.data[i]) return false;
    }
    return true;
}

set<customset> default_coverage;
set<pair<set<size_t>, set<size_t>>> default_found_coverages;

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
    BitMatrix() : height(0), width(0), chunks(0) {}

    BitMatrix(istream &in, size_t n, size_t m) : height(n), width(m) {
        ull bit;
        chunks = width / CHUNK_SIZE + 1 - (width % CHUNK_SIZE == 0);

        for (int i = 0; i < height; i++) {
            vector<ull> row;
            row.clear();
            for (int j = 0; j < chunks; j++) {
                ull chunk = 0;
                for (int k = 1;
                     (k <= CHUNK_SIZE) && (j * CHUNK_SIZE + k <= width); k++) {
                    if (!(in >> bit)) {
                        cerr << "Incorrect size of input matrix" << endl;
                        throw length_error("");
                    }
                    if (bit > 1) {
                        cerr << "Incorrect value encountered in input stream: "
                             << bit << endl;
                        throw out_of_range("");
                    }
                    chunk |= (bit << (CHUNK_SIZE - k));
                }
                row.push_back(chunk);
            }
            matrix.push_back(row);
        }
    }

    BitMatrix(const string &filename, size_t n, size_t m) : height(n),
                                                            width(m) {
        ifstream in;
        in.open(filename);
        if (!in.is_open()) {
            cerr << "Error: " << strerror(errno) << endl;
            cerr << "Failed to open input file" << endl;
            throw bad_exception();
        }
        ull bit;
        chunks = width / CHUNK_SIZE + 1 - (width % CHUNK_SIZE == 0);

        for (int i = 0; i < height; i++) {
            vector<ull> row;
            row.clear();
            for (int j = 0; j < chunks; j++) {
                ull chunk = 0;
                for (size_t k = 1;
                     (k <= CHUNK_SIZE) && (j * CHUNK_SIZE + k <= width); k++) {
                    if (!(in >> bit)) {
                        cerr << "Incorrect size of input matrix" << endl;
                        throw length_error("");
                    }
                    if (bit > 1) {
                        cerr << "Incorrect value encountered in input stream: "
                             << bit << endl;
                        throw out_of_range("");
                    }
                    chunk |= (bit << (CHUNK_SIZE - k));
                }
                row.push_back(chunk);
            }
            matrix.push_back(row);
        }
        in.close();
    }

    ~BitMatrix() {
        for (auto a: matrix) a.clear();
        matrix.clear();
    }

    size_t getHeight() const {
        return height;
    }

    size_t getWidth() const {
        return width;
    }

    size_t getChunks() const {
        return chunks;
    }

    const vector<vector<ull>> &getMatrix() const {
        return matrix;
    }

    int at(size_t i, size_t j) const {
        size_t shift = CHUNK_SIZE - 1 - j % CHUNK_SIZE;
        return int((matrix[i][j / CHUNK_SIZE] & (1ULL << shift)) >> shift);
    }

    friend ostream &operator<<(ostream &os, const BitMatrix &bm);
};

ostream &operator<<(ostream &os, const BitMatrix &bm) {
    cout << "Bit matrix with shape (" << bm.height << ", " << bm.width << ")"
         << endl << endl;
    for (size_t i = 0; i < bm.height; i++) {
        for (size_t j = 0; j < bm.width; j++) {
            cout << bm.at(i, j) << ' ';
        }
        cout << endl;
    }
    cout << endl;
    return os;
}

class PartialBitMatrix : public BitMatrix {
    set<size_t> available_rows;
    set<size_t> available_cols;
    set<size_t> selected_cols;
    size_t cur_height;
    size_t cur_width;

    void delete_zero_columns() {
        vector<ull> rows_disjunction;

        for (size_t j = 0; j < this->getChunks(); j++)
            rows_disjunction.push_back(0);

        for (size_t i: available_rows) {
            for (size_t j = 0; j < this->getChunks(); j++) {
                rows_disjunction[j] |= this->getMatrix()[i][j];
            }
        }
        size_t shift;
        for (size_t j: available_cols) {
            shift = CHUNK_SIZE - 1 - j % CHUNK_SIZE;
            if ((rows_disjunction[j / CHUNK_SIZE] & (1ULL << shift)) >> shift ==
                0) {
                //cout << "Deleted zero column: " << j << endl;
                delete_column(j, false);
            }
        }
    }

    void delete_wider_rows() {
        bool fl; //fl is true if line1 >= line2
        size_t shift;
        set<size_t> wider_rows;
        for (size_t i: available_rows) {
            for (size_t j: available_rows) {
                if (i == j) continue;
                fl = true;
                for (size_t k: available_cols) {
                    shift = CHUNK_SIZE - 1 - k % CHUNK_SIZE;
                    if ((this->getMatrix()[i][k / CHUNK_SIZE] & (1ULL << shift))
                        < (this->getMatrix()[j][k / CHUNK_SIZE] &
                           (1ULL << shift))) {
                        fl = false;
                        break;
                    }
                }
                if (fl) {
                    //cout << "Deleted wider row: " << i << endl;
                    delete_row(i, false);
                }
            }
        }
    }

public:
    PartialBitMatrix(istream &in, size_t n, size_t m) : BitMatrix(in, n, m) {
        for (size_t i = 0; i < n; i++) available_rows.insert(i);
        for (size_t i = 0; i < m; i++) available_cols.insert(i);
        cur_height = this->getHeight();
        cur_width = this->getWidth();

        update_matrix();
    }

    PartialBitMatrix(const string &filename, size_t n, size_t m) : BitMatrix(
            filename, n, m) {
        for (size_t i = 0; i < n; i++) available_rows.insert(i);
        for (size_t i = 0; i < m; i++) available_cols.insert(i);
        cur_height = this->getHeight();
        cur_width = this->getWidth();

        update_matrix();
    }

    PartialBitMatrix() : cur_width(0), cur_height(0) {}

    void delete_column(size_t col, bool outside = true) {
        if (col >= this->getWidth()) {
            cerr << "Invalid column number" << endl;
            throw out_of_range("");
        }
        if (cur_width <= 0)
            cur_height = 0;
        else if (available_cols.find(col) != available_cols.end()) cur_width--;

        available_cols.erase(col);
        if (outside) selected_cols.insert(col);
    }

/**
 * @param row
 * @return
 */
    void delete_row(size_t row, bool outside = true) {
        if (row >= this->getHeight()) {
            cerr << "Invalid row number" << endl;
            throw out_of_range("");
        }
        if (cur_height <= 0)
            cur_height = 0;
        else if (available_rows.find(row) != available_rows.end()) cur_height--;

        available_rows.erase(row);
    }

    void update_matrix() {
        delete_zero_columns();
        delete_wider_rows();
    }

    size_t getCur_height() const {
        return cur_height;
    }

    size_t getCur_width() const {
        return cur_width;
    }

    const set<size_t> &getAvailable_rows() const {
        return available_rows;
    }

    const set<size_t> &getAvailable_cols() const {
        return available_cols;
    }

    const set<size_t> &getSelected_cols() const {
        return selected_cols;
    }

    pair<size_t, size_t> getLightestRow() const {
        ull tmp;
        size_t count, min = getHeight() + 1, min_n = 0;
        for (auto row: available_rows) {
            count = 0;
            for (int i = 0; i < getChunks(); i++) {
                tmp = getMatrix()[row][i];
                while (tmp != 0) {
                    tmp = tmp & (tmp - 1);
                    count++;
                }
            }
            if (count < min) {
                min = count;
                min_n = row;
            }
        }
        return {min_n, min};
    }

    friend ostream &operator<<(ostream &os, const PartialBitMatrix &bm);
};

ostream &operator<<(ostream &os, const PartialBitMatrix &bm) {
    cout << "Partial bit matrix with initial shape (" << bm.getHeight() << ", "
         << bm.getWidth() << ")" << endl;
    cout << "Current shape (" << bm.cur_height << ", " << bm.cur_width << ")"
         << endl << endl;
    for (auto i: bm.available_rows) {
        for (auto j: bm.available_cols) {
            cout << bm.at(i, j) << ' ';
        }
        cout << endl;
    }
    cout << "Available cols: [";
    for (auto entry: bm.getAvailable_cols()) cout << entry << ",";
    cout << ']' << endl;
    cout << "Available rows: [";
    for (auto entry: bm.getAvailable_rows()) cout << entry << ",";
    cout << ']' << endl;
    cout << "Selected cols: [";
    for (auto entry: bm.getSelected_cols()) cout << entry << ",";
    cout << ']' << endl;
    cout << endl;
    return os;
}

bool check_support_rows(PartialBitMatrix &L,
                        map<size_t, set<size_t>> &supporting_rows, size_t col) {
    set<size_t> one_rows;
    for (size_t i = 0; i < L.getHeight(); i++) {
        if (L.at(i, col)) one_rows.insert(i);
    }
    for (auto &entry: supporting_rows) {
        if (includes(one_rows.begin(), one_rows.end(), entry.second.begin(),
                     entry.second.end()))
            return false;
    }
    return true;
}

map<size_t, set<size_t>> update_support_rows(PartialBitMatrix &L,
                                             map<size_t, set<size_t>> &supporting_rows,
                                             size_t col) {
    set<size_t> one_rows;
    for (size_t i = 0; i < L.getHeight(); i++) {
        if (L.at(i, col)) one_rows.insert(i);
    }
    map<size_t, set<size_t>> copy = supporting_rows;
    for (auto &entry: copy) {
        for (auto row: entry.second) {
            if (one_rows.find(row) != one_rows.end()) {
                entry.second.erase(row);
                one_rows.erase(row);
            }
        }
    }
    copy[col] = one_rows;
    return copy;
}

void D1_dualization(PartialBitMatrix &L1, PartialBitMatrix &L2, \
    map<size_t, set<size_t>> &supporting_rows1,
                    map<size_t, set<size_t>> &supporting_rows2,
                    bool weights = false, \
    bool save = false,
                    set<pair<set<size_t>, set<size_t>>> &found_coverages = default_found_coverages) {
    //cout << "FIRST:" << L1 << "SECOND:" << L2 << endl << endl;
    PartialBitMatrix L_new;
    map<size_t, set<size_t>> supporting_rows_new;
    bool L1_empty, L2_empty, first = false;
    size_t row_number;

    L1_empty = L1.getCur_height() == 0;
    L2_empty = L2.getCur_height() == 0;
    if (L1_empty && L2_empty) {
        if (save) {
            found_coverages.insert(
                    {L1.getSelected_cols(), L2.getSelected_cols()});
        } else {
            cov_count++;
            printf("{");
            for (auto entry: L1.getSelected_cols()) {
                printf("%ld ", entry);
            }
            printf("}  {");
            for (auto entry: L2.getSelected_cols()) {
                printf("%ld ", entry);
            }
            printf("}\n");
        }
        return;
    }
    if (weights) {
        pair<size_t, size_t> res1, res2;
        if (!L1_empty) res1 = L1.getLightestRow();
        if (!L2_empty) res2 = L2.getLightestRow();
        if (res1.second >= res2.second) {
            row_number = res1.first;
            first = true;
        } else {
            row_number = res2.first;
            first = false;
        }
    } else {
        if (!L1_empty) {
            row_number = *L1.getAvailable_rows().begin();
            first = true;
            if (!L2_empty && (row_number > *L2.getAvailable_rows().begin())) {
                row_number = *L2.getAvailable_rows().begin();
                first = false;
            }
        } else {
            row_number = *L2.getAvailable_rows().begin();
        }
    }
    if (first) {
        for (size_t col: L1.getAvailable_cols()) {
            if (L1.at(row_number, col) && (L2.getSelected_cols().find(col) ==
                                           L2.getSelected_cols().end())) {
                if (check_support_rows(L1, supporting_rows1, col)) {
                    supporting_rows_new = update_support_rows(L1,
                                                              supporting_rows1,
                                                              col);
                    L_new = L1;
                    for (size_t i = 0; i < L1.getHeight(); i++) {
                        if (L1.at(i, col)) L_new.delete_row(i);
                    }
                    L_new.delete_column(col);
                    L_new.update_matrix();
                    D1_dualization(L_new, L2, supporting_rows_new,
                                   supporting_rows2, weights, save,
                                   found_coverages);
                }
            }
        }

    } else {
        for (size_t col: L2.getAvailable_cols()) {
            if (L2.at(row_number, col) && (L1.getSelected_cols().find(col) ==
                                           L1.getSelected_cols().end())) {
                if (check_support_rows(L2, supporting_rows2, col)) {
                    supporting_rows_new = update_support_rows(L2,
                                                              supporting_rows2,
                                                              col);
                    L_new = L2;
                    for (size_t i = 0; i < L2.getHeight(); i++) {
                        if (L2.at(i, col)) L_new.delete_row(i);
                    }
                    L_new.delete_column(col);
                    L_new.update_matrix();
                    D1_dualization(L1, L_new, supporting_rows1,
                                   supporting_rows_new, weights, save,
                                   found_coverages);
                }
            }
        }
    }
}

void D2_dualization(PartialBitMatrix &L1, PartialBitMatrix &L2, \
    map<size_t, set<size_t>> &supporting_rows1,
                    map<size_t, set<size_t>> &supporting_rows2,
                    bool weights = false, \
    bool save = false,
                    set<pair<set<size_t>, set<size_t>>> &found_coverages = default_found_coverages) {

    // cout << "FIRST:" << L1 << "SECOND:" << L2 << endl << endl;
    PartialBitMatrix L1_new, L2_new;
    map<size_t, set<size_t>> supporting_rows1_new, supporting_rows2_new;
    bool L1_empty, L2_empty;
    size_t row_number1, row_number2;
    set<size_t> one_cols1, one_cols2;

    L1_empty = L1.getCur_height() == 0;
    L2_empty = L2.getCur_height() == 0;

    if (L1_empty && L2_empty) {
        if (save) {
            found_coverages.insert(
                    {L1.getSelected_cols(), L2.getSelected_cols()});
        } else {
            cov_count++;
            printf("{");
            for (auto entry: L1.getSelected_cols()) {
                printf("%ld ", entry);
            }
            printf("}  {");
            for (auto entry: L2.getSelected_cols()) {
                printf("%ld ", entry);
            }
            printf("}\n");
        }
    }

    if (weights) {
        if (!L1_empty) {
            row_number1 = L1.getLightestRow().first;
            for (size_t col: L1.getAvailable_cols()) {
                if (L1.at(row_number1, col)) one_cols1.insert(col);
            }
        }
        if (!L2_empty) {
            row_number2 = L2.getLightestRow().first;
            for (size_t col: L2.getAvailable_cols()) {
                if (L2.at(row_number2, col)) one_cols2.insert(col);
            }
        }
    } else {
        if (!L1_empty) {
            row_number1 = *L1.getAvailable_rows().begin();
            for (size_t col: L1.getAvailable_cols()) {
                if (L1.at(row_number1, col)) one_cols1.insert(col);
            }
        }
        if (!L2_empty) {
            row_number2 = *L2.getAvailable_rows().begin();
            for (size_t col: L2.getAvailable_cols()) {
                if (L2.at(row_number2, col)) one_cols2.insert(col);
            }
        }
    }
    if (!L1_empty && !L2_empty) {
        for (size_t col1: one_cols1) {
            if ((L2.getSelected_cols().find(col1) !=
                 L2.getSelected_cols().end()) ||
                !check_support_rows(L1, supporting_rows1, col1))
                continue;
            for (size_t col2: one_cols2) {
                if ((col1 == col2) || (L1.getSelected_cols().find(col2) !=
                                       L1.getSelected_cols().end()))
                    continue;
                if (check_support_rows(L2, supporting_rows2, col2)) {
                    supporting_rows1_new = update_support_rows(L1,
                                                               supporting_rows1,
                                                               col1);
                    supporting_rows2_new = update_support_rows(L2,
                                                               supporting_rows2,
                                                               col2);
                    L1_new = L1;
                    L2_new = L2;
                    for (size_t i = 0; i < L1.getHeight(); i++) {
                        if (L1.at(i, col1)) L1_new.delete_row(i);
                    }
                    for (size_t i = 0; i < L2.getHeight(); i++) {
                        if (L2.at(i, col2)) L2_new.delete_row(i);
                    }
                    L1_new.delete_column(col1);
                    L1_new.update_matrix();
                    L2_new.delete_column(col2);
                    L2_new.update_matrix();
                    D2_dualization(L1_new, L2_new, supporting_rows1_new,
                                   supporting_rows2_new, weights, save,
                                   found_coverages);
                }
            }
        }
    } else {
        if (!L1_empty) {
            for (size_t col1: one_cols1) {
                if (L2.getSelected_cols().find(col1) !=
                    L2.getSelected_cols().end())
                    continue;
                if (check_support_rows(L1, supporting_rows1, col1)) {
                    supporting_rows1_new = update_support_rows(L1,
                                                               supporting_rows1,
                                                               col1);
                    L1_new = L1;
                    for (size_t i = 0; i < L1.getHeight(); i++) {
                        if (L1.at(i, col1)) L1_new.delete_row(i);
                    }
                    L1_new.delete_column(col1);
                    L1_new.update_matrix();
                    D2_dualization(L1_new, L2, supporting_rows1_new,
                                   supporting_rows2, weights, save,
                                   found_coverages);
                }
            }
        } else {
            for (size_t col2: one_cols2) {
                if (L1.getSelected_cols().find(col2) !=
                    L1.getSelected_cols().end())
                    continue;
                if (check_support_rows(L2, supporting_rows2, col2)) {
                    supporting_rows2_new = update_support_rows(L2,
                                                               supporting_rows2,
                                                               col2);
                    L2_new = L2;
                    for (size_t i = 0; i < L2.getHeight(); i++) {
                        if (L2.at(i, col2)) L2_new.delete_row(i);
                    }
                    L2_new.delete_column(col2);
                    L2_new.update_matrix();
                    D2_dualization(L1, L2_new, supporting_rows1,
                                   supporting_rows2_new, weights, save,
                                   found_coverages);
                }
            }
        }
    }
}

void print_results(set<pair<set<size_t>, set<size_t>>> &found_coverages) {
    for (auto &cov_pair: found_coverages) {
        printf("{");
        for (auto entry: cov_pair.first) {
            printf("%ld ", entry);
        }
        printf("}  {");
        for (auto entry: cov_pair.second) {
            printf("%ld ", entry);
        }
        printf("}\n");
    }
}

void
dualization(PartialBitMatrix &L1, map<size_t, set<size_t>> &supporting_rows1,
            bool weights = false, \
    bool save = false, set<customset> &coverages = default_coverage) {
    //cout << "FIRST:" << L1 << endl << endl;
    PartialBitMatrix L_new(L1);
    map<size_t, set<size_t>> supporting_rows_new;
    bool L1_empty;
    size_t row_number;

    L1_empty = L1.getCur_height() == 0;
    if (L1_empty) {
        if (save) {
            coverages.insert(customset(L1.getSelected_cols(), L1.getWidth()));
        } else {
            printf("{");
            for (auto entry: L1.getSelected_cols()) {
                printf("%ld ", entry);
            }
            printf("}\n");
        }
        return;
    }
    if (weights) {
        row_number = L1.getLightestRow().first;
    } else {
        row_number = *L1.getAvailable_rows().begin();
    }
    for (size_t col: L1.getAvailable_cols()) {
        if (L1.at(row_number, col)) {
            if (check_support_rows(L1, supporting_rows1, col)) {
                supporting_rows_new = update_support_rows(L1, supporting_rows1,
                                                          col);
                L_new = L1;
                for (size_t i = 0; i < L1.getHeight(); i++) {
                    if (L1.at(i, col)) L_new.delete_row(i);
                }
                L_new.delete_column(col);
                L_new.update_matrix();
                dualization(L_new, supporting_rows_new, weights, save,
                            coverages);
            }
        }
    }
}

void combine(set<customset> &cov1, set<customset> &cov2, bool print = true) {
    for (auto &set1: cov1) {
        for (auto &set2: cov2) {
            if (check_intersection(set1, set2)) {
                cov_count++;
                if (print) {
                    printf("{");
                    for (size_t i = 0; i < set1.sz; i++) {
                        if (set1.in(i)) printf("%ld ", i);
                    }
                    printf("}  {");
                    for (size_t i = 0; i < set2.sz; i++) {
                        if (set2.in(i)) printf("%ld ", i);
                    }
                    printf("}\n");
                }
                if (MODE) break;
            }
        }
    }
}

void generate_matrix(size_t n, size_t m, const string &filename,
                     double density = 0.5, int seed = -1) {
    if (seed != -1) srand(seed);
    ofstream out;
    out.open(filename);
    //cout << "Opened file:" << out.is_open() << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            out << 1 - int(rand() > RAND_MAX * density) << ' ';
        }
        out << endl;
    }
    out.close();
}

#endif //DUALIZATION_DUALIZATION_H