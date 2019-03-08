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

enum {
    CHUNK_SIZE = sizeof(ull)
};

class customset {
    size_t sz;
    ull *data;
public:
    explicit customset(size_t size = 0): sz(size) {
        size_t chunks = (size + CHUNK_SIZE - 1) / CHUNK_SIZE;
        data = (ull*)calloc(chunks, sizeof(ull));
    }

    customset (const customset& right) {
        sz = right.sz;
        size_t chunks = (right.sz + CHUNK_SIZE - 1) / CHUNK_SIZE;
        data = (ull*)calloc(chunks, sizeof(ull));
        for (size_t i = 0; i < chunks; i++) {
            data[i] = right.data[i];
        }
    }

    ~customset() {
        free(data);
    }

    void set(ull k) {
        data[k / CHUNK_SIZE] |= ((ull)1) << (k % CHUNK_SIZE);
    }

    void clear(ull k) {
        data[k / CHUNK_SIZE] &= ~(1 << (k % CHUNK_SIZE));
    }

    bool in(ull k) const {
        return bool(1 & (data[k / CHUNK_SIZE] >> (k % CHUNK_SIZE)));
    }

    size_t getSize() const {
        return sz;
    }

    customset& operator=(const customset& right) {
        size_t chunks = (right.sz + CHUNK_SIZE - 1) / CHUNK_SIZE;
        free(data);
        data = (ull*)calloc(chunks, sizeof(ull));
        for (size_t i = 0; i < chunks; i++) {
            data[i] = right.data[i];
        }
        return *this;
    }
};

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
    BitMatrix(): height(10), width(10), chunks(10) {}

    BitMatrix(istream& in, size_t n, size_t m):height(n), width(m) {
        ull bit;
        chunks = width / CHUNK_SIZE + 1 - (width % CHUNK_SIZE == 0);

        for (int i = 0; i < height; i++) {
            vector <ull> row;
            row.clear();
            for (int j = 0; j < chunks; j++) {
                ull chunk = 0;
                for (int k = 1; (k <= CHUNK_SIZE) && (j * CHUNK_SIZE + k <= width); k++) {
                    if (!(in >> bit)) {
                        cerr << "Incorrect size of input matrix" << endl;
                        throw length_error("");
                    }
                    if (bit > 1) {
                        cerr << "Incorrect value encountered in input stream: " << bit << endl;
                        throw out_of_range("");
                    }
                    chunk |= (bit << (CHUNK_SIZE - k));
                }
                row.push_back(chunk);
            }
            matrix.push_back(row);
        }
    }

    BitMatrix(const string& filename, size_t n, size_t m):height(n), width(m) {
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
            vector <ull> row;
            row.clear();
            for (int j = 0; j < chunks; j++) {
                ull chunk = 0;
                for (int k = 1; (k <= CHUNK_SIZE) && (j * CHUNK_SIZE + k <= width); k++) {
                    if (!(in >> bit)) {
                        cerr << "Incorrect size of input matrix" << endl;
                        throw length_error("");
                    }
                    if (bit > 1) {
                        cerr << "Incorrect value encountered in input stream: " << bit << endl;
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

    int at (size_t i, size_t j) const {
        size_t shift = CHUNK_SIZE - 1 - j % CHUNK_SIZE;
        return int((matrix[i][j / CHUNK_SIZE] & (1ULL << shift)) >> shift);
    }

    friend ostream& operator<<(ostream& os, const BitMatrix& bm);
};

ostream& operator<<(ostream& os, const BitMatrix& bm) {
    cout << "Bit matrix with shape (" << bm.height << ", " << bm.width << ")" << endl << endl;
    for (size_t i = 0; i < bm.height; i++) {
        for (size_t j = 0; j < bm.width; j++) {
            cout << bm.at(i, j) << ' ';
        }
        cout << endl;
    }
    cout << endl;
    return os;
}

class PartialBitMatrix: public BitMatrix {
    customset available_rows;
    customset available_cols;
    customset selected_cols;
    size_t cur_height;
    size_t cur_width;

    void delete_zero_columns() {
        vector <ull> rows_disjunction;

        for (size_t j = 0; j < this->getChunks(); j++) rows_disjunction.push_back(0);

        for (size_t i = 0; i < available_rows.getSize(); i++) {
            if (!available_rows.in(i)) continue;
            for (size_t j = 0; j < this->getChunks(); j++) {
                rows_disjunction[j] |= this->getMatrix()[i][j];
            }
        }
        size_t shift;
        for (size_t j = 0; j < available_cols.getSize(); j++) {
            if (!available_cols.in(j)) continue;
            shift = CHUNK_SIZE - 1 - j % CHUNK_SIZE;
            if ((rows_disjunction[j / CHUNK_SIZE] & (1ULL << shift)) >> shift == 0) {
                //cout << "Deleted zero column: " << j << endl;
                delete_column(j, false);
            }
        }
    }

    void delete_wider_rows() {
        bool fl; //fl is true if line1 >= line2
        size_t shift;
        for (size_t i = 0; i < available_rows.getSize(); i++) {
            if (!available_rows.in(i)) continue;
            for (size_t j = 0; j < available_rows.getSize(); j++) {
                if (!available_rows.in(j) || (i == j)) continue;
                fl = true;
                for (size_t k = 0; k < available_cols.getSize(); k++) {
                    if (!available_cols.in(k)) continue;
                    shift = CHUNK_SIZE - 1 - k % CHUNK_SIZE;
                    if ((this->getMatrix()[i][k / CHUNK_SIZE] & (1ULL << shift))
                        < (this->getMatrix()[j][k / CHUNK_SIZE] & (1ULL << shift))) {
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
    PartialBitMatrix(istream& in, size_t n, size_t m): BitMatrix(in, n, m), available_rows(n), available_cols(m), selected_cols(m) {
        for (size_t i = 0; i < n; i++) available_rows.set(i);
        for (size_t i = 0; i < m; i++) available_cols.set(i);
        cur_height = this->getHeight();
        cur_width = this->getWidth();

        update_matrix();
    }

    PartialBitMatrix(const string& filename, size_t n, size_t m): BitMatrix(filename, n, m), available_rows(n), available_cols(m), selected_cols(m) {
        for (size_t i = 0; i < n; i++) available_rows.set(i);
        for (size_t i = 0; i < m; i++) available_cols.set(i);
        cur_height = this->getHeight();
        cur_width = this->getWidth();

        update_matrix();
    }

    PartialBitMatrix():cur_width(0), cur_height(0), available_rows(0), available_cols(0), selected_cols(0) {}

    void delete_column(size_t col, bool outside = true) {
        if (col >= this->getWidth()) {
            cerr << "Invalid column number" << endl;
            throw out_of_range("");
        }
        if (cur_width <= 0)
            cur_height = 0;
        else
            if (available_cols.in(col)) cur_width--;

        available_cols.clear(col);
        if (outside) selected_cols.set(col);
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
        if (cur_height <= 1)
            cur_width = 0;
        else
            if (available_rows.in(row)) cur_height--;

        available_rows.clear(row);
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

    const customset &getAvailable_rows() const {
        return available_rows;
    }

    const customset &getAvailable_cols() const {
        return available_cols;
    }

    const customset &getSelected_cols() const {
        return selected_cols;
    }

    friend ostream& operator<<(ostream& os, const PartialBitMatrix& bm);
};

ostream& operator<<(ostream& os, const PartialBitMatrix& bm) {
    cout << "Partial bit matrix with initial shape (" << bm.getHeight() << ", " << bm.getWidth() << ")" << endl;
    cout << "Current shape (" << bm.cur_height << ", " << bm.cur_width << ")" << endl << endl;

    for (size_t i = 0; i < bm.available_rows.getSize(); i++) {
        if (!bm.available_rows.in(i)) continue;
        for (size_t j = 0; j < bm.available_cols.getSize(); j++) {
            if (!bm.available_cols.in(j)) continue;
            cout << bm.at(i, j) << ' ';
        }
        cout << endl;
    }
    cout << "Available cols: [";
    for (size_t j = 0; j < bm.available_cols.getSize(); j++) {
        if (!bm.available_cols.in(j)) continue;
        cout << j << ",";
    }
    cout << ']' << endl;
    cout << "Available rows: [";
    for (size_t j = 0; j < bm.available_rows.getSize(); j++) {
        if (!bm.available_rows.in(j)) continue;
        cout << j << ",";
    }
    cout << ']' << endl;
    cout << "Selected cols: [";
    for (size_t j = 0; j < bm.selected_cols.getSize(); j++) {
        if (!bm.selected_cols.in(j)) continue;
        cout << j << ",";
    }
    cout << ']' << endl;
    cout << endl;
    return os;
}

bool check_support_rows(PartialBitMatrix& L, map<size_t, customset>& supporting_rows, size_t col) {
    customset one_rows(L.getHeight());
    bool fl;
    for (size_t i = 0; i < L.getHeight(); i++) {
        if (L.at(i, col)) one_rows.set(i);
    }
    for (auto& entry: supporting_rows) {
        fl = true;
        for (size_t j = 0; j < L.getHeight(); j++) {
            if (entry.second.in(j) && !one_rows.in(j)) {
                fl = false;
                break;
            }
        }
        if (fl) return false;
    }
    return true;
}

map<size_t, customset> update_support_rows(PartialBitMatrix& L, map<size_t, customset>& supporting_rows, size_t col) {
    customset one_rows(L.getHeight());
    for (size_t i = 0; i < L.getHeight(); i++) {
        if (L.at(i, col)) one_rows.set(i);
    }
    map<size_t, customset> copy = supporting_rows;
    for (auto& entry: copy) {
        for (size_t row = 0; row < entry.second.getSize(); row++) {
            if (!entry.second.in(row)) continue;
            if (one_rows.in(row)) {
                entry.second.clear(row);
                one_rows.clear(row);
            }
        }
    }
    copy[col] = one_rows;
    return copy;
}

void D1_dualization (PartialBitMatrix& L1, PartialBitMatrix& L2, \
    map<size_t, customset> &supporting_rows1, map<size_t, customset> &supporting_rows2) {

    cout << "FIRST:" << L1 << "SECOND:" << L2 << endl << endl;
    PartialBitMatrix L_new;
    map<size_t, customset> supporting_rows_new;
    bool L1_empty, L2_empty, first = false;
    size_t row_number1 = 0, row_number2 = 0;

    L1_empty = L1.getCur_height() == 0;
    L2_empty = L2.getCur_height() == 0;
    if (L1_empty && L2_empty) {
        cout << "{ ";
        for (size_t entry = 0; entry < L1.getSelected_cols().getSize(); entry++) {
            if (!L1.getSelected_cols().in(entry)) continue;
            cout << entry << ' ';
        }
        cout << '}' << "    { ";
        for (size_t entry = 0; entry < L2.getSelected_cols().getSize(); entry++) {
            if (!L2.getSelected_cols().in(entry)) continue;
            cout << entry << ' ';
        }
        cout << '}' << endl;
        return;
    }
    if (!L1_empty) {
        for (size_t i = 0; i < L1.getAvailable_rows().getSize(); i++) {
            if (L1.getAvailable_rows().in(i)) row_number1 = i;
        }
    }
    if (!L2_empty) {
        for (size_t i = 0; i < L2.getAvailable_rows().getSize(); i++) {
            if (L2.getAvailable_rows().in(i)) row_number2 = i;
        }
    }
    if (row_number1 <= row_number2) {
        first = true;
    }
    if (first) {
        for (size_t col = 0; col < L1.getAvailable_cols().getSize(); col++) {
            if (!L1.getAvailable_cols().in(col)) continue;
            if (L1.at(row_number1, col) && (!L2.getSelected_cols().in(col))) {
                if (check_support_rows(L1, supporting_rows1, col)) {
                    supporting_rows_new = update_support_rows(L1, supporting_rows1, col);
                    L_new = L1;
                    for (size_t i = 0; i < L1.getHeight(); i++) {
                        if (L1.at(i, col)) L_new.delete_row(i);
                    }
                    L_new.delete_column(col);
                    L_new.update_matrix();
                    D1_dualization(L_new, L2, supporting_rows_new, supporting_rows2);
                }
            }
        }

    } else {
        for (size_t col = 0; col < L2.getAvailable_cols().getSize(); col++) {
            if (!L2.getAvailable_cols().in(col)) continue;
            if (L2.at(row_number2, col) && (!L1.getSelected_cols().in(col))) {
                if (check_support_rows(L2, supporting_rows2, col)) {
                    supporting_rows_new = update_support_rows(L2, supporting_rows2, col);
                    L_new = L2;
                    for (size_t i = 0; i < L2.getHeight(); i++) {
                        if (L2.at(i, col)) L_new.delete_row(i);
                    }
                    L_new.delete_column(col);
                    L_new.update_matrix();
                    D1_dualization(L1, L_new, supporting_rows1, supporting_rows_new);
                }
            }
        }
    }
}

void D2_dualization (PartialBitMatrix& L1, PartialBitMatrix& L2, \
    map<size_t, customset> &supporting_rows1, map<size_t, customset> &supporting_rows2) {

    cout << "FIRST:" << L1 << "SECOND:" << L2 << endl << endl;
    PartialBitMatrix L1_new, L2_new;
    map<size_t, customset> supporting_rows1_new, supporting_rows2_new;
    bool L1_empty, L2_empty;
    size_t row_number1 = 0, row_number2 = 0;
    customset one_cols1(L1.getWidth()), one_cols2(L2.getWidth());

    L1_empty = L1.getCur_height() == 0;
    L2_empty = L2.getCur_height() == 0;

    if (L1_empty && L2_empty) {
        cout << "{ ";
        for (size_t entry = 0; entry < L1.getSelected_cols().getSize(); entry++) {
            if (!L1.getSelected_cols().in(entry)) continue;
            cout << entry << ' ';
        }
        cout << '}' << "    { ";
        for (size_t entry = 0; entry < L2.getSelected_cols().getSize(); entry++) {
            if (!L2.getSelected_cols().in(entry)) continue;
            cout << entry << ' ';
        }
        cout << '}' << endl;
        return;
    }

    if (!L1_empty) {
        for (size_t i = 0; i < L1.getAvailable_rows().getSize(); i++) {
            if (L1.getAvailable_rows().in(i)) row_number1 = i;
        }
        for (size_t col = 0; col < L1.getAvailable_cols().getSize(); col++) {
            if (!L1.getAvailable_cols().in(col)) continue;
            if (L1.at(row_number1, col)) one_cols1.set(col);
        }
    }
    if (!L2_empty) {
        for (size_t i = 0; i < L2.getAvailable_rows().getSize(); i++) {
            if (L2.getAvailable_rows().in(i)) row_number2 = i;
        }
        for (size_t col = 0; col < L2.getAvailable_cols().getSize(); col++) {
            if (!L2.getAvailable_cols().in(col)) continue;
            if (L2.at(row_number2, col)) one_cols2.set(col);
        }
    }
    if (!L1_empty && !L2_empty) {
        for (size_t col1 = 0; col1 < one_cols1.getSize(); col1++) {
            if (!one_cols1.in(col1)) continue;
            if ((L2.getSelected_cols().in(col1)) || !check_support_rows(L1, supporting_rows1, col1)) continue;
            for (size_t col2 = 0; col2 < one_cols2.getSize(); col2++) {
                if (!one_cols2.in(col2)) continue;
                if ((col1 == col2) || (L1.getSelected_cols().in(col2))) continue;
                if (check_support_rows(L2, supporting_rows2, col2)) {
                    supporting_rows1_new = update_support_rows(L1, supporting_rows1, col1);
                    supporting_rows2_new = update_support_rows(L2, supporting_rows2, col2);
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
                    D2_dualization(L1_new, L2_new, supporting_rows1_new, supporting_rows2_new);
                }
            }
        }
    } else {
        if (!L1_empty) {
            for (size_t col1 = 0; col1 < one_cols1.getSize(); col1++) {
                if (!one_cols1.in(col1)) continue;
                if (L2.getSelected_cols().in(col1)) continue;
                if (check_support_rows(L1, supporting_rows1, col1)) {
                    supporting_rows1_new = update_support_rows(L1, supporting_rows1, col1);
                    L1_new = L1;
                    for (size_t i = 0; i < L1.getHeight(); i++) {
                        if (L1.at(i, col1)) L1_new.delete_row(i);
                    }
                    L1_new.delete_column(col1);
                    L1_new.update_matrix();
                    D2_dualization(L1_new, L2, supporting_rows1_new, supporting_rows2);
                }
            }
        } else {for (size_t col2 = 0; col2 < one_cols2.getSize(); col2++) {
                if (!one_cols2.in(col2)) continue;
                if (L1.getSelected_cols().in(col2)) continue;
                if (check_support_rows(L2, supporting_rows2, col2)) {
                    supporting_rows2_new = update_support_rows(L2, supporting_rows2, col2);
                    L2_new = L2;
                    for (size_t i = 0; i < L2.getHeight(); i++) {
                        if (L2.at(i, col2)) L2_new.delete_row(i);
                    }
                    L2_new.delete_column(col2);
                    L2_new.update_matrix();
                    D2_dualization(L1, L2_new, supporting_rows1, supporting_rows2_new);
                }
            }
        }
    }
}

void generate_matrix(size_t n, size_t m, const string& filename, unsigned int seed = 6, double density = 0.5) {
    srand(seed);
    ofstream out;
    out.open(filename);
    cout << out.is_open() << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            out << 1 - int(rand() > RAND_MAX*density) << ' ';
        }
        out << endl;
    }
    out.close();
}

int main() {
    size_t n1 = 2, m1 = 4;
    size_t n2 = 2, m2 = 4;
    //generate_matrix(n1, m1, "matrix1.txt", 6, 0.5);
    //generate_matrix(n2, m2, "matrix2.txt", 6, 0.5);
    PartialBitMatrix matr1 = PartialBitMatrix("matrix1.txt", n1, m1);
    PartialBitMatrix matr2 = PartialBitMatrix("matrix2.txt", n2, m2);
    cout << matr1 << matr2 << endl;
    map<size_t, customset> supporting_rows1, supporting_rows2;

    clock_t start = clock();

    D1_dualization(matr1, matr2, supporting_rows1, supporting_rows2);

    clock_t stop = clock();
    double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
    printf("\nTime elapsed: %.5f s\n", elapsed);
    return 0;
}