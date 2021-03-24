#include<iostream>
#include<vector>
#include<time.h>
#include<math.h>
using namespace std;

namespace MatrixLib {
    class Matrix {

        vector<vector<float>> matrix;

    public:


        Matrix(int rows, int cols, bool isRandom = false, float initVal = 0) {
            matrix.resize(rows, vector<float>(cols, initVal));
            if (isRandom) {
                srand(unsigned(time(0)));
                for (int r = 0; r < rows; ++r) {
                    for (int c = 0; c < cols; ++c) {
                        matrix[r][c] = (float)rand() / 32767.f;
                    }
                }
            }
        }
        
        Matrix(const Matrix& another) {
            this->matrix = another.matrix;
        }

        Matrix(const vector<vector<float>>& vecvec) {
            matrix = vecvec;
        }

        vector<float>& operator[](int row) {
            return matrix[row];
        }

        void operator+=(const Matrix& another) {
            for (int r = 0; r < min(matrix.size(), another.matrix.size()); ++r) {
                for (int c = 0; c < min(matrix[0].size(), another.matrix[0].size()); ++c) {
                    this->matrix[r][c] += another.matrix[r][c];
                }
            }
        }
        void operator+=(float num) {
            for (int r = 0; r < matrix.size(); ++r) {
                for (int c = 0; c < matrix[0].size(); ++c) {
                    this->matrix[r][c] += num;
                }
            }
        }
        void operator-=(const Matrix& another) {
            for (int r = 0; r < min(matrix.size(), another.matrix.size()); ++r) {
                for (int c = 0; c < min(matrix[0].size(), another.matrix[0].size()); ++c) {
                    this->matrix[r][c] -= another.matrix[r][c];
                }
            }
        }
        void operator-=(float num) {
            for (int r = 0; r < matrix.size(); ++r) {
                for (int c = 0; c < matrix[0].size(); ++c) {
                    this->matrix[r][c] -= num;
                }
            }
        }
        void operator*=(float num) {
            for (int r = 0; r < matrix.size(); ++r) {
                for (int c = 0; c < matrix[0].size(); ++c) {
                    this->matrix[r][c] *= num;
                }
            }
        }
        void operator*=(const Matrix& another) {
            matrix = (*this * another).matrix;
        }
        void operator/=(float num) {
            if (num == 0) return;
            for (int r = 0; r < matrix.size(); ++r) {
                for (int c = 0; c < matrix[0].size(); ++c) {
                    this->matrix[r][c] /= num;
                }
            }
        }

        Matrix& operator=(const Matrix& another) {
            this->matrix = another.matrix;
            return *this;
        }
        Matrix operator+(const Matrix& another) {
            Matrix res(*this);
            for (int r = 0; r < min(matrix.size(), another.matrix.size()); ++r) {
                for (int c = 0; c < min(matrix[0].size(), another.matrix[0].size()); ++c) {
                    res.matrix[r][c] += another.matrix[r][c];
                }
            }
            return res;
        }
        Matrix operator+(float num) {
            Matrix res(*this);
            for (int r = 0; r < matrix.size(); ++r) {
                for (int c = 0; c < matrix[0].size(); ++c) {
                    res.matrix[r][c] += num;
                }
            }
            return res;
        }
        Matrix operator-(const Matrix& another) {
            Matrix res(*this);
            for (int r = 0; r < min(matrix.size(), another.matrix.size()); ++r) {
                for (int c = 0; c < min(matrix[0].size(), another.matrix[0].size()); ++c) {
                    res.matrix[r][c] -= another.matrix[r][c];
                }
            }
            return res;
        }
        Matrix operator-(float num) {
            return *this + (-num);
        }
        Matrix operator*(float num) {
            Matrix res(*this);
            for (int r = 0; r < matrix.size(); ++r) {
                for (int c = 0; c < matrix[0].size(); ++c) {
                    res.matrix[r][c] *= num;
                }
            }
            return res;
        }
        Matrix operator*(const Matrix& another) {
            int R1 = matrix.size();
            if (R1 == 0) return *this;
            int C1 = matrix[0].size();
            int R2 = another.matrix.size();
            if (R2 == 0) return *this;
            int C2 = another.matrix[0].size();

            if (C1 != R2) return *this;

            Matrix res(R1, C2);

            for (int r = 0; r < R1; ++r) {
                for (int c = 0; c < C2; ++c) {
                    double sum = 0;
                    for (int i = 0; i < C1; ++i) {
                        sum += this->matrix[r][i] * another.matrix[i][c];
                    }

                    res[r][c] = sum;
                }
            }

            return res;
        }   
        Matrix operator/(float num) {
            if (num == 0) return *this;
            return *this * (1.f / num);
        }

        vector<float>& getRows(int row) {
            return matrix[row];
        }
        vector<float> getCols(int col) {
            vector<float> colVec; colVec.reserve(matrix.size());
            for (int r = 0; r < matrix.size(); ++r) {
                colVec.emplace_back(matrix[r][col]);
            }
            return colVec;
        }

        void addRow(vector<float>& rowVec, int row = -1) {
            if (row == -1) row = matrix.size();
            auto it = matrix.insert(matrix.begin() + row, rowVec);
            if (matrix.size()) it->resize(matrix[0].size(), 0);

        }
        void addCol(vector<float>& colVec, int col = -1) {
            if (matrix.empty()) {
                for (int num : colVec) {
                    matrix.push_back({});
                    matrix.back().emplace_back(num);
                }
                return;
            }
            if (col == -1) col = matrix[0].size();

            for (int r = 0; r < matrix.size(); ++r) {
                matrix[r].insert(matrix[r].begin()+col, r < colVec.size() ? colVec[r] : 0);
            }

        }
        void setRow(vector<float>& rowVec, int row) {
            for (int c = 0; c < matrix[0].size(); ++c) {
                matrix[row][c] = (c < rowVec.size() ? rowVec[c] : 0);
            }
        }
        void setCol(vector<float>& colVec, int col) {
            for (int r = 0; r < matrix.size(); ++r) {
                matrix[r][col] = (r < colVec.size() ? colVec[r] : 0);
            }
        }
        void extendRow(int amount, int initVal = 0) {
            if (matrix.empty()) return;
            while (amount-- > 0) matrix.emplace_back(vector<float>(matrix[0].size(), initVal));
        }
        void extendCol(int amount, int initVal = 0) {
            if (matrix.empty()) return;
            for (int r = 0; r < matrix.size(); ++r) {
                matrix[r].resize(matrix[r].size() + amount, initVal);
            }
        }
        void shrinkRow(int amount) {
            while (amount-- > 0) matrix.pop_back();
        }
        void shrinkCol(int amount) {
            for (int r = 0; r < matrix.size(); ++r) {
                for (int t = 0; t < amount; ++t) matrix[r].pop_back();
            }
        }


        float norm() const {
            double res = 0;
            for (int r = 0; r < matrix.size(); ++r) {
                for (int c = 0; c < matrix[0].size(); ++c) {
                    res += matrix[r][c] * matrix[r][c];
                }
            }
            return sqrt(res);
        }
        void transpose() {
            vector<vector<float>> res(matrix[0].size(), vector<float>(matrix.size()));
            for (int r = 0; r < matrix.size(); ++r) {
                for (int c = 0; c < matrix[0].size(); ++c) {
                    res[c][r] = matrix[r][c];
                }
            }
            matrix = res;
        }
        float dot(const Matrix& another) const {
            double res = 0;
            for (int r = 0; r < min(matrix.size(), another.matrix.size()); ++r) {
                for (int c = 0; c < min(matrix[0].size(), another.matrix[0].size()); ++c) {
                    res += another.matrix[r][c] * this->matrix[r][c];
                }
            }
            return res;
        }


        int rowSize() const {
            return matrix.size();
        }
        int colSize() const {
            return rowSize() ? matrix[0].size() : 0;
        }
        float at(int row, int col) const {
            return matrix[row][col];
        }


        void print() const {
            cout << "Matrix===========================================" << endl;
            cout << "Rows = " << matrix.size() << " , Cols = " << (matrix.size() ? matrix[0].size() : 0) << endl;
            if (!matrix.size()) return;
            cout << '\\' << '\t';
            for (int c = 0; c < matrix[0].size(); ++c) {
                cout << c << '\t';
            }
            cout << endl;

            for (int r = 0; r < matrix.size(); ++r) {

                int count = 0;
                for (int c = 0; c < matrix[0].size(); ++c) {
                    if (count == 0) {
                        cout << r << '\t';
                        ++count;
                    }
                    cout << matrix[r][c] << '\t';
                }
                cout << endl;
            }

            cout << endl;
        }
    };


    Matrix unitMatrix(int n) {
        Matrix res(n, n);
        for (int i = 0; i < n; ++i) {
            res[i][i] = 1;
        }
        return res;
    }

    float norm(const Matrix& m) {
        return m.norm();
    }

    float dot(const Matrix& m1, const Matrix& m2) {
        return m1.dot(m2);
    }

    Matrix transpose(const Matrix& m) {
        Matrix res(m.colSize(), m.rowSize());
        for (int r = 0; r < m.rowSize(); ++r) {
            for (int c = 0; c < m.colSize(); ++c) {
                res[c][r] = m.at(r, c);
            }
        }
        return res;
    }


};




int main () {

    MatrixLib::Matrix m1({
        {-2, 0, 3}
    });
    cout << m1.norm() << endl;
    

    m1.print();

    system("pause");
}