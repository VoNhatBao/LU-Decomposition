#include <vector>
#include <stdexcept>
#include <cmath>

class LUdcmp {
private:
    int n;
    std::vector<std::vector<double>> lu;
    std::vector<int> indx;
    double d;
    static constexpr double TINY = 1.0e-40;

public:
    // Constructor
    LUdcmp(const std::vector<std::vector<double>>& a)
        : n(static_cast<int>(a.size())),
        lu(a),
        indx(static_cast<size_t>(n)),
        d(1.0) {

        if (a.empty() || a[0].size() != static_cast<size_t>(n)) {
            throw std::runtime_error("LUdcmp: non-square matrix");
        }

        std::vector<double> vv(static_cast<size_t>(n));

        // Loop over rows to get implicit scaling information
        for (int i = 0; i < n; i++) {
            double big = 0.0;
            for (int j = 0; j < n; j++) {
                double temp = std::abs(lu[i][j]);
                if (temp > big) big = temp;
            }
            if (big == 0.0) throw std::runtime_error("Singular matrix in LUdcmp");
            vv[i] = 1.0 / big;
        }

        // Crout's algorithm implementation
        for (int k = 0; k < n; k++) {
            double big = 0.0;
            int imax = k;

            for (int i = k; i < n; i++) {
                double temp = vv[i] * std::abs(lu[i][k]);
                if (temp > big) {
                    big = temp;
                    imax = i;
                }
            }

            if (k != imax) {
                for (int j = 0; j < n; j++) {
                    double temp = lu[imax][j];
                    lu[imax][j] = lu[k][j];
                    lu[k][j] = temp;
                }
                d = -d;
                vv[imax] = vv[k];
            }

            indx[k] = imax;
            if (lu[k][k] == 0.0) lu[k][k] = TINY;

            for (int i = k + 1; i < n; i++) {
                lu[i][k] /= lu[k][k];
                for (int j = k + 1; j < n; j++) {
                    lu[i][j] -= lu[i][k] * lu[k][j];
                }
            }
        }
    }
    // Modified solve method for multiple right-hand sides
    void solve(const std::vector<std::vector<double>>& b, std::vector<std::vector<double>>& x) {
        int nrhs = b[0].size();  // Number of right-hand sides

        // Check dimensions
        if (b.size() != n || x.size() != n || (x[0].size() != nrhs)) {
            throw std::runtime_error("LUdcmp::solve bad sizes");
        }

        // Solve each right-hand side
        std::vector<double> xx(n);
        for (int j = 0; j < nrhs; j++) {
            // Copy the j-th right-hand side
            for (int i = 0; i < n; i++) {
                xx[i] = b[i][j];
            }

            // Solve for this right-hand side
            solve(xx, xx);

            // Copy solution to output
            for (int i = 0; i < n; i++) {
                x[i][j] = xx[i];
            }
        }
    }

    // Single right-hand side solver
    void solve(std::vector<double>& b, std::vector<double>& x) {
        // Original single vector solve implementation
        if (b.size() != n || x.size() != n) {
            throw std::runtime_error("LUdcmp::solve bad sizes");
        }

        x = b;
        int ii = 0;

        // Forward substitution
        for (int i = 0; i < n; i++) {
            int ip = indx[i];
            double sum = x[ip];
            x[ip] = x[i];
            if (ii != 0) {
                for (int j = ii - 1; j < i; j++) {
                    sum -= lu[i][j] * x[j];
                }
            }
            else if (sum != 0.0) {
                ii = i + 1;
            }
            x[i] = sum;
        }

        // Back substitution
        for (int i = n - 1; i >= 0; i--) {
            double sum = x[i];
            for (int j = i + 1; j < n; j++) {
                sum -= lu[i][j] * x[j];
            }
            x[i] = sum / lu[i][i];
        }
    }

    // Calculate determinant
    double det() {
        double dd = d;
        for (int i = 0; i < n; i++) {
            dd *= lu[i][i];
        }
        return dd;
    }

    // Calculate inverse
    void inverse(std::vector<std::vector<double>>& ainv) {
        ainv.resize(n, std::vector<double>(n, 0.0));
        std::vector<double> b(n, 0.0);
        std::vector<double> x(n);

        for (int i = 0; i < n; i++) {
            b[i] = 1.0;
            solve(b, x);
            for (int j = 0; j < n; j++) {
                ainv[j][i] = x[j];
            }
            b[i] = 0.0;
        }
    }
};
