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
    LUdcmp(const std::vector<std::vector<double>>& a) : n(a.size()), lu(a), indx(n), d(1.0) {
        // Constructor implementation remains the same
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