#include <vector>
#include <stdexcept>
#include <cmath>
using namespace std;

class LUdcmp
{
    private:
        int n;
        vector<vector<double>> lu;
        vector<int> index;
        double d;
        static constexpr double TINY = 1.0e-40;

    public:
        // Constructor
        LUdcmp(const vector<vector<double>>& a): n(static_cast<int>(a.size())),lu(a),index(static_cast<size_t>(n)),d(1.0) 
        {
            //Check square matrix
            if (a.empty() || a[0].size() != static_cast<size_t>(n)) 
                throw runtime_error("LUdcmp: non-square matrix");
            

            vector<double> vv(static_cast<size_t>(n));

            // Loop over rows to get implicit scaling information
            for (int i = 0; i < n; i++) 
            {
                double big = 0.0;
                for (int j = 0; j < n; j++) big = max(big, abs(lu[i][j]));

                if (big == 0.0) throw runtime_error("Singular matrix in LUdcmp");
                vv[i] = 1.0 / big;
            }

            // Crout's algorithm implementation
            //loop each column of matrix
            for (int k = 0; k < n; k++) 
            {
                double big = 0.0;
                int imax = k;
                //loop row from k to n-1
                for (int i = k; i < n; i++) 
                {
                    double temp = vv[i] * abs(lu[i][k]);
                    if (temp > big)
                    {
                        big = temp;
                        imax = i;
                    }
                }

                if (k != imax) 
                {
                    for (int j = 0; j < n; j++) swap(lu[imax][j], lu[k][j]);
                    d = -d;
                    vv[imax] = vv[k];
                }

                index[k] = imax;
                if (lu[k][k] == 0.0) lu[k][k] = TINY;

                for (int i = k + 1; i < n; i++) 
                {
                    //create matrix L
                    lu[i][k] /= lu[k][k];
                    //create matrix U
                    for (int j = k + 1; j < n; j++) lu[i][j] -= lu[i][k] * lu[k][j];
                    
                }
            }
        }
        // Modified solve method for multiple right-hand sides
        void solve(const vector<vector<double>>& b, vector<vector<double>>& x)
        {
            int nrhs = b[0].size();  // Number of right-hand sides

            // Check dimensions
            if (b.size() != n || x.size() != n || (x[0].size() != nrhs)) 
                throw runtime_error("LUdcmp::solve bad sizes");
            

            // Solve each right-hand side
            vector<double> xx(n);
            for (int j = 0; j < nrhs; j++) 
            {
                // Copy the j-th right-hand side
                for (int i = 0; i < n; i++) xx[i] = b[i][j];
                
                // Solve for this right-hand side
                solve(xx, xx);

                // Copy solution to output
                for (int i = 0; i < n; i++)  x[i][j] = xx[i];
                
            }
        }

        // Single right-hand side solver
        void solve(vector<double>& b, vector<double>& x) 
        {
            // Original single vector solve implementation
            if (b.size() != n || x.size() != n) 
                throw runtime_error("LUdcmp::solve bad sizes");
           

            x = b;
            int ii = 0;

            // Forward substitution
            for (int i = 0; i < n; i++) 
            {
                int ip = index[i];
                double sum = x[ip];
                x[ip] = x[i];
                if (ii != 0) 
                {
                    for (int j = ii - 1; j < i; j++) sum -= lu[i][j] * x[j];
                }
                else if (sum != 0.0) ii = i + 1;
                
                x[i] = sum;
            }

            // Back substitution
            for (int i = n - 1; i >= 0; i--) 
            {
                double sum = x[i];
                for (int j = i + 1; j < n; j++) sum -= lu[i][j] * x[j];
                
                x[i] = sum / lu[i][i];
            }
        }

        // Calculate determinant
        double det() 
        {
            double dd = d;
            for (int i = 0; i < n; i++) dd *= lu[i][i];
            return dd;
        }

        // Calculate inverse
        void inverse(vector<vector<double>>& ainv) 
        {
            ainv.resize(n, vector<double>(n, 0.0));
            vector<double> b(n, 0.0);
            vector<double> x(n);

            for (int i = 0; i < n; i++) 
            {
                b[i] = 1.0;
                solve(b, x);
                for (int j = 0; j < n; j++) ainv[j][i] = x[j];
                b[i] = 0.0;
            }
        }
};
