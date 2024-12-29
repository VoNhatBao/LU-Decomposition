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
            if (a.empty() || a[0].size() != static_cast<size_t>(n))  throw runtime_error("LUdcmp: non-square matrix");
               

            vector<double> vv(static_cast<size_t>(n));

            // Loop over rows to get implicit scaling information
            for (int i = 0; i < n; i++) 
            {
                double big = 0.0;
                for (int j = 0; j < n; j++) big = max(big, abs(lu[i][j]));

                if (big == 0.0) throw runtime_error("Singular matrix in LUdcmp");
                vv[i] = 1.0 / big;
            }

            // Crout algorithm implementation
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
                    swap(vv[imax], vv[k]);
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
        void solve(vector<double>& b, vector<double>& x) 
        {
            if (b.size() != n || x.size() != n) throw runtime_error("LUdcmp::solve bad sizes");
            x = b;

            // Ly = b  (Forward substitution)
            for (int i = 0; i < n; i++) 
            {
                // swap row
                swap(x[i], x[index[i]]);
                for (int j = 0; j < i; j++) x[i] -= lu[i][j] * x[j];
            }

            // Ux = y  (Back substitution)
            for (int i = n - 1; i >= 0; i--) 
            {
                for (int j = i + 1; j < n; j++) x[i] -= lu[i][j] * x[j];
                x[i] /= lu[i][i];
            }
        }

        void solve(const vector<vector<double>>& b, vector<vector<double>>& x) 
        {
            int nrhs = b[0].size();  

            if (b.size() != n || x.size() != n || x[0].size() != nrhs) throw runtime_error("LUdcmp::solve bad sizes");
                
            vector<double> single_b(n);
            vector<double> single_x(n);

            for (int k = 0; k < nrhs; k++) 
            {
                // Lấy cột k từ ma trận b
                for (int i = 0; i < n; i++) single_b[i] = b[i][k];
           
                solve(single_b, single_x);

                // Lưu kết quả vào cột k của ma trận x
                for (int i = 0; i < n; i++) x[i][k] = single_x[i];
                
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
        void getLU(vector<vector<double>>& L, vector<vector<double>>& U) const {
            L.resize(n, vector<double>(n, 0.0));
            U.resize(n, vector<double>(n, 0.0));

            // Fill L matrix (including diagonal of 1's)
            for (int i = 0; i < n; i++) {
                L[i][i] = 1.0; // Diagonal elements of L are 1
                for (int j = 0; j < i; j++) {
                    L[i][j] = lu[i][j];
                }
            }

            // Fill U matrix
            for (int i = 0; i < n; i++) {
                for (int j = i; j < n; j++) {
                    U[i][j] = lu[i][j];
                }
            }
        }

        void printLU() const {
            vector<vector<double>> L, U;
            getLU(L, U);

            cout << "\nL Matrix:\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    cout << setw(12) << fixed << setprecision(6) << L[i][j];
                }
                cout << endl;
            }

            cout << "\nU Matrix:\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    cout << setw(12) << fixed << setprecision(6) << U[i][j];
                }
                cout << endl;
            }
        }
};
