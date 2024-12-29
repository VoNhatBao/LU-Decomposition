#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include "ludcmp.h"

void clearInputBuffer() {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
}

void getDimensions(int& n) {
    while (true) {
        cout << "Enter matrix dimension (n x n): ";
        if (cin >> n && n > 0) {
            clearInputBuffer();
            break;
        }
        cout << "Invalid input. Please enter a positive integer.\n";
        clearInputBuffer();
    }
}

vector<vector<double>> getMatrix(int n) {
    vector<vector<double>> matrix(n, vector<double>(n));
    cout << "\nEnter matrix elements row by row (separate elements by spaces, press Enter after each row):\n";

    for (int i = 0; i < n; i++) {
        while (true) {
            string line;
            getline(cin, line);
            istringstream iss(line);
            vector<double> row;
            double value;

            
            while (iss >> value) {
                row.push_back(value);
            }

           
            if (row.size() != n) {
                cout << "Error: Please enter exactly " << n << " numbers for row " << i + 1 << ".\n";
                continue;
            }

           
            for (int j = 0; j < n; j++) {
                matrix[i][j] = row[j];
            }
            break;
        }
    }
    return matrix;
}

vector<double> getRightHandSide(int n) {
    vector<double> b(n);
    cout << "\nEnter right-hand side vector elements (separate by spaces): ";

    while (true) {
        string line;
        clearInputBuffer();
        getline(cin, line);
        istringstream iss(line);
        vector<double> values;
        double value;

        while (iss >> value) {
            values.push_back(value);
        }

        if (values.size() != n) {
            cout << "Error: Please enter exactly " << n << " numbers.\n";
            continue;
        }

        b = values;
        break;
    }
    return b;
}


int getMenuChoice() {
    int choice;
    cout << "\nLU Decomposition Calculator\n";
    cout << "1. Solve linear system (Ax = b)\n";
    cout << "2. Calculate matrix inverse\n";
    cout << "3. Calculate determinant\n";
    cout << "4. Exit\n";
    cout << "5. Print LU decomposition" << "\n";
    cout << "Enter your choice (1-4): ";
    while (!(cin >> choice) || choice < 1 || choice > 5) {
        cout << "Invalid input. Please enter a number between 1 and 4: ";
        clearInputBuffer();
    }
    return choice;
}

int main() {
    try {
        int n;
        getDimensions(n);
        auto A = getMatrix(n);
        LUdcmp solver(A);
        
        while (true) {
            int choice = getMenuChoice();
            switch (choice) {
            case 1: {
                auto b = getRightHandSide(n);
                vector<double> x(n);
                solver.solve(b, x);
                cout << "\nSolution:\n";
                for (int i = 0; i < n; i++) {
                    cout << "x[" << i + 1 << "] = " << fixed
                        << setprecision(6) << x[i] << endl;
                }
                break;
            }
            case 2: {
                vector<vector<double>> inverse;
                solver.inverse(inverse);
                cout << "\nInverse Matrix:\n";
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        cout << setw(12) << fixed
                            << setprecision(6) << inverse[i][j];
                    }
                    cout << endl;
                }
                break;
            }
            case 3: {
                cout << "\nDeterminant = " << fixed
                    << setprecision(6) << solver.det() << endl;
                break;
            }
            case 4:{
                cout << "Thank you for using LU Decomposition Calculator!\n";
                return 0;
            }
            case 5:
                solver.printLU();
                break;
            }
        }
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
