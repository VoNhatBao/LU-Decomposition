#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include "ludcmp.h"

// Function to clear input buffer
void clearInputBuffer() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

// Function to get matrix dimensions
void getDimensions(int& n) {
    while (true) {
        std::cout << "Enter matrix dimension (n x n): ";
        if (std::cin >> n && n > 0) {
            break;
        }
        std::cout << "Invalid input. Please enter a positive integer.\n";
        clearInputBuffer();
    }
}

// Function to input matrix elements
std::vector<std::vector<double>> getMatrix(int n) {
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));

    std::cout << "\nEnter matrix elements row by row:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            while (true) {
                std::cout << "Enter element [" << i + 1 << "][" << j + 1 << "]: ";
                if (std::cin >> matrix[i][j]) {
                    break;
                }
                std::cout << "Invalid input. Please enter a number.\n";
                clearInputBuffer();
            }
        }
    }
    return matrix;
}

// Function to input right-hand side vector
std::vector<double> getRightHandSide(int n) {
    std::vector<double> b(n);

    std::cout << "\nEnter right-hand side vector elements:\n";
    for (int i = 0; i < n; i++) {
        while (true) {
            std::cout << "Enter element [" << i + 1 << "]: ";
            if (std::cin >> b[i]) {
                break;
            }
            std::cout << "Invalid input. Please enter a number.\n";
            clearInputBuffer();
        }
    }
    return b;
}

// Function to display menu and get user choice
int getMenuChoice() {
    int choice;
    std::cout << "\nLU Decomposition Calculator\n";
    std::cout << "1. Solve linear system (Ax = b)\n";
    std::cout << "2. Calculate matrix inverse\n";
    std::cout << "3. Calculate determinant\n";
    std::cout << "4. Exit\n";
    std::cout << "Enter your choice (1-4): ";

    while (!(std::cin >> choice) || choice < 1 || choice > 4) {
        std::cout << "Invalid input. Please enter a number between 1 and 4: ";
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
                // Solve linear system
                auto b = getRightHandSide(n);
                std::vector<double> x(n);
                solver.solve(b, x);

                std::cout << "\nSolution:\n";
                for (int i = 0; i < n; i++) {
                    std::cout << "x[" << i + 1 << "] = " << std::fixed
                        << std::setprecision(6) << x[i] << std::endl;
                }
                break;
            }

            case 2: {
                // Calculate inverse
                std::vector<std::vector<double>> inverse;
                solver.inverse(inverse);

                std::cout << "\nInverse Matrix:\n";
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        std::cout << std::setw(12) << std::fixed
                            << std::setprecision(6) << inverse[i][j];
                    }
                    std::cout << std::endl;
                }
                break;
            }

            case 3: {
                // Calculate determinant
                std::cout << "\nDeterminant = " << std::fixed
                    << std::setprecision(6) << solver.det() << std::endl;
                break;
            }

            case 4:
                std::cout << "Thank you for using LU Decomposition Calculator!\n";
                return 0;
            }
        }

    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
