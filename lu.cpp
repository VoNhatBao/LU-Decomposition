#include <iostream>
#include "ludcmp.h"
int main() {
    // Test matrix 2x2
    std::vector<std::vector<double>> A = {
        {4.0, 3.0},
        {6.0, 3.0}
    };

    // Create right-hand side
    std::vector<std::vector<double>> b = {
        {10.0, 5.0},
        {12.0, 6.0}
    };

    try {
        LUdcmp solver(A);

        // Prepare solution matrix
        std::vector<std::vector<double>> x(2, std::vector<double>(2));

        // Solve system
        solver.solve(b, x);

        // Print solution
        std::cout << "Solution:" << std::endl;
        for (const auto& row : x) {
            for (double val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }

    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}