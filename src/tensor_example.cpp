#include "../include/tensors/tensor_example.h"
#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>

void TensorExample::run() {
    std::cout << "Running tensor examples...\n";
    tensor_creation_example();
    tensor_contraction_example();
}

void TensorExample::tensor_creation_example() {
    // Create a 3x3x3 tensor
    Eigen::Tensor<double, 3> tensor(3, 3, 3);
    
    // Initialize with sequential values
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                tensor(i, j, k) = i + j * 3 + k * 9;
            }
        }
    }
    
    // Print some values
    std::cout << "\n3D Tensor example:\n";
    std::cout << "tensor(0,0,0) = " << tensor(0,0,0) << "\n";
    std::cout << "tensor(1,1,1) = " << tensor(1,1,1) << "\n";
    std::cout << "tensor(2,2,2) = " << tensor(2,2,2) << "\n";
}

void TensorExample::tensor_contraction_example() {
    // Create two 2D tensors (matrices)
    Eigen::Tensor<double, 2> a(2, 2);
    Eigen::Tensor<double, 2> b(2, 2);
    
    // Initialize
    a.setValues({{1, 2}, {3, 4}});
    b.setValues({{5, 6}, {7, 8}});
    
    // Contract over second dimension of a and first dimension of b
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair = {Eigen::IndexPair<int>(1, 0)};
    Eigen::Tensor<double, 2> result = a.contract(b, contraction_pair);
    
    std::cout << "\nTensor contraction example (matrix multiplication):\n";
    std::cout << "Matrix A:\n" << a << "\n";
    std::cout << "Matrix B:\n" << b << "\n";
    std::cout << "Result (A * B):\n" << result << "\n";
}