#ifndef TENSOR_EXAMPLE_H
#define TENSOR_EXAMPLE_H

#include <unsupported/Eigen/CXX11/Tensor>

/**
 * @class TensorExample
 * @brief Example class demonstrating Eigen tensor operations
 * 
 * This class provides examples of creating and manipulating tensors using
 * Eigen's unsupported tensor module, including basic operations and contractions.
 */
class TensorExample {
public:
    /**
     * @brief Default constructor
     */
    TensorExample() = default;
    
    /**
     * @brief Default destructor
     */
    ~TensorExample() = default;
    
    /**
     * @brief Run all tensor examples
     * 
     * Executes both tensor creation and contraction examples,
     * demonstrating various tensor operations.
     */
    void run();

private:
    /**
     * @brief Demonstrates basic tensor creation and initialization
     * 
     * Creates a 3D tensor (3x3x3), initializes it with sequential values,
     * and displays some sample elements to demonstrate tensor indexing.
     */
    void tensor_creation_example();
    
    /**
     * @brief Demonstrates tensor contraction operations
     * 
     * Creates two 2D tensors (matrices), performs tensor contraction
     * (equivalent to matrix multiplication), and displays the results.
     * Shows how to use Eigen's contraction functionality for tensor operations.
     */
    void tensor_contraction_example();
};

#endif // TENSOR_EXAMPLE_H