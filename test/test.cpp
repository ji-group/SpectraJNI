
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/MatOp/DenseGenRealShiftSolve.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <Eigen/Dense>
#include <iostream>



using namespace Spectra;

using DenseOp = DenseGenRealShiftSolve<double>;
using SparseOp = SparseGenRealShiftSolve<double>;
//constexpr bool is_dense = std::is_same<MatType, Matrix>::value;
//using OpType = typename std::conditional<is_dense, DenseOp, SparseOp>::type;

//#include "catch.cpp"

using Eigen::Index;

//TEMPLATE_TEST_CASE("matrix operations [100x100]", "[SparseGenMatProd]", float, double)
int main( int argc, const char* argv[] )
{
std::srand(123);
constexpr Index n = 100;

std::vector<Eigen::Triplet<double>> tripletVector;
for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
    {
//        auto v_ij = dist(gen);
        double v_ij = 1.4 / (i + j + 1.0);
        if (v_ij < 0.5)
        {
            // if larger than treshold, insert it
            tripletVector.push_back(Eigen::Triplet<double>(i, j, v_ij));
        }
    }
Eigen::SparseMatrix<double> mat1(n, n);
// create the matrix
mat1.setFromTriplets(tripletVector.begin(), tripletVector.end());

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat2 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Random(n, n);

SparseGenMatProd<double> sparse1(mat1);
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> xs = sparse1 * mat2;
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ys = mat1 * mat2;

std::cout<<"The matrix-matrix product must be the same as in eigen.";
std::cout<<xs.isApprox(ys);
mat1.coeff(45, 22) == sparse1(45, 22);

    int k = 1;
    int m = 50;
    double sigma = 100.0;

    SparseOp op(mat1);
    GenEigsRealShiftSolver<SparseOp> eigs(op, k, m, sigma);

    eigs.init();
    // maxit = 200 to reduce running time for failed cases
    int nconv = eigs.compute(SortRule::SmallestMagn, 200);
    int niter = eigs.num_iterations();
    int nops = eigs.num_operations();

    std::cout << "nconv = " << nconv << std::endl;
    std::cout << "niter = " << niter << std::endl;
    std::cout << "nops  = " << nops << std::endl;

    Eigen::VectorXcd evals = eigs.eigenvalues();
    Eigen::MatrixXcd evecs = eigs.eigenvectors();

//    Eigen::VectorXcd resid = mat1 * evecs - evecs * evals.asDiagonal();

    std::cout << evals << std::endl;
    std::cout << evecs << std::endl;

return 0;
}
