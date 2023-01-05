
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
    constexpr Index n = 12;

    std::vector<Eigen::Triplet<double>> tripletVector;
    tripletVector.push_back(Eigen::Triplet<double>(0, 0, -1.0330024148108399));
    tripletVector.push_back(Eigen::Triplet<double>(0, 3, 0.4024684733029246));
    tripletVector.push_back(Eigen::Triplet<double>(0, 4, 0.2817279313120472));
    tripletVector.push_back(Eigen::Triplet<double>(0, 5, 0.34880601019586804));

    tripletVector.push_back(Eigen::Triplet<double>(1, 1, -1.0598336463643683));
    tripletVector.push_back(Eigen::Triplet<double>(1, 6, 0.4024684733029246));
    tripletVector.push_back(Eigen::Triplet<double>(1, 7, 0.30855916286557555));
    tripletVector.push_back(Eigen::Triplet<double>(1, 8, 0.34880601019586804));

    tripletVector.push_back(Eigen::Triplet<double>(2, 2, -0.9927555674805474));
    tripletVector.push_back(Eigen::Triplet<double>(2, 9, 0.4024684733029246));
    tripletVector.push_back(Eigen::Triplet<double>(2, 10, 0.30855916286557555));
    tripletVector.push_back(Eigen::Triplet<double>(2, 11, 0.2817279313120472));

    tripletVector.push_back(Eigen::Triplet<double>(3, 0, 0.30855916286557555));
    tripletVector.push_back(Eigen::Triplet<double>(3, 1, 0.2817279313120472));
    tripletVector.push_back(Eigen::Triplet<double>(3, 2, 0.34880601019586804));
    tripletVector.push_back(Eigen::Triplet<double>(3, 3, -0.9390931043734909));

    tripletVector.push_back(Eigen::Triplet<double>(4, 4, -1.0598336463643683));
    tripletVector.push_back(Eigen::Triplet<double>(4, 6, 0.4024684733029246));
    tripletVector.push_back(Eigen::Triplet<double>(4, 7, 0.30855916286557555));
    tripletVector.push_back(Eigen::Triplet<double>(4, 8, 0.34880601019586804));

    tripletVector.push_back(Eigen::Triplet<double>(5, 5, -0.9927555674805474));
    tripletVector.push_back(Eigen::Triplet<double>(5, 9, 0.4024684733029246));
    tripletVector.push_back(Eigen::Triplet<double>(5, 10, 0.30855916286557555));
    tripletVector.push_back(Eigen::Triplet<double>(5, 11, 0.2817279313120472));

    tripletVector.push_back(Eigen::Triplet<double>(6, 0, 0.30855916286557555));
    tripletVector.push_back(Eigen::Triplet<double>(6, 1, 0.2817279313120472));
    tripletVector.push_back(Eigen::Triplet<double>(6, 2, 0.34880601019586804));
    tripletVector.push_back(Eigen::Triplet<double>(6, 6, -0.9390931043734909));

    tripletVector.push_back(Eigen::Triplet<double>(7, 3, 0.4024684733029246));
    tripletVector.push_back(Eigen::Triplet<double>(7, 4, 0.2817279313120472));
    tripletVector.push_back(Eigen::Triplet<double>(7, 5, 0.34880601019586804));
    tripletVector.push_back(Eigen::Triplet<double>(7, 7, -1.0330024148108399));

    tripletVector.push_back(Eigen::Triplet<double>(8, 8, -0.9927555674805474));
    tripletVector.push_back(Eigen::Triplet<double>(8, 9, 0.4024684733029246));
    tripletVector.push_back(Eigen::Triplet<double>(8, 10, 0.30855916286557555));
    tripletVector.push_back(Eigen::Triplet<double>(8, 11, 0.2817279313120472));

    tripletVector.push_back(Eigen::Triplet<double>(9, 0, 0.30855916286557555));
    tripletVector.push_back(Eigen::Triplet<double>(9, 1, 0.2817279313120472));
    tripletVector.push_back(Eigen::Triplet<double>(9, 2, 0.34880601019586804));
    tripletVector.push_back(Eigen::Triplet<double>(9, 9, -0.9390931043734909));

    tripletVector.push_back(Eigen::Triplet<double>(10, 3, 0.4024684733029246));
    tripletVector.push_back(Eigen::Triplet<double>(10, 4, 0.2817279313120472));
    tripletVector.push_back(Eigen::Triplet<double>(10, 5, 0.34880601019586804));
    tripletVector.push_back(Eigen::Triplet<double>(10, 10, -1.0330024148108399));

    tripletVector.push_back(Eigen::Triplet<double>(11, 6, 0.4024684733029246));
    tripletVector.push_back(Eigen::Triplet<double>(11, 7, 0.30855916286557555));
    tripletVector.push_back(Eigen::Triplet<double>(11, 8, 0.34880601019586804));
    tripletVector.push_back(Eigen::Triplet<double>(11, 11, -1.0598336463643683));





//    for (int i = 0; i < n; ++i)
//    for (int j = 0; j < n; ++j)
//    {
////        auto v_ij = dist(gen);
//        double v_ij = 1.4 / (i + j + 1.0);
//        if (v_ij < 0.5)
//        {
//            // if larger than treshold, insert it
//            tripletVector.push_back(Eigen::Triplet<double>(i, j, v_ij));
//        }
//    }
    Eigen::SparseMatrix<double> mat1(n, n);
// create the matrix
    mat1.setFromTriplets(tripletVector.begin(), tripletVector.end());

    mat1 = mat1.transpose();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat2 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Random(n, n);

    SparseGenMatProd<double> sparse1(mat1);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> xs = sparse1 * mat2;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ys = mat1 * mat2;

    std::cout<<"Matrix =" << std::endl;
    std::cout<<mat1<<std::endl;

    int k = n - 2;
    k = 1;
    int m = k + 2;
    double sigma = 1E-7;

    SparseOp op(mat1);
    GenEigsRealShiftSolver<SparseOp> eigs(op, k, m, sigma);

    eigs.init();
    // maxit = 200 to reduce running time for failed cases
    int nconv = eigs.compute(SortRule::SmallestReal, 200);
    int niter = eigs.num_iterations();
    int nops = eigs.num_operations();

    std::cout << "nconv = " << nconv << std::endl;
    std::cout << "niter = " << niter << std::endl;
    std::cout << "nops  = " << nops << std::endl;

    Eigen::VectorXcd evals = eigs.eigenvalues();

    std::cout << "Eigen values:\n" << evals << std::endl;

    Eigen::MatrixXcd evecs = eigs.eigenvectors();

//    Eigen::VectorXcd resid = mat1 * evecs - evecs * evals.asDiagonal();


    std::cout << "Eigen vectors:\n" << evecs << std::endl;

    return 0;
}
