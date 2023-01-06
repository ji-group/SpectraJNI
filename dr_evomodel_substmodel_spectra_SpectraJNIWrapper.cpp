#include "dr_evomodel_substmodel_spectra_SpectraJNIWrapper.h"
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/MatOp/DenseGenRealShiftSolve.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <Eigen/Dense>
#include <iostream>

typedef Eigen::SparseMatrix<double> SpMatrix;
using SparseOp = Spectra::SparseGenRealShiftSolve<double>;

class SpectraImpl {

private:
    SpMatrix* matrices;
    int kState;
    enum SpectraReturnCodes {
        SPECTRA_SUCCESS = 0,
    };

public:
    int createInstance(int matrixCount,
                       int stateCount) {
        matrices = new SpMatrix [matrixCount];
        kState = stateCount;
        for (int i = 0; i < matrixCount; i++) {
            SpMatrix matrix(stateCount, stateCount);
            matrices[i] = matrix;
        }
        return SPECTRA_SUCCESS;
    };

    int setMatrix(int matrix,
                  int* indices,
                  double* values,
                  int nonZeroCount){
        std::vector<Eigen::Triplet<double>> tripletVector;
        for (int i = 0; i < nonZeroCount; i++) {
            tripletVector.push_back(Eigen::Triplet<double>(indices[2 * i], indices[2 * i + 1], values[i]));
        }
        matrices[matrix].setFromTriplets(tripletVector.begin(), tripletVector.end());
        return SPECTRA_SUCCESS;
    };

    int getEigenVectors(int matrix,
                        int numEigenValues,
                        double sigma,
                        double* eigenValues,
                        double* eigenVectors) {
        SparseOp op(matrices[matrix]);
        Spectra::GenEigsRealShiftSolver<SparseOp> eigs(op, numEigenValues, numEigenValues + 2, sigma);
        eigs.init();
        int nconv = eigs.compute(Spectra::SortRule::SmallestReal, 200);
        int niter = eigs.num_iterations();
        int nops = eigs.num_operations();

        Eigen::VectorXcd evals = eigs.eigenvalues();
        Eigen::MatrixXcd evecs = eigs.eigenvectors();

        memcpy(eigenValues, evals.data(), 2 * numEigenValues);
        memcpy(eigenVectors, evecs.data(), 2 * numEigenValues * kState);
        return SPECTRA_SUCCESS;
    }

    static SpectraImpl* factory(int matrixCount,
                                int stateCount,
                                int* errorCode) {
        SpectraImpl* impl = new SpectraImpl();
        try {
            *errorCode = impl->createInstance(matrixCount, stateCount);
            if (*errorCode == SPECTRA_SUCCESS) {
                return impl;
            }
            delete impl;
            return NULL;
        }
        catch(...) {
            delete impl;
            throw;
        }

        delete impl;
        return NULL;
    };


    ~SpectraImpl() {
        free(matrices);
    };

};

typedef std::shared_ptr<SpectraImpl> InstancePtr;
std::vector<InstancePtr> instances;

/*
 * Class:     dr_evomodel_substmodel_spectra_SpectraJNIWrapper
 * Method:    createInstance
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_dr_evomodel_substmodel_spectra_SpectraJNIWrapper_createInstance
        (JNIEnv * env, jobject obj, jint matrixCount, jint stateCount) {
    jint errorCode;
    instances.emplace_back(
            SpectraImpl::factory(matrixCount, stateCount, &errorCode)
    );
    return errorCode;
}


/*
 * Class:     dr_evomodel_substmodel_spectra_SpectraJNIWrapper
 * Method:    setMatrix
 * Signature: (I[I[DI)I
 */
JNIEXPORT jint JNICALL Java_dr_evomodel_substmodel_spectra_SpectraJNIWrapper_setMatrix
        (JNIEnv * env, jobject obj, jint matrix, jintArray inIndices, jdoubleArray inValues, jint nonZeroCount) {
    jint *indices = env->GetIntArrayElements(inIndices, NULL);
    jdouble *values = env->GetDoubleArrayElements(inValues, NULL);

    jint errCode = (jint) instances[0]->setMatrix(matrix, indices, values, nonZeroCount);

    env->ReleaseDoubleArrayElements(inValues, values, JNI_ABORT);
    env->ReleaseIntArrayElements(inIndices, indices, JNI_ABORT);
    return errCode;
}


/*
 * Class:     dr_evomodel_substmodel_spectra_SpectraJNIWrapper
 * Method:    getEigenVectors
 * Signature: (IID[D[D)I
 */
JNIEXPORT jint JNICALL Java_dr_evomodel_substmodel_spectra_SpectraJNIWrapper_getEigenVectors
        (JNIEnv *env, jobject obj, jint matrix, jint numEigenValues, jdouble sigma, jdoubleArray outEigenValues, jdoubleArray outEigenVectors) {
    jdouble *eigenValues = env->GetDoubleArrayElements(outEigenValues, NULL);
    jdouble *eigenVectors = env->GetDoubleArrayElements(outEigenVectors, NULL);

    jint errCode = (jint) instances[0]->getEigenVectors(matrix, numEigenValues, sigma, eigenValues, eigenVectors);

    env->ReleaseDoubleArrayElements(outEigenValues, eigenValues, 0);
    env->ReleaseDoubleArrayElements(outEigenVectors, eigenVectors, 0);

    return errCode;
}

/*
 * Class:     dr_evomodel_substmodel_spectra_SpectraJNIWrapper
 * Method:    getVersion
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_dr_evomodel_substmodel_spectra_SpectraJNIWrapper_getVersion
        (JNIEnv *env, jobject obj) {
    return env->NewStringUTF("0.0.1");
}

