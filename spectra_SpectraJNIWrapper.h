/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class spectra_SpectraJNIWrapper */

#ifndef _Included_spectra_SpectraJNIWrapper
#define _Included_spectra_SpectraJNIWrapper
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     spectra_SpectraJNIWrapper
 * Method:    createInstance
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_spectra_SpectraJNIWrapper_createInstance
  (JNIEnv *, jobject, jint, jint);

/*
 * Class:     spectra_SpectraJNIWrapper
 * Method:    setMatrix
 * Signature: (I[I[DI)I
 */
JNIEXPORT jint JNICALL Java_spectra_SpectraJNIWrapper_setMatrix
  (JNIEnv *, jobject, jint, jintArray, jdoubleArray, jint);

/*
 * Class:     spectra_SpectraJNIWrapper
 * Method:    getEigenVectors
 * Signature: (IID[D[D)I
 */
JNIEXPORT jint JNICALL Java_spectra_SpectraJNIWrapper_getEigenVectors
  (JNIEnv *, jobject, jint, jint, jdouble, jdoubleArray, jdoubleArray);

#ifdef __cplusplus
}
#endif
#endif
