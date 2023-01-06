/*
 * SpectraJNIWrapper.java
 *
 */

package spectra;

/*
 * SpectraJNIWrapper
 *
 * @author Xiang Ji
 *
 */

public class SpectraJNIWrapper {

    public SpectraJNIWrapper () {
    }

    static {
        System.load("/usr/local/lib/libspectra-jni.jnilib");
    }

    public native int createInstance(int matrixCount,
                                     int stateCount);

    public native int setMatrix(int matrix,
                                int[] indices,
                                double[] values,
                                int nonZeroCount);

    public native String getVersion();

    public native int getEigenVectors(int matrix,
                                      int numEigenValues,
                                      double sigma,
                                      double[] eigenValues,
                                      double[] eigenVectors);

}