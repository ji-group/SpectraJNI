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

    private SpectraJNIWrapper () {

    }

    public native int createInstance(int matrixCount,
                                     int stateCount);

    public native int setMatrix(int matrix,
                                int[] indices,
                                double[] values,
                                int nonZeroCount);

    public native int getEigenVectors(int matrix,
                                      int numEigenValues,
                                      double sigma,
                                      double[] eigenValues,
                                      double[] eigenVectors);

}