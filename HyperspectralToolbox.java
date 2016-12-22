import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.Arrays;

/**
 * Created by ARPAN on 12-11-2016.
 */
class Point{
    double x, y;

    Point(double x,double y){
        this.x=x;
        this.y=y;
    }

    public String toString() {
        return "("+x + "," + y+")";
    }

}

public class HyperspectralToolbox {

    public static double cross(Point O, Point A, Point B) {
        return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
    }

    public static Point[] convexHull(Point[] P) {
        if (P.length > 1) {
            int n = P.length, k = 0;
            Point[] H = new Point[2 * n];

            //Arrays.sort(P);

            // Build upper hull
            for (int i = 0; i <n; i++) {
                while (k >= 2 && cross(H[k - 2], H[k - 1], P[i]) >= 0)
                    k--;
                H[k++] = P[i];
            }

            //Last element must be a part of convex hull
            if(H[k]!=P[P.length-1])
                H[k++]=P[P.length-1];

            if (k > 1) {
                H = Arrays.copyOfRange(H, 0, k - 1); // remove non-hull vertices after k; remove k - 1 which is a duplicate
            }
            return H;
        } else if (P.length <= 1) {
            return P;
        } else{
            return null;
        }
    }

    public static double[] interpLinear(double[] x, double[] y, double[] xi) throws IllegalArgumentException {

        if (x.length != y.length) {
            throw new IllegalArgumentException("X and Y must be the same length");
        }
        if (x.length == 1) {
            throw new IllegalArgumentException("X must contain more than one value");
        }
        double[] dx = new double[x.length - 1];
        double[] dy = new double[x.length - 1];
        double[] slope = new double[x.length - 1];
        double[] intercept = new double[x.length - 1];

        // Calculate the line equation (i.e. slope and intercept) between each point
        for (int i = 0; i < x.length - 1; i++) {
            dx[i] = x[i + 1] - x[i];
            if (dx[i] == 0) {
                throw new IllegalArgumentException("X must be montotonic. A duplicate " + "x-value was found");
            }
            if (dx[i] < 0) {
                throw new IllegalArgumentException("X must be sorted");
            }
            dy[i] = y[i + 1] - y[i];
            slope[i] = dy[i] / dx[i];
            intercept[i] = y[i] - x[i] * slope[i];
        }

        // Perform the interpolation here
        double[] yi = new double[xi.length];
        for (int i = 0; i < xi.length; i++) {
            if ((xi[i] > x[x.length - 1]) || (xi[i] < x[0])) {
                yi[i] = Double.NaN;
            }
            else {
                int loc = Arrays.binarySearch(x, xi[i]);
                if (loc < -1) {
                    loc = -loc - 2;
                    yi[i] = slope[loc] * xi[i] + intercept[loc];
                }
                else {
                    yi[i] = y[loc];
                }
            }
        }

        return yi;
    }

    public static double[] removeContinuum(Point[] data,Point[] convexHull,int bands){
        double[] x=new double[convexHull.length];
        double[] y=new double[convexHull.length];

        for(int i=0;i<convexHull.length;i++){
            x[i]=convexHull[i].x;
            y[i]=convexHull[i].y;
        }

        double[] xi=new double[bands];
        double[] yi=new double[bands];
        for(int i=0;i<bands;i++){
            xi[i]=data[i].x;
        }

        yi=interpLinear(x,y,xi);
        double[] CR=new double[bands];

        for(int i=0;i<bands;i++){
            CR[i]=data[i].y/yi[i];

            if(CR[i]>1)
                CR[i]=1;
        }
        return CR;
    }

    public static double[] performContinuumRemoval(double[] data,double[] wavelength,int nDim){
        Point[] pixel=new Point[nDim];
        double[] CRdata=new double[nDim];


            for (int i = 0; i < nDim; i++) {
                pixel[i]=new Point(wavelength[i],data[i]);
            }

            Point[] hull = convexHull(pixel).clone();
            CRdata=removeContinuum(pixel,hull,nDim);


        return CRdata;
    }

    public static double[][] runningAvgFilter(double[][] CRdata,int nData,int nDim,int pt){
        int newDim=nDim-pt+1;
        double[][] CRSdata=new double[nData][newDim];
        double sum;

        for(int i=0;i<nData;i++){
            for(int j=0;j<newDim;j++){
                sum=0;
                //System.out.println(i+" "+j+"\t"+CRdata[i][j]);
                for(int k=j;k<j+pt;k++){
                    sum+=CRdata[i][k];
                }
                //System.out.println(sum);
                CRSdata[i][j]=sum/pt;
                //System.out.println(i+" "+j+"\t"+CRSdata[i][j]);

                if(CRSdata[i][j]>1)
                    CRSdata[i][j]=1;

                CRSdata[i][j]=1-CRSdata[i][j];

            }
        }
        return CRSdata;
    }

    public static double[][] SAMIntensityImage(double[][] data,double[] librarySpectra,int nData, int nDim,int imgDim1,int imgDim2)
    throws IllegalArgumentException{
        if(data[0].length != librarySpectra.length){
            throw new IllegalArgumentException("Dimensions of data and library spectra must match");
        }

        RealMatrix Data=new Array2DRowRealMatrix(data);
        RealVector LibrarySpectra=new ArrayRealVector(librarySpectra);
        double[] SAMvec=new double[nData];

        for(int i=0;i<nData;i++){
            SAMvec[i]=Data.getRowVector(i).cosine(LibrarySpectra);
        }

        return reshape(SAMvec,imgDim1,imgDim2);

    }

    public static double[] invert(double[] data){
        for(int i=0;i<data.length;i++){
            if(data[i]>1){
                data[i]=1;
            }
            else if(data[i]<0){
                data[i]=0;
            }
            data[i]=1-data[i];
        }
        return  data;
    }

    public static int[][] reshape(int[] index,int imgDim1,int imgDim2){
        int count=0;
        int[][] classificationMat=new int[imgDim1][imgDim2];
        for(int j=0;j<imgDim2;j++)
            for(int i=0;i<imgDim1;i++){
                classificationMat[i][j]=index[count];count++;
            }
        return classificationMat;
    }

    public static double[][] reshape(double[] index,int imgDim1,int imgDim2){
        int count=0;
        double[][] classificationMat=new double[imgDim1][imgDim2];
        for(int j=0;j<imgDim2;j++)
            for(int i=0;i<imgDim1;i++){
                classificationMat[i][j]=index[count];count++;
            }
        return classificationMat;
    }
}
