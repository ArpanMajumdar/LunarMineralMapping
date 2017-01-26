

/**
 * Created by ARPAN on 29-10-2016.
 */
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Scanner;
import org.apache.commons.math3.stat.StatUtils;

public class BandParameterDetermination {
    public static double[] scale(double[] data){
        double minVal=StatUtils.min(data),maxVal=StatUtils.max(data);
        double[] scaledData=new double[data.length];
        for(int i=0;i<data.length;i++){
            scaledData[i]=(data[i]-minVal)/(maxVal-minVal);
        }
        return scaledData;
    }

    public static int minIndex(double[] data){
        double minVal=data[0];
        int minValIndex=0;
        for(int i=1;i<data.length;i++){
            if(data[i]<minVal){
                minVal=data[i];
                minValIndex=i;
            }
        }
        return minValIndex;
    }

    public static int maxIndex(double[] data){
        double maxVal=data[0];
        int maxValIndex=0;
        for(int i=1;i<data.length;i++){
            if(data[i]>maxVal){
                maxVal=data[i];
                maxValIndex=i;
            }
        }
        return maxValIndex;
    }

    public static double trapezoidalIntegration(double[]wav,double[]refVal){
        double area=0;

        for(int i=0;i<wav.length-1;i++){
            area+=.5*(refVal[i]+refVal[i+1])*(wav[i+1]-wav[i]);
        }
        return area;
    }


    public static double[][] extractBandParameters(double[][] data,int nData,int nDim,double[] wavelength){

        double[][] region1000nm=new double[nData][35];
        double[][] region2000nm=new double[nData][38];
        double[] scaledWavelength;



        for(int i=0;i<nData;i++) {
            for (int j = 0; j < nDim; j++) {
                if (j < 35) {
                    region1000nm[i][j] = data[i][j];
                } else {
                    region2000nm[i][j - 35] = data[i][j];
                }
            }
        }



        //scaledWavelength=scale(wavelength);

        //Initialize band parameters
        double[] bandStrength1000nm=new double[nData];
        double[] bandStrength2000nm=new double[nData];
        double[] bandCenter1000nm=new double[nData];
        double[] bandCenter2000nm=new double[nData];
        int[] bandCenterIndex1000nm=new int[nData];
        int[] bandCenterIndex2000nm=new int[nData];
        double[] bandArea1000nm=new double[nData];
        double[] bandArea2000nm=new double[nData];
        double[] bandAreaRatio=new double[nData];
        double[] bandStrength=new double[nData];
        double[] bandCurvature=new double[nData];
        double[] bandTilt=new double[nData];

        double[] wavelengthRegion1000nm=Arrays.copyOfRange(wavelength, 0, 35);
        double[] wavelengthRegion2000nm=Arrays.copyOfRange(wavelength, 35, 73);

        double[][] featureVector=new double[nData][8];

        for(int i=0;i<nData;i++){
            //Calculate band center
            bandCenterIndex1000nm[i]=maxIndex(region1000nm[i]);
            bandCenterIndex2000nm[i]=maxIndex(region2000nm[i])+35;

            //Calculate band strength
            bandStrength1000nm[i]=region1000nm[i][bandCenterIndex1000nm[i]];
            bandStrength2000nm[i]=region2000nm[i][bandCenterIndex2000nm[i]-35];

            //Calculate band center index
            bandCenter1000nm[i]=wavelength[bandCenterIndex1000nm[i]];
            bandCenter2000nm[i]=wavelength[bandCenterIndex2000nm[i]];

            //Calculate band area
            bandArea1000nm[i]=trapezoidalIntegration(wavelengthRegion1000nm,region1000nm[i]);
            bandArea2000nm[i]=trapezoidalIntegration(wavelengthRegion2000nm,region2000nm[i]);

            //Calculate band area ratio
            bandAreaRatio[i]=bandArea2000nm[i]/bandArea1000nm[i];

            //Make feature vector from calculated band parameters
            featureVector[i][0]=bandStrength1000nm[i];
            featureVector[i][1]=bandStrength2000nm[i];
            //Olivine to pyroxene band ratio
            featureVector[i][2]=bandStrength2000nm[i]/bandStrength1000nm[i];
            featureVector[i][3]=bandCenter1000nm[i];
            featureVector[i][4]=bandCenter2000nm[i];
            featureVector[i][5]=bandArea1000nm[i];
            featureVector[i][6]=bandArea2000nm[i];
            featureVector[i][7]=bandAreaRatio[i];

            /*
            featureVector[i][8]=bandCenterIndex1000nm[i];
            featureVector[i][9]=bandCenterIndex2000nm[i];


            double CRdata750nm=(data[i][6]);
            double CRdata900nm=(data[i][13]+data[i][14])/2;
            double CRdata1000nm=(data[i][18]+data[i][19])/2;
            //double CRdata2000nm=(data[i][57]+data[i][58])/2;

            //Band strength
            featureVector[i][10]=CRdata1000nm/CRdata750nm;

            //Band curvature
            featureVector[i][11]=(CRdata750nm/CRdata900nm)+(CRdata1000nm/CRdata900nm);

            //Band tilt
            featureVector[i][12]=CRdata900nm/CRdata1000nm;

            //
            */

        }


        //IO.writeData(featureVector,nData,8,"J:\\LUNAR MINERAL MAPPING\\IISC\\BandParameterFeatureVector.txt");


        return featureVector;
    }
}
