import org.apache.commons.math3.stat.StatUtils;

import java.io.File;
import java.util.Arrays;

/**
 * Created by ARPAN on 10-12-2016.
 */
public class MainSalinasDataset {
    public static void main(String args[]){
        // Initialize variables
        int imgDim1=512,
                imgDim2=217,
                nData=imgDim1*imgDim2,
                nDim=204;

        double[][] data;double[] SalinasWavelength;

        System.out.println("Reading data...");
        File dataFile=new File("./data/Salinas/SalinasRawData_111104x204.txt");
        File wavFile=new File("./data/Salinas/SalinasWavelength_1x204.txt");

        //Read data from text file
        data=IO.readDoubleMat(dataFile,nData,nDim);
        SalinasWavelength=IO.readDoubleVec(wavFile,nDim);

        //double[][] CRdata=data;
        System.out.println("Reading done!");

        /*
        //Perform continuum removal of data
        System.out.println("Performing continuum removal...");
        double[][] CRdata=new double[nData][nDim];

        for(int i=0;i<nData;i++){
            CRdata[i]=HyperspectralToolbox.performContinuumRemoval(data[i],SalinasWavelength,nDim);
        }

        IO.writeData(CRdata,nData,nDim,"./data/Salinas/SalinasCRData_111104x204.txt");
        System.out.println("Continuum removal done!");
        */

        /*
        //Perform de-noising of spectra by running average filter
        System.out.println("Performing de-noising of spectra by running average filter...");
        int pt=3;
        double[][] CRSdata=HyperspectralToolbox.runningAvgFilter(CRdata,nData,nDim,pt);
        nDim=nDim-pt+1;
        IO.writeData(CRSdata,nData,nDim,"./data/LunarData/CRSdata_450000x73.txt");
        System.out.println("De-noising done!");


        //Extract band parameters
        double[][] BandParameterFeatureVector=BandParameterDetermination.extractBandParameters(CRSdata,nData,nDim,M3wavelength);
        */

        //Perform initial labelling
        System.out.println("Obtaining initial labelling by Iterative ORASIS algorithm...");
        double minThresholdAngle=.99,
                maxThresholdAngle=.996,
                stepSize=.001,
                minThresholdAbundance=0.005;
        String imgPath="./data/Salinas/";
        int[] init=EndmemberExtraction.assignInitialLabels(data,nData,nDim,imgDim1,imgDim2,minThresholdAngle,maxThresholdAngle,stepSize,minThresholdAbundance,imgPath);
        IO.writeData(init,nData,"./data/Salinas/IndianPinesORASISLabel.txt");
        System.out.println("Initial labels obtained!");

        int nClass= 15;
        //Randomly initializing labels for rejected pixels
        System.out.println("Initializing the rejected pixels with random labels");
        int[] EMinit=new int[nData];
        for(int i=0;i<nData;i++){
            if(init[i]==-1){
                EMinit[i]=(int)(Math.random()*nClass);
            }
            else{
                EMinit[i]=init[i];
            }
        }

        //Running Expectation-Maximization algorithm
        System.out.println("Running Expectation-Maximization algorithm...");
        int[] finalLabel=ExpectationMaximization.EM(data,EMinit,nData,nDim,nClass);
        IO.display(finalLabel,"EMLabel",nData);
        //Display exemplar classification image
        File EMOutputImg=new File("./data/Salinas/EMOutputImg.png");
        ImageProc.imagesc(EMOutputImg, reshape(finalLabel,imgDim1,imgDim2),nClass , imgDim1, imgDim2);
        System.out.println("EM agorithm finished!");
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
}
