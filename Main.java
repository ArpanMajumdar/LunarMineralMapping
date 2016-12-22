import java.io.*;
import java.util.Arrays;

/**
 * Created by ARPAN on 12-11-2016.
 */
public class Main {

    public static void main(String args[]){

        // Initialize variables
        int imgDim1=1500,
                imgDim2=300,
                nData=imgDim1*imgDim2,
                nDim=75;

        double[][] data;double[] M3wavelength;

        System.out.println("Reading data...");
        File file_in=new File("./data/LunarData/RawData_450000x75.txt");
        File wavFile=new File("./data/LunarData/M3wavelength_1x75.txt");

        //Read data from binary file
        //data=IO.readBinaryData(file_in,nData,nDim);

        //Read data from text file
        data=IO.readDoubleMat(file_in,nData,nDim);
        M3wavelength=IO.readDoubleVec(wavFile,nDim);
        //double[][] CRdata=data;
        System.out.println("Reading done!");

        //Perform continuum removal of data
        System.out.println("Performing continuum removal...");
        double[][] CRdata=new double[nData][nDim];

        for(int i=0;i<nData;i++){
            CRdata[i]=HyperspectralToolbox.performContinuumRemoval(data[i],M3wavelength,nDim);
        }

        IO.writeData(CRdata,nData,nDim,"./data/LunarData/CRdata_450000x75.txt");
        System.out.println("Continuum removal done!");


        //Perform de-noising of spectra by running average filter
        System.out.println("Performing de-noising of spectra by running average filter...");
        int pt=3;
        double[][] CRSdata=HyperspectralToolbox.runningAvgFilter(CRdata,nData,nDim,pt);
        nDim=nDim-pt+1;
        IO.writeData(CRSdata,nData,nDim,"./data/LunarData/CRSdata_450000x73.txt");
        System.out.println("De-noising done!");

        /*
        //Extract band parameters
        double[][] BandParameterFeatureVector=BandParameterDetermination.extractBandParameters(CRSdata,nData,nDim,M3wavelength);
        */

        //Perform initial labelling
        System.out.println("Obtaining initial labelling by Iterative ORASIS algorithm...");
        double minThresholdAngle=.8,
                maxThresholdAngle=.9,
                stepSize=.01,
                minThresholdAbundance=.01;
        String imgPath="./data/LunarData/";
        int[] init=EndmemberExtraction.assignInitialLabels(CRSdata,nData,nDim,imgDim1,imgDim2,minThresholdAngle,maxThresholdAngle,stepSize,minThresholdAbundance,imgPath);
        IO.writeData(init,nData,"./data/LunarData/InitialLabel.txt");
        System.out.println("Initial labels obtained!");

        //Randomly initializing labels for rejected pixels
        System.out.println("Initializing the rejected pixels with random labels");
        int[] EMinit=new int[nData];
        for(int i=0;i<nData;i++){
            if(init[i]==-1){
                EMinit[i]=(int)(Math.random()*3);
            }
        }

        System.out.println("Running Expectation-Maximization algorithm...");
        int[] finalLabel=ExpectationMaximization.EM(CRSdata,EMinit,nData,nDim,3);
        IO.display(finalLabel,"Final Label",nData);
        //Display exemplar classification image
        File EMOutputImg=new File("./data/LunarData/EMOutputImg.png");
        ImageProc.imagesc(EMOutputImg, reshape(finalLabel,imgDim1,imgDim2),3 , imgDim1, imgDim2);
        System.out.println("EM agorithm finished!");



        //Read Library spectra
        File ref=new File("./data/LunarData/Library Spectra/HighCaPyroxene.txt");
        File libwav=new File("./data/LunarData/Library Spectra/HighCaPyroxeneWavelength.txt");

        double[] libSpectra=IO.readDoubleVec(ref);
        double[] libWavelength=IO.readDoubleVec(libwav);
        double[] M3wavelength_73= Arrays.copyOfRange(M3wavelength,0,nDim);

        double[] interpolatedLibSpectra=HyperspectralToolbox.interpLinear(libWavelength,libSpectra,M3wavelength_73);

        //Perform continuum removal of library spectra
        double[] CRinterpolatedLibSpectra=HyperspectralToolbox.performContinuumRemoval(interpolatedLibSpectra,M3wavelength_73,nDim);

        //Invert library spectra for matching
        CRinterpolatedLibSpectra=HyperspectralToolbox.invert(CRinterpolatedLibSpectra);

        //Multiplying data by 10 to match order and performing sanity check
        for(int i=0;i<nData;i++){
            for(int j=0;j<nDim;j++){
                CRSdata[i][j]=CRSdata[i][j]*10;
                if(CRSdata[i][j]>1){
                    CRSdata[i][j]=1;
                }
                else if(CRSdata[i][j]<0){
                    CRSdata[i][j]=0;
                }
            }
        }

        //Calculate SAM angle for each pixel and displaying as SAM intensity image
        double[][] SAMmat=HyperspectralToolbox.SAMIntensityImage(CRSdata,CRinterpolatedLibSpectra,nData,nDim,imgDim1,imgDim2);
        IO.writeData(SAMmat,imgDim1,imgDim2,"./data/LunarData/SAM Intensity Image.txt");
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


