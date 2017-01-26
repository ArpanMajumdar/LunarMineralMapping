import java.io.File;

/**
 * Created by ARPAN on 24-01-2017.
 */
public class MainBandParametersFeatureVector {
    public static  void main(String args[]){
        // Initialize variables
        int imgDim1=1500,
                imgDim2=300,
                nData=imgDim1*imgDim2,
                nDim=73;

        double[][] data;double[] M3wavelength;

        System.out.println("Reading data...");
        File dataFile=new File("./data/LunarData/CRSdata_450000x73_MATLAB.txt");
        File wavFile=new File("./data/LunarData/M3wavelength_1x75.txt");

        data=IO.readDoubleMat(dataFile,nData,nDim);
        M3wavelength=IO.readDoubleVec(wavFile,nDim);


        //Extract band parameters
        double[][] BandParameterFeatureVector=BandParameterDetermination.extractBandParameters(data,nData,nDim,M3wavelength);
        BandParameterFeatureVector=HyperspectralToolbox.normalize(BandParameterFeatureVector,nData,8);
        IO.writeData(BandParameterFeatureVector,nData,BandParameterFeatureVector[0].length,"J:\\LUNAR MINERAL MAPPING\\IISC\\BandParameterFeatureVector.txt");

    }
}
