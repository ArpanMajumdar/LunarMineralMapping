import java.io.File;
import java.util.Arrays;

/**
 * Created by ARPAN on 10-12-2016.
 */
public class Test {
    public static void main(String args[]){
           File ref=new File("./data/LunarData/Library Spectra/HighCaPyroxene.txt");
           File wav=new File("./data/LunarData/Library Spectra/HighCaPyroxeneWavelength.txt");
           File M3wav=new File("./data/LunarData/M3wavelength_1x75.txt");
           File data=new File("./data/LunarData/CRSdata_450000x73_MATLAB.txt");

           int nData=450000,nDim=73,imgDim1=1500,imgDim2=300;
           double[] libSpectra=IO.readDoubleVec(ref);
           double[] libWavelength=IO.readDoubleVec(wav);
           double[] M3wavelength=IO.readDoubleVec(M3wav,75);
           double[][] CRSdata=IO.readDoubleMat(data,nData,nDim);
           double[] M3wavelength_73= Arrays.copyOfRange(M3wavelength,0,73);

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


       System.out.println("Library reflectance");
       for(int i=0;i<libSpectra.length;i++){
           System.out.print(libSpectra[i]+"\t");
       }
       System.out.println();

        System.out.println("Library wavelength");
        for(int i=0;i<libWavelength.length;i++){
            System.out.print(libWavelength[i]+"\t");
        }
        System.out.println();

        System.out.println("M3 wavelength");
        for (double i:M3wavelength) {
            System.out.print(i+",");
        }
        System.out.println();


        System.out.println("Interpolated reflectance");
            for (double i:interpolatedLibSpectra) {
                System.out.print(i+",");
            }
        System.out.println();

        System.out.println("Continuum Removed Interpolated reflectance");
        for (double i:CRinterpolatedLibSpectra) {
            System.out.print(i+",");
        }
        System.out.println();

    }
}
