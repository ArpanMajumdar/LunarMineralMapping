import java.io.File;

/**
 * Created by ARPAN on 10-12-2016.
 */
public class Test {
    public static void main(String args[]){
       File ref=new File("./data/LunarData/Library Spectra/HighCaPyroxene.txt");
       File wav=new File("./data/LunarData/Library Spectra/HighCaPyroxeneWavelength.txt");
       File M3wav=new File("./data/LunarData/M3wavelength_1x75.txt");

       double[] reflectance=IO.readDoubleVec(ref);
       double[] wavelength=IO.readDoubleVec(wav);

       /*
       for(int i=0;i<reflectance.length;i++){
           System.out.print(reflectance[i]+"\t");
       }
       System.out.println();

        for(int i=0;i<wavelength.length;i++){
            System.out.print(wavelength[i]+"\t");
        }
        System.out.println();
        */



    }
}
