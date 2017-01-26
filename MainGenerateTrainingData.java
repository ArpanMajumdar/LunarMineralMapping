import java.io.File;
import java.util.Arrays;

/**
 * Created by ARPAN on 10-12-2016.
 */
public class MainGenerateTrainingData {
    public static void main(String args[]){


          File dataFile=new File("J:\\LUNAR MINERAL MAPPING\\IISC\\IrisDataset.txt");
          File labelFile=new File("J:\\LUNAR MINERAL MAPPING\\IISC\\IrisDatasetLabel.txt");

          int nData=150,nDim=4,nClass=3;
          double[][] data=IO.readDoubleMat(dataFile,nData,nDim);
          int[] label=IO.readIntVec(labelFile,nData);
          HyperspectralToolbox.labelCount(label,3);

          //This means that 20 members from Class 0,20 from Class 1 and 20 from Class 2.
          int[] membercount={20,20,20};
          int[] trainingLabel=new int[60];//20+20+20
          double[][] trainingData=HyperspectralToolbox.generateTrainingData(data,nData,nDim,label,nClass,membercount,trainingLabel);

        System.out.println("Training data:");
          for(int i=0;i<trainingData.length;i++){
              for(int j=0;j<4;j++){
                  System.out.print(trainingData[i][j]+",");
              }
              System.out.print("\t"+trainingLabel[i]);
              System.out.println();
          }
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
