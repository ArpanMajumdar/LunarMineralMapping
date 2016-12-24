/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hyperspectral.workflow;

import static hyperspectral.workflow.ImageProc.imagesc;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

/**
 *
 * @author Arpan
 */
public class SpatialInformationExtraction {
    
    public static double[][][] readSpectralCube(File file,int imgDim1,int imgDim2,int bands){
        Scanner sc;
        double[][][] dataCube=new double[1500][300][73];        
        
        try{
            sc=new Scanner(file);
            sc.useDelimiter(",|\\n");
            
            for(int j=0;j<imgDim2;j++){
                for(int i=0;i<imgDim1;i++){
                    for(int k=0;k<bands;k++){
                        dataCube[i][j][k]=Double.parseDouble(sc.next());
                    }
                }
            }            
        }
        catch(FileNotFoundException e){
            e.printStackTrace();
        }
        
        return dataCube;
    }
    
    public static double[][] SAMmeanWindowFilter(double[][][] dataCube,int imgDim1,int imgDim2,int bands){
        ArrayList<ArrayRealVector> arrVecList;
        ArrayRealVector refVec=null;
        ArrayRealVector neighbourVec;
        double[][] arrVec=new double[9][bands];
        double[] SAMangle=new double[8];
        double[][] SAMsimilarityImg=new double[imgDim1][imgDim2];
        
        for(int i=1;i<imgDim1-1;i++){
            for(int j=1;j<imgDim2-1;j++){
                for(int k=1;k<bands;k++){
                    arrVec[0][k]=dataCube[i-1][j-1][k];
                    arrVec[1][k]=dataCube[i-1][j][k];
                    arrVec[2][k]=dataCube[i-1][j+1][k];
                    arrVec[3][k]=dataCube[i][j-1][k];
                    arrVec[4][k]=dataCube[i][j][k];
                    arrVec[5][k]=dataCube[i][j+1][k];
                    arrVec[6][k]=dataCube[i+1][j-1][k];
                    arrVec[7][k]=dataCube[i+1][j][k];
                    arrVec[8][k]=dataCube[i+1][j+1][k];
                }
                
                arrVecList=new ArrayList<>();
                
                for(int k=0;k<9;k++){
                    if(k==4)
                        refVec=new ArrayRealVector(arrVec[k]);
                    else
                        arrVecList.add(new ArrayRealVector(arrVec[k]));
                }                
               
                for(int k=0;k<arrVecList.size();k++){
                    neighbourVec=arrVecList.get(k);
                    SAMangle[k]=Math.toDegrees(Math.acos(refVec.dotProduct(neighbourVec)/(refVec.getNorm()*neighbourVec.getNorm())));
                }
                
                SAMsimilarityImg[i][j]=StatUtils.mean(SAMangle);
            }
        }
        return SAMsimilarityImg;
    }
    
    public static double[][] meanFilter(double[][] SAMsimilarityImg,int imgDim1,int imgDim2,int bands){
        
        double[][] avgSAMsimilarityImg=new double[imgDim1][imgDim2];
        double[] window=new double[9];
        for(int i=1;i<imgDim1-1;i++){
            for(int j=1;j<imgDim2-1;j++){
                window[0]=SAMsimilarityImg[i-1][j-1];
                window[1]=SAMsimilarityImg[i-1][j];
                window[2]=SAMsimilarityImg[i-1][j+1];
                window[3]=SAMsimilarityImg[i][j-1];
                window[4]=SAMsimilarityImg[i][j];
                window[5]=SAMsimilarityImg[i][j+1];
                window[6]=SAMsimilarityImg[i+1][j-1];
                window[7]=SAMsimilarityImg[i+1][j];
                window[8]=SAMsimilarityImg[i+1][j+1];
                
                int nnz=0;
                for(double k:window){
                    if(k!=0.0)
                        nnz++;
                }
                avgSAMsimilarityImg[i][j]=StatUtils.sum(window)/nnz;
            }            
        }
        return avgSAMsimilarityImg;
    }
    
    public static void main(String args[]){
        int imgDim1=1500,
            imgDim2=300,
            bands=73;
        
        String root="J:\\Java programs\\Hyperspectral Workflow\\Data files\\";
        File file=new File(root+"Data Cube.txt");
        double[][][] dataCube=readSpectralCube(file,imgDim1,imgDim2,bands);
        double[][] SAMsimilarityImg=SAMmeanWindowFilter(dataCube,imgDim1,imgDim2,bands);
        double[][] avgSAMsimilarityImg1=meanFilter(SAMsimilarityImg,imgDim1,imgDim2,bands);
        double[][] avgSAMsimilarityImg2=meanFilter(avgSAMsimilarityImg1,imgDim1,imgDim2,bands);
        
        
        for(int i=0;i<imgDim1;i++){
            System.out.print(i+" : ");
            for(int j=0;j<imgDim2;j++){
                System.out.print(avgSAMsimilarityImg2[i][j]+"\t");
            }
            System.out.println();
        }
        
        int[][] thresholdedImg=new int[imgDim1][imgDim2];
        
        for(int i=0;i<imgDim1;i++){
            for(int j=0;j<imgDim2;j++){
                if(avgSAMsimilarityImg2[i][j]>0.0 && avgSAMsimilarityImg2[i][j]<9.0)
                    thresholdedImg[i][j]=1;
                else
                    thresholdedImg[i][j]=0;
            }
        }    
            
        File binaryConnCompImg=new File(root+"BinaryConnCompImg.png");
        imagesc(binaryConnCompImg,thresholdedImg,2,imgDim1,imgDim2);
        
    }    
}