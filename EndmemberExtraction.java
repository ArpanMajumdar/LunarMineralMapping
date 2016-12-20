/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


import org.apache.commons.math3.linear.ArrayRealVector;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

/**
 *
 * @author Jay
 */
public class EndmemberExtraction {

    public static ArrayList<Integer> ORASIS(double[][] data,int nData,int nDim,double threshold,int[] exemplarLabel){
        ArrayRealVector vec;
        ArrayList<ArrayRealVector> X=new ArrayList<>();
        ArrayList<ArrayRealVector> E=new ArrayList<>();
        ArrayList<Integer> exemplarIndex=new ArrayList<>();

        for(int i=0;i<nData;i++){
            vec=new ArrayRealVector(data[i]);
            vec.unitize();
            X.add(vec);
        }

        E.add(X.get(0));
        exemplarIndex.add(0);
        double t=Math.sqrt(2*(1-threshold));

        //Add first element of test spectra to set of exemplar spectra
        exemplarLabel[0]=0;

        boolean flag;
        double maxCos,sigmaMin,sigmaMax,dotXR,dotER,cosTheta;

        double[] vecR=new double[nDim];
        for(int i=0;i<nDim;i++){
            vecR[i]=1/Math.sqrt(nDim);
        }
        ArrayRealVector R=new ArrayRealVector(vecR);

        ArrayRealVector exemplarSpec,testSpec;

        for(int i=0;i<X.size();i++){
            if(i==0 || exemplarLabel[i]==-1){
                continue;
            }



            flag=false;
            maxCos=0;
            testSpec=X.get(i);
            dotXR=testSpec.dotProduct(R);
            sigmaMin=dotXR-t;
            sigmaMax=dotXR+t;

            for(int j=0;j<E.size();j++){
                exemplarSpec=E.get(j);
                dotER=exemplarSpec.dotProduct(R);

                if(dotER<sigmaMax && dotER>sigmaMin){
                    cosTheta=testSpec.dotProduct(exemplarSpec);

                    if(cosTheta>threshold){
                        //Test spectra is similar to one of the exemplar spectra
                        if(cosTheta>maxCos){
                            maxCos=cosTheta;
                            exemplarLabel[i]=j;
                            //System.out.println("Count: "+i+"\texemplarLabel: "+exemplarLabel[i]);
                            flag=true;
                        }
                    }
                }
            }

            if(!flag){
                //Test spectra is unique, add it to set of exemplars
                E.add(testSpec);
                exemplarIndex.add(i);
                exemplarLabel[i]=E.size()-1;
                //System.out.println("Count: "+i+"\texemplarLabel: "+exemplarLabel[i]);
            }

        }
        return exemplarIndex;
    }

    public static void exemplarFrequency(int[] exemplarLabel,double curThresholdAbundance){
        int[] exemplarFreq=new int[100];
        HashSet<Integer> uniqueExemplars=new HashSet<>();
        int totalExemplars=0;
        for(int i=0;i<exemplarLabel.length;i++){
            if(exemplarLabel[i]!=-1){
                exemplarFreq[exemplarLabel[i]]++;
                totalExemplars++;
                uniqueExemplars.add(exemplarLabel[i]);
            }
        }
        
        int nExemplar=uniqueExemplars.size();
        double[] fracAbundance=new double[nExemplar];
        boolean[] isRejected=new boolean[nExemplar];
        for(boolean i:isRejected){
            i=false;
        }
        //System.out.println("nExemplar="+nExemplar);
            
        //System.out.println("Fractional abundance:");
        for(int i=0;i<nExemplar;i++){
            fracAbundance[i]=(double)exemplarFreq[i]/totalExemplars;
            
            if(fracAbundance[i]<curThresholdAbundance){
                isRejected[i]=true;
            }
            //System.out.println(i+" : "+(fracAbundance[i]*100));
        }
        //System.out.println();
        int count=0;
        for(int i=0;i<exemplarLabel.length;i++){
            if(exemplarLabel[i]!=-1){
                if(isRejected[exemplarLabel[i]]){
                    exemplarLabel[i]=-1;
                }                
            }
        }
    }   

    public static void iterativeORASIS(double[][] data,int nData,int nDim,double minThresholdAngle,double maxThresholdAngle,double stepsize,double minThresholdAbundance,int[] exemplarLabel){
        int itr=1;
        double curThresholdAngle=minThresholdAngle,curThresholdAbundance=minThresholdAbundance;
        ArrayList<Integer> exemplarIndex;
        
        
        while(curThresholdAngle<=maxThresholdAngle){
                exemplarIndex=ORASIS(data,nData,nDim,curThresholdAngle,exemplarLabel);
                
                //System.out.println("Itr="+itr);
                exemplarFrequency(exemplarLabel,curThresholdAbundance);
                
                curThresholdAngle+=stepsize;
                //curThresholdAbundance+=.005;
                itr++;
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
    
    public static int[] assignInitialLabels(double[][] data,int nData,int nDim,int imgDim1,int imgDim2,double minThresholdAngle,double maxThresholdAngle,double stepSize,double minThresholdAbundance,String filepath){

        int[] exemplarLabel=new int[nData];
        int[][] exemplarMat;

        
        iterativeORASIS(data,nData,nDim,minThresholdAngle,maxThresholdAngle,stepSize,minThresholdAbundance,exemplarLabel);
        HashSet<Integer> uniqueExemplars=new HashSet<>();
        for(int i=0;i<nData;i++){
            uniqueExemplars.add(exemplarLabel[i]);
        }
        
        
        //System.out.println("Unique exemplars:");
        HashMap<Integer,Integer> exemplars=new HashMap<>();
        int count=0;
        for(int i:uniqueExemplars){
            if(i!=-1){
                
                exemplars.put(i, count);
                //System.out.print(i+"\t");
                count++;
            }            
        }
        System.out.println("Exemplar Count:"+count);
        //System.out.println();
        
        
        for(int i=0;i<nData;i++){
            if(exemplarLabel[i]!=-1){
                exemplarLabel[i]=exemplars.get(exemplarLabel[i]);
            }
        }
        exemplarMat=reshape(exemplarLabel,imgDim1,imgDim2);
        IO.writeData(exemplarMat,imgDim1,imgDim2,filepath+"exemplarMat.txt");
        
        /*
        for(int i=0;i<imgDim1;i++){
            System.out.print(i+" : ");
            for(int j=0;j<imgDim2;j++){
                System.out.print(exemplarMat[i][j]+"\t");
            }
            System.out.println();
        }
        */
        
        //Display exemplar classification image
        File exemplarImg=new File(filepath+"ExemplarImg.png");
        ImageProc.imagesc(exemplarImg, exemplarMat,uniqueExemplars.size() , imgDim1, imgDim2);

        return exemplarLabel;
        
    }
}
