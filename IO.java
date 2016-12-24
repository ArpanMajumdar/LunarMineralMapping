


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.io.*;
import java.util.Scanner;
/**
 *
 * @author ARPAN
 */
public class IO {

    public static double[][] readDoubleMat(File file,int row,int col){
        double[][] data=new double[row][col];
        try{
            Scanner sc=new Scanner(file);
            sc.useDelimiter(",|\\n");

            for(int i=0;i<row;i++)
                for(int j=0;j<col;j++){
                    data[i][j]=Double.parseDouble(sc.next());
                }

        }
        catch(FileNotFoundException e){
            e.printStackTrace();
        }
        return data;
    }

    public static int[] readIntVec(File file,int nData){
        int[] data=new int[nData];
        try {
            Scanner sc = new Scanner(file);
            sc.useDelimiter(",|\\n");

            for (int i = 0; i < nData; i++) {
                data[i] = sc.nextInt();
            }
        }
        catch(FileNotFoundException e){
            e.printStackTrace();
        }
        return data;
    }

    public static double[] readDoubleVec(File file){
        double[] data;
        try{
            int nData=countLines(file.getAbsolutePath());
            data=new double[nData];

            Scanner sc=new Scanner(file);
            sc.useDelimiter(",|\\n");
            for(int i=0;i<nData;i++){
                data[i]=sc.nextDouble();
            }
            return data;
        }
        catch(IOException e){
            e.printStackTrace();
            data=new double[0];
            return data;
        }
    }

    public static double[][] readBinaryData(File file,int nData,int nDim){
        FileInputStream inStream;
        double[][] data=new double[nData][nDim];
        try{
            inStream= new FileInputStream (file);
            DataInputStream input = new DataInputStream (inStream);

            for(int i=0;i<nData;i++)
                for(int j=0;j<nDim;j++){
                    data[i][j]=input.readDouble();
                }
        }
        catch (FileNotFoundException e){
            e.printStackTrace();
        }
        catch (IOException e){
            e.printStackTrace();
        }

        return data;
    }

    public static double[] readDoubleVec(File file,int nDim){
		Scanner sc;
		double[] wavelength=new double[nDim];
		try{
			sc=new Scanner(file);
            sc.useDelimiter(",|\\n");
			for(int i=0;i<nDim;i++){
				wavelength[i]=Double.parseDouble(sc.next());
				//System.out.print(wavelength[i]+"\t");
			}
		}
		catch(FileNotFoundException e){
			System.out.println("Wavelength file not found\n"+e);
			System.exit(0);
		}
		return wavelength;
	}

    public static double[][][] readDouble3DMat(File file,int imgDim1,int imgDim2,int bands){
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

	public static void writeData(double[][] data,int nData,int nDim,String filepath){
        File file=new File(filepath);
        PrintWriter writer;

        try{

            writer=new PrintWriter(file);
            for(int i=0;i<nData;i++){
                for(int j=0;j<nDim;j++){
                    if(j==nDim-1){
                        writer.print(data[i][j]);
                    }
                    else{
                        writer.print(data[i][j]+",");
                    }

                }
                writer.println();
            }
            writer.close();
        }
        catch(FileNotFoundException e){
            e.printStackTrace();
        }
    }

    public static void writeData(int[][] data,int nData,int nDim,String filepath){
        File file=new File(filepath);
        PrintWriter writer;

        try{

            writer=new PrintWriter(file);
            for(int i=0;i<nData;i++){
                for(int j=0;j<nDim;j++){
                    if(j==nDim-1){
                        writer.print(data[i][j]);
                    }
                    else{
                        writer.print(data[i][j]+",");
                    }
                }
                writer.println();
            }
            writer.close();
        }
        catch(FileNotFoundException e){
            e.printStackTrace();
        }
    }

    public static void writeData(int[] data,int nData,String filepath){
        File file=new File(filepath);
        PrintWriter writer;

        try{

            writer=new PrintWriter(file);
            for(int i=0;i<nData;i++){
                if(i==nData-1){
                    writer.print(data[i]);
                }
                else{
                    writer.print(data[i]+",");
                }
            }
            writer.close();
        }
        catch(FileNotFoundException e){
            e.printStackTrace();
        }
    }


    public static void display(double[][] mat,String name,int row,int col){
        System.out.println(name+":");
        for(int i=0;i<row;i++){
            for(int j=0;j<col;j++){
                System.out.print(mat[i][j]+",");
            }
            System.out.println();
            System.out.println();
        }
    }

    public static void display(int[] vec,String name,int len){
        System.out.println(name+":");
        for(int i=0;i<len;i++){
            System.out.print(vec[i]+",");
        }
        System.out.println();
        System.out.println();
    }

    public static void display(RealMatrix mat,String name){

        int row=mat.getRowDimension();
        int col=mat.getColumnDimension();

        System.out.println(name+":"+row+"x"+col);
        for(int i=0;i<row;i++){
            for(int j=0;j<col;j++){
                System.out.print(mat.getEntry(i,j)+",");
            }
            System.out.println();
        }
    }

    public static void display(RealVector vec, String name){
        int len=vec.getDimension();
        System.out.println(name+":"+len);
        for(int i=0;i<len;i++){
            System.out.print(vec.getEntry(i)+",");
        }
        System.out.println();
    }

    public static int countLines(String filename) throws IOException {
        InputStream is = new BufferedInputStream(new FileInputStream(filename));
        try {
            byte[] c = new byte[1024];
            int count = 0;
            int readChars = 0;
            boolean empty = true;
            while ((readChars = is.read(c)) != -1) {
                empty = false;
                for (int i = 0; i < readChars; ++i) {
                    if (c[i] == '\n') {
                        ++count;
                    }
                }
            }
            return (count == 0 && !empty) ? 1 : count;
        } finally {
            is.close();
        }
    }

}


