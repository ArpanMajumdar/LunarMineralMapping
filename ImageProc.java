

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Jay
 */
public class ImageProc {
    public static void imwrite(File file,int[][] classificationMat,int imgDim1,int imgDim2){
         BufferedImage image = new BufferedImage(imgDim2, imgDim1, BufferedImage.TYPE_INT_RGB);

        for (int i = 0; i < imgDim1; i++) {
            for (int j = 0; j < imgDim2; j++) {
                
                switch(classificationMat[i][j]){
                        case 1: image.setRGB(j, i, Color.BLUE.getRGB());
                        break;
                        case 2: image.setRGB(j, i, Color.CYAN.getRGB());
                        break;
                        case 3: image.setRGB(j, i, Color.GREEN.getRGB());
                        break;
                        case 4: image.setRGB(j, i, Color.ORANGE.getRGB());
                        break;
                        case 5: image.setRGB(j, i, Color.RED.getRGB());
                        break;
                }      
            }
        }
        try{
        ImageIO.write(image, "png", file);
        }
        catch(IOException e){
          e.printStackTrace();
        }
    }
    
    public static void imagesc(File file,int[][]classificationMat,int nClass,int imgDim1,int imgDim2){
        BufferedImage image = new BufferedImage(imgDim2, imgDim1, BufferedImage.TYPE_INT_RGB);
        int index,rgb;
        float[][] jet=colormapJet(nClass);
        for (int i = 0; i < imgDim1; i++) {
            for (int j = 0; j < imgDim2; j++) {
                index=classificationMat[i][j];
                if(index==-1){
                    image.setRGB(j, i, Color.BLACK.getRGB()); 
                }
                else{
                    rgb=new Color(jet[index][0],jet[index][1],jet[index][2]).getRGB();
                    image.setRGB(j, i, rgb);
                }
            }
        }
        try{
        ImageIO.write(image, "png", file);
        }
        catch(IOException e){
          e.printStackTrace();
        }
    }
    
    public static float[][] colormapJet(int nClass){
        int n;
        n=(int)Math.ceil((float)nClass/4);
        ArrayRealVector u,r,g,b;
        RealVector R,G,B;
        
        u=new ArrayRealVector(3*n-1);
        for(int i=0;i<n;i++){
            u.setEntry(i, (i+1.0)/n);
            u.setEntry(u.getDimension()-i-1, (i+1.0)/n);
        }
        u.setSubVector(n,new ArrayRealVector(n,1));
              
        g=new ArrayRealVector(u.getDimension());
        
        float m;
        m=(float)Math.ceil((float)n/2);
        if(nClass%4==1)
            m=m-1;
        
        for(int i=0;i<g.getDimension();i++){
            g.setEntry(i, (i+1+m));
        }
        
        r=g.add(new ArrayRealVector(g.getDimension(),n));
        b=g.subtract(new ArrayRealVector(g.getDimension(),n));
        
        
        R=new ArrayRealVector();
        G=new ArrayRealVector();
        B=new ArrayRealVector();
        
        for(int i=0;i<r.getDimension();i++){
            if(r.getEntry(i)<=nClass)
                R=R.append(r.getEntry(i));
        }
        
        for(int i=0;i<g.getDimension();i++){
            if(g.getEntry(i)<=nClass)
                G=G.append(g.getEntry(i));
        }
        
        for(int i=0;i<b.getDimension();i++){
            if(b.getEntry(i)>=1)
                B=B.append(b.getEntry(i));
        }
                
        float[][] J=new float[nClass][3];
        int index;
        for(int i=0;i<R.getDimension();i++){
            index=(int)R.getEntry(i);
            J[index-1][0]=(float)u.getEntry(i);
        }
        
        for(int i=0;i<G.getDimension();i++){
            index=(int)G.getEntry(i);            
            J[index-1][1]=(float)u.getEntry(i);
        }
        
        for(int i=u.getDimension()-B.getDimension(),j=0;i<u.getDimension();i++,j++){
            index=(int)B.getEntry(j);
            J[index-1][2]=(float)u.getEntry(i);
        }
        
        /*
        System.out.println("\nRGB Matrix:");
        for(int i=0;i<nClass;i++){
            for(int j=0;j<3;j++){
                System.out.print(J[i][j]+"\t");
            }
            System.out.println();
        }
        */
        return J;
    }
}
