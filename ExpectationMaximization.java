import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.stat.StatUtils;

import java.io.File;

import static javafx.scene.input.KeyCode.R;

/**
 * Created by ARPAN on 17-11-2016.
 */
class Model{
    RealMatrix mu;
    RealMatrix[] sigma;
    RealVector w;
}

public class ExpectationMaximization {
    static final double pi=3.14159265;


    public static int[] EM(double[][] data,int[] label,int nData,int nDim,int init){
        double tol=1E-6;
        int maxitr=10;
        double[] llh=new double[maxitr];

        for(double i:llh)
            i=Double.NEGATIVE_INFINITY;

        RealMatrix X=new Array2DRowRealMatrix(data);
        X=X.transpose();
        //IO.display(X,"X");

        RealMatrix R=initialization(nData,label,init);
        //IO.display(R,"R");


        for(int i=1;i<maxitr;i++){
            System.out.print("itr="+i+"\t"+":");

            label=getLabels(R,nData,init);
            IO.display(label,"Label",nData);

            Model model=maximization(X,R,nData,nDim,init);
            R=expectation(X,model,nData,nDim,init,llh,i);

            if(Math.abs(llh[i]-llh[i-1])<tol*Math.abs(llh[i])){
                break;
            }

        }
        return label;
    }

    public static RealMatrix initialization(int nData,int init){
        int[] label=new int[nData];

        for(int i:label)
            i=(int)Math.ceil(init*Math.random());

        double[][] R=new double[nData][init];

        for(int i=0;i<nData;i++){
            R[i][label[i]]=1;
        }
        return new Array2DRowRealMatrix(R);
    }

    public static RealMatrix initialization(int nData,int[] label,int init){
        RealMatrix R=new Array2DRowRealMatrix(nData,init);

        for(int i=0;i<nData;i++){
            R.setEntry(i,label[i],1);
        }
        return R;
    }

    public static int[] getLabels(RealMatrix R,int nData,int k){

        double max;int maxindex;
        int[] label=new int[nData];
        for(int i=0;i<nData;i++){
            max=R.getEntry(i,0);
            maxindex=0;
            for(int j=1;j<k;j++){
                if(R.getEntry(i,j)>max){
                    max=R.getEntry(i,j);
                    maxindex=j;
                }
            }
            label[i]=maxindex;
        }
        return label;
    }

    public static Model maximization(RealMatrix X,RealMatrix R,int nData,int nDim,int k){
        RealVector nk=new ArrayRealVector(k);

        for(int i=0;i<k;i++){
            nk.setEntry(i,StatUtils.sum(R.getColumn(i)));
        }


        RealVector w=new ArrayRealVector(k);
        for(int i=0;i<k;i++){
            w.setEntry(i,nk.getEntry(i)/nData);
        }

        RealMatrix XtimesR=R.preMultiply(X);

        //mu: nDim x k
        RealMatrix mu=new Array2DRowRealMatrix(nDim,k);
        for(int i=0;i<nDim;i++){
            RealVector row=new ArrayRealVector(XtimesR.getRowVector(i));
            mu.setRowVector(i,row.ebeDivide(nk));

        }

        //IO.display(mu,"mu");

        RealMatrix[] sigma=new RealMatrix[k];
        for(int i=0;i<k;i++){
            sigma[i]=new Array2DRowRealMatrix(nData,nDim);
        }

        //r: nData x k
        RealMatrix r=new Array2DRowRealMatrix(nData,k);
        for(int i=0;i<nData;i++){
            for(int j=0;j<k;j++){
                //r[i][j]=Math.sqrt(R.getEntry(i,j));
                r.setEntry(i,j,Math.sqrt(R.getEntry(i,j)));
            }
        }

        //X: nDim x nData
        RealMatrix Xo=new Array2DRowRealMatrix(nDim,nData);
        for(int i=0;i<k;i++){
            RealVector cur_mu=mu.getColumnVector(i);

            for(int j=0;j<nData;j++){
                Xo.setColumnVector(j,X.getColumnVector(j).subtract(cur_mu));
            }

            RealVector cur_r=r.getColumnVector(i);
            for(int j=0;j<nDim;j++){
                Xo.setRowVector(j,Xo.getRowVector(j).ebeMultiply(cur_r));
            }


            //IO.display(Xo,"Xo");

            RealMatrix Sigma=Xo.transpose().preMultiply(Xo);
            for(int j=0;j<nDim;j++){
                for(int l=0;l<nDim;l++){
                    Sigma.setEntry(j,l,(Sigma.getEntry(j,l)/nk.getEntry(i))+1E-6);
                }
            }

            //IO.display(Sigma,"Sigma");

            sigma[i]=Sigma;
        }

        Model model=new Model();
        model.mu=mu;
        model.sigma=sigma;
        model.w=w;

        return model;

    }

    public static RealMatrix expectation(RealMatrix X,Model model,int nData,int nDim,int k,double[] llh,int iter){
        RealMatrix mu=model.mu;
        RealMatrix[] sigma=model.sigma;
        RealVector w=model.w;

        RealMatrix R1=new Array2DRowRealMatrix(nData,k);

        for(int i=0;i<k;i++){
            R1.setColumnVector(i,loggausspdf(X,mu.getColumnVector(i),sigma[i],nData,nDim));
        }

       for(int i=0;i<k;i++){
           w.setEntry(i,Math.log(w.getEntry(i)));
       }

       for(int i=0;i<nData;i++){
           R1.setRowVector(i,R1.getRowVector(i).add(w));
       }

        RealVector T=logsumexp(R1,nData,nDim,k);
        llh[iter]=StatUtils.sum(T.toArray())/nData;

        for(int i=0;i<k;i++){
            R1.setColumnVector(i,R1.getColumnVector(i).subtract(T));
        }

        for(int i=0;i<nData;i++){
            for(int j=0;j<k;j++){
                R1.setEntry(i,j,Math.exp(R1.getEntry(i,j)));
            }
        }

        return R1;
    }

    public static RealVector loggausspdf(RealMatrix X,RealVector mu,RealMatrix sigma,int nData,int nDim){

        RealMatrix X_mu=new Array2DRowRealMatrix(nDim,nData);
        for(int i=0;i<nData;i++){
            X_mu.setColumnVector(i,X.getColumnVector(i).subtract(mu));
        }
        //IO.display(X_mu,"X_mu");

        CholeskyDecomposition cholesky=new CholeskyDecomposition(sigma);
        RealMatrix U=cholesky.getL();
        //IO.display(U,"U");

        RealMatrix invU = new LUDecomposition(U).getSolver().getInverse();
        //IO.display(invU,"invU");

        RealMatrix Q=X_mu.preMultiply(invU);
        //IO.display(Q,"Q");

        RealVector q=new ArrayRealVector(nData);

        for(int i=0;i<nData;i++){
            q.setEntry(i,Q.getColumnVector(i).dotProduct(Q.getColumnVector(i)));
        }
        //IO.display(q,"q");

        double diagsum=0;
        for(int i=0;i<nDim;i++){
            diagsum+=Math.log(U.getEntry(i,i));
        }

        RealVector c=new ArrayRealVector(nData,nDim*Math.log(2*pi)+(2*diagsum));

        //IO.display(c,"c");

        RealVector y=c.add(q);
        y.mapMultiplyToSelf(-0.5);
        //IO.display(y,"y");
        return y;
    }

    public static RealVector logsumexp(RealMatrix R,int nData,int nDim,int k){
        RealMatrix y=new Array2DRowRealMatrix(nData,k);
        double rowmax=Double.NEGATIVE_INFINITY,rowsum=0;
        RealVector s=new ArrayRealVector(nData);

        for(int i=0;i<nData;i++){
            rowmax=StatUtils.max(R.getRow(i));
            rowsum=0;
            for(int j=0;j<k;j++){
                y.setEntry(i,j,Math.exp(R.getEntry(i,j)-rowmax));
                rowsum+=y.getEntry(i,j);
            }
            s.setEntry(i,Math.log(rowsum)+rowmax);
        }

        return s;

    }




    public static void main(String args[]){
        int nData=150,nDim=4;

        File file_in=new File("./data/LunarData/iris_dataset.txt");
        File init_label=new File("./data/LunarData/iris_label.txt");
        double[][] data=IO.readDoubleMat(file_in,nData,nDim);
        int[] init=IO.readIntVec(init_label,nData);

        //IO.display(data,"Iris dataset",nData,nDim);
        //IO.display(init,"Initial label",nData);

        int[] EMlabel=EM(data,init,nData,nDim,3);

        IO.display(EMlabel,"EM Label",150);
    }
}
