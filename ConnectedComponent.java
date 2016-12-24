

// Java program to fnd number of islands using Disjoint
// Set data structure.
import java.io.*;
import java.util.*;
 
public class ConnectedComponent
{
    public static void main(String[] args)throws IOException
    {
        int imgDim1=1500,imgDim2=300;
        int[][] connCompLabel=new int[imgDim1][imgDim2];
        
        Scanner sc;
        String root="J:\\Java programs\\Hyperspectral Workflow\\Data files\\";
        File file=new File(root+"ThresholdedImg.txt");
        int[][] data=new int[imgDim1][imgDim2];
        
        try{
            sc=new Scanner(file);
            sc.useDelimiter(",|\\n");
            
            for(int i=0;i<imgDim1;i++){
                for(int j=0;j<imgDim2;j++){
                    data[i][j]=Integer.parseInt(sc.next());
                }
            }
            sc.close();
            
        }
        catch(FileNotFoundException e){
            e.printStackTrace();
        }
        
        System.out.println("Number of Islands is: " +
                           countIslands(data,connCompLabel));
        
        HashSet<Integer> uniqueLabels=new HashSet<>();
        for(int i=0;i<imgDim1;i++){
            for(int j=0;j<imgDim2;j++){
                uniqueLabels.add(connCompLabel[i][j]);
            }            
        }
        
        
        //System.out.println("Unique exemplars:");
        HashMap<Integer,Integer> label=new HashMap<>();
        int count=1;
        for(int i:uniqueLabels){
            if(i!=0){                
                label.put(i, count);
                //System.out.print(i+"\t");
                count++;
            }            
        }
        
        for(int i=0;i<imgDim1;i++){
            for(int j=0;j<imgDim2;j++){
                if(connCompLabel[i][j]!=0){                
                    connCompLabel[i][j]=label.get(connCompLabel[i][j]);                    
                }
            }
        }
        
        
        File labeledImg=new File(root+"LabeledImg.txt");
        PrintWriter writer=new PrintWriter(labeledImg);
        System.out.println("Conneceted components label:");        
        for(int i=0;i<imgDim1;i++){
            //System.out.print(i+" : ");
            for(int j=0;j<imgDim2;j++){
                System.out.print(connCompLabel[i][j]+"\t");
                writer.print(connCompLabel[i][j]+" ");
            }
            System.out.println();
            writer.println();
        }
        writer.close();
        
        File connCompLabeledImg=new File(root+"ConnCompLabeledImg.png");
        ImageProc.imagesc(connCompLabeledImg, connCompLabel, uniqueLabels.size(), imgDim1, imgDim2);
        
        
        
     }
 
     // Returns number of islands in a[][]
     static int countIslands(int a[][],int[][] connCompLabel)
     {
        int n = a.length;
        int m = a[0].length;
 
        DisjointUnionSets dus = new DisjointUnionSets(n*m);
 
        /* The following loop checks for its neighbours
           and unites the indexes  if both are 1. */
        
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                // If cell is 0, nothing to do
                if (a[j][k] == 0)
                    continue;
 
                // Check all 8 neighbours and do a union
                // with neighbour's set if neighbour is 
                // also 1
                if (j+1 < n && a[j+1][k]==1){
                    dus.union(j*(m)+k, (j+1)*(m)+k);
                }
                if (j-1 >= 0 && a[j-1][k]==1){
                    dus.union(j*(m)+k, (j-1)*(m)+k);                    
                }
                if (k+1 < m && a[j][k+1]==1){
                    dus.union(j*(m)+k, (j)*(m)+k+1);                    
                }
                if (k-1 >= 0 && a[j][k-1]==1){
                    dus.union(j*(m)+k, (j)*(m)+k-1);                    
                }
                if (j+1<n && k+1<m && a[j+1][k+1]==1){
                    dus.union(j*(m)+k, (j+1)*(m)+k+1);                    
                }
                if (j+1<n && k-1>=0 && a[j+1][k-1]==1){
                    dus.union(j*m+k, (j+1)*(m)+k-1);                    
                }
                if (j-1>=0 && k+1<m && a[j-1][k+1]==1){
                    dus.union(j*m+k, (j-1)*m+k+1);                    
                }
                if (j-1>=0 && k-1>=0 && a[j-1][k-1]==1){
                    dus.union(j*m+k, (j-1)*m+k-1);                    
                }
            }
        }
        
        /*
        //Merge labels
        for(int j=0;j<n;j++){
            for(int k=0;k<m;k++){
                if(connCompLabel[j][k]==0)
                    continue;
                
                if (j+1 < n && connCompLabel[j+1][k]!=0 && connCompLabel[j+1][k]!=connCompLabel[j][k]){
                    dus.union(connCompLabel[j][k], connCompLabel[j+1][k]);
                    connCompLabel[j+1][k]=dus.find(connCompLabel[j][k]);
                }
                if (j-1 >= 0 && connCompLabel[j-1][k]!=0 && connCompLabel[j-1][k]!=connCompLabel[j][k]){
                    dus.union(connCompLabel[j][k], connCompLabel[j-1][k]);
                    connCompLabel[j-1][k]=dus.find(connCompLabel[j][k]);
                }
                if (k+1 < m && connCompLabel[j][k+1]!=0 && connCompLabel[j][k+1]!=connCompLabel[j][k]){
                    dus.union(connCompLabel[j][k], connCompLabel[j][k+1]);
                    connCompLabel[j][k+1]=dus.find(connCompLabel[j][k]);
                }
                if (k-1 >= 0 && connCompLabel[j][k-1]!=0 && connCompLabel[j][k-1]!=connCompLabel[j][k]){
                    dus.union(connCompLabel[j][k], connCompLabel[j][k-1]);
                    connCompLabel[j][k-1]=dus.find(connCompLabel[j][k]);
                }
                if (j+1<n && k+1<m && connCompLabel[j+1][k+1]!=0 && connCompLabel[j+1][k+1]!=connCompLabel[j][k]){
                    dus.union(connCompLabel[j][k], connCompLabel[j+1][k+1]);
                    connCompLabel[j+1][k+1]=dus.find(connCompLabel[j][k]);
                }
                if (j+1<n && k-1>=0 && connCompLabel[j+1][k-1]!=0 && connCompLabel[j+1][k-1]!=connCompLabel[j][k]){
                    dus.union(connCompLabel[j][k], connCompLabel[j+1][k-1]);
                    connCompLabel[j+1][k-1]=dus.find(connCompLabel[j][k]);
                }
                if (j-1>=0 && k+1<m && connCompLabel[j-1][k+1]!=0 && connCompLabel[j-1][k+1]!=connCompLabel[j][k]){
                    dus.union(connCompLabel[j][k], connCompLabel[j-1][k+1]);
                    connCompLabel[j-1][k+1]=dus.find(connCompLabel[j][k]);
                }
                if (j-1>=0 && k-1>=0 && connCompLabel[j-1][k-1]!=0 && connCompLabel[j-1][k-1]!=connCompLabel[j][k]){
                    dus.union(connCompLabel[j][k], connCompLabel[j-1][k-1]);
                    connCompLabel[j-1][k-1]=dus.find(connCompLabel[j][k]);
                }                
            }
        }
        */
        // Array to note down frequency of each set
        int[] c = new int[n*m];
        int numberOfIslands = 0;
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                if (a[j][k]==1)
                {
 
                    int x = dus.find(j*m+k);
                    connCompLabel[j][k]=x;
                    // If frequency of set is 0, 
                    // increment numberOfIslands
                    if (c[x]==0)
                    {
                        numberOfIslands++;
                        c[x]++;
                    }
 
                    else
                        c[x]++;
                }
            }
        }
        return numberOfIslands;
    }
}
 
// Class to represent Disjoint Set Data structure
class DisjointUnionSets
{
    int[] rank, parent;
    int n;
 
    public DisjointUnionSets(int n)
    {
        rank = new int[n];
        parent = new int[n];
        this.n = n;
        makeSet();
    }
 
    void makeSet()
    {
        // Initially, all elements are in their
        // own set.
        for (int i=0; i<n; i++)
            parent[i] = i;
    }
 
    // Finds the representative of the set that x
    // is an element of
    int find(int x)
    {
        if (parent[x] != x)
        {
            // if x is not the parent of itself,
            // then x is not the representative of
            // its set.
            // so we recursively call Find on its parent
            // and move i's node directly under the
            // representative of this set
            return find(parent[x]);
        }
 
        return x;
    }
 
    // Unites the set that includes x and the set
    // that includes y
    void union(int x, int y)
    {
        // Find the representatives (or the root nodes)
        // for x an y
        int xRoot = find(x);
        int yRoot = find(y);
 
        // Elements are in the same set, no need
        // to unite anything.
        if (xRoot == yRoot)
            return;
 
        // If x's rank is less than y's rank
        // Then move x under y  so that depth of tree
        // remains less
        if (rank[xRoot] < rank[yRoot])
            parent[xRoot] = yRoot;
 
        // Else if y's rank is less than x's rank
        // Then move y under x so that depth of tree
        // remains less
        else if(rank[yRoot]<rank[xRoot])
            parent[yRoot] = xRoot;
 
        else  // Else if their ranks are the same
        {
            // Then move y under x (doesn't matter
            // which one goes where)
            parent[yRoot] = xRoot;
 
            // And increment the the result tree's
            // rank by 1
            rank[xRoot] = rank[xRoot] + 1;
        }
    }
}