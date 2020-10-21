import java.io.*;
import java.util.*;
import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Queue;
import java.io.FilenameFilter;


public class FileIO
{
    public static int[][] readMatrixFromFile(String filename, String split, int ncol)
    {
        Scanner scn; 
        try 
        { scn = new Scanner(new File(filename)); }
        catch(IOException e)  
        { System.err.println("Invalid filename");  return null; }
 
        ArrayList<String> rows = new ArrayList<String>();

        while(true)
        {  
            try 
            { rows.add(scn.nextLine()); }
            catch(NoSuchElementException e) 
            { break; }
        }

        int n        = rows.size();
        if (ncol < 1) { ncol = n; }
        String[] bla = rows.toArray(new String[ncol]);
        int[][] res  = new int[n][ncol];


        for(int i=0; i< n; i++)
        {
            String[] numbers = bla[i].split(split);  

            for(int j = 0; j < ncol; j++)
            {       
                res[i][j] = Integer.parseInt(numbers[j]);
            }
        }

        return res;
    }

    public static void exportSubpaths(ArrayList<ArrayList<Integer>> results, String filename)
    {
        try 
        { 
            File file = new File(filename);
            
            if (!file.exists()) { file.createNewFile(); }
    
            FileWriter fw     = new FileWriter(file.getAbsoluteFile(), true);
            BufferedWriter bw = new BufferedWriter(fw);

            for(ArrayList<Integer> al : results)
            {
                for(Integer i : al)
                {
                    bw.write(i + " ");
                }
                bw.newLine();
            }
                    
            bw.close();
        }
        catch (IOException e) 
        {
            e.printStackTrace();
        }
    }
    public static void export(String results, String filename)
    {
        try 
        { 
            File file = new File(filename);
            
            if (!file.exists()) { file.createNewFile(); }
    
            FileWriter fw     = new FileWriter(file.getAbsoluteFile(), true);
            BufferedWriter bw = new BufferedWriter(fw);

            bw.write(results);
            bw.newLine();
            bw.close();
        }
        catch (IOException e) 
        {
            e.printStackTrace();
        }
    } 
}
