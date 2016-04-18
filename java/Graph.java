import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Stack;
import java.util.concurrent.*;


public class Graph 
{
    private int N;
    public int size; 
	public ArrayList<Node> nodes;
	public ArrayList<Integer> al;
	public ArrayList<ArrayList<Integer>> results;
	private int[][] adjMatrix; 	  
    long startTime;

    public Graph(int[][] matrix, long sTime)
    {
        size    = 0;
    	nodes   = new ArrayList<Node>();
    	al      = new ArrayList<Integer>();
        results = new ArrayList<ArrayList<Integer>>();
        startTime = sTime;
        createAdjMatrix(matrix);
    }

    public void addNode(Node n)
    {
		//
        nodes.add(n);
    }

	public ArrayList<ArrayList<Integer>> getResults()
	{
        //
		return results;
	}

	public void createAdjMatrix(int[][] matrix)
	{
		if(adjMatrix==null)
		{
            for(int j=0; j<matrix.length; j++)
            {
                addNode(new Node(j));
            }
			N         = nodes.size();
			adjMatrix = new int[N][N];
			adjMatrix = matrix;
		}
	}
	
	public void dfs(int src, int dst, int a, int b, int limit) 
	{
        float t = TimeUnit.NANOSECONDS.toSeconds(System.nanoTime() - startTime);
        if (t > limit) return;

        al.add(src);
        size++;
        nodes.get(src).visited = true;
        if ( src == dst ) 
        {
            if ( size >= a)
            {
                ArrayList<Integer> r = new ArrayList<Integer>();
                for(Integer i : al) { r.add(i); }
                results.add(r);
            }
            
            return;
        }

        if ( size >= b ) return;

        for ( int i = 0; i < N; i++ ) 
        {
            if ( adjMatrix[src][i] == 1 )
            {
                if ( nodes.get(i).visited == false )
                {
                    dfs(i, dst, a, b, limit);
                    nodes.get(i).visited = false;   
                    size--;
                    al.remove(size);
                }
            }
        }
    }
}

class Node 
{
    public int label;
    public boolean visited=false;
    public Node(int l)
    {
        this.label = l;
    }
}
