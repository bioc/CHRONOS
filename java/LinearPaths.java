import java.util.Queue;
import java.util.ArrayList;
import java.io.IOException;
import java.util.concurrent.*;

public class LinearPaths implements Runnable
{
	static int a;
	static int b;
	static int threshold;
	static int limit;
	static String matFile;
	static String parFile;
	static String outFile;
	static String logFile;

	private static int[][] matrix;
	private static int[][] validPairs;
	private static ArrayList<ArrayList<Integer>> results;
	static String[] queue;
	
	public static int getLinearPaths(String inputFile, String validPairsFile, String outFile, boolean writeToDisk, int a, int b, int threshold, int limit)
	{
		queue = new String[2];

		// Read the adjacency matrix and thw valid source/destination pairs from file 
		matrix     = FileIO.readMatrixFromFile(inputFile, ",", 0);

		// Read valid source/destination pairs
		validPairs = FileIO.readMatrixFromFile(validPairsFile, "	", 2);
		if (validPairs == null) { return 0; }
		int rows = validPairs.length;
		int cols = validPairs[0].length;
		int subpaths = 0;

		long startTime = System.nanoTime();
		results   = new ArrayList<ArrayList<Integer>>();

		for (int i=0; i < rows; i++)
		{
			// Create the graph based on the given adjacency matrix
			Graph g = new Graph(matrix, startTime);
			// Find all possible paths between source and destination using dfs			
			int src = validPairs[i][0];
			int dst = validPairs[i][1];
			g.dfs(src-1, dst-1, a, b, limit);
			int writeBlock = 50;

			// Add patrial results for current source/destination pair to results
			for (ArrayList<Integer> r : g.getResults())
			{
				// Translate results to one based alphabet.
				for (int idx = 0; idx < r.size(); idx++) { r.set(idx, r.get(idx) + 1); }
				results.add(r);
								
				if(results.size() == writeBlock)
				{
					if (writeToDisk) 
					{
						// Output results to disk
						FileIO.exportSubpaths(results, outFile);
						results.clear();
					}
				}
				subpaths = subpaths + 1;
				if (subpaths > threshold) { break; }
			}
			if (subpaths > threshold) { break; }
		}
		
		if (writeToDisk) { FileIO.exportSubpaths(results, outFile); }
		subpaths = subpaths + results.size();

		return subpaths;
	}

	public static void main(String[] args)  throws Exception 
	{
	 	String s        = args[0];
	 	String org      = args[1];
	 	String baseDir  = args[2];

	 	matFile   = baseDir + "mat/" + org + "/" + s;
	 	parFile   = baseDir + "par/" + org + "/" + s;
	 	outFile   = baseDir + "out/" + org + "/" + s;
	 	logFile   = baseDir + "log/" + org + "/" + s;
	 	a         = Integer.parseInt(args[3]);
	 	b         = Integer.parseInt(args[4]);
	 	threshold = Integer.parseInt(args[5]);
	 	limit     = Integer.parseInt(args[6]);
	 	LinearPaths lp = new LinearPaths();

	 	ExecutorService executor = Executors.newSingleThreadExecutor();
	 	Callable<String> callable = new Callable<String>() {
	 	    public String call() 
	 	    {
	 	    	try {
	 	    	    lp.run();
	 	    	} finally {
	 	    	    String out = lp.queue[0] + "-" + lp.queue[1];
	 	    	    return out;
	 	    	}
	 	    }
	 	};
	 	Future<String> future = executor.submit(callable);
	 	String out = future.get();
	 	executor.shutdown();
	 	FileIO.export(out, logFile);
	}

	public void run() 
	{
		long startTime = System.nanoTime();
		int n = getLinearPaths(matFile, parFile, outFile, true, a, b, threshold, limit);
		long stopTime = System.nanoTime();
		float t = (int)TimeUnit.NANOSECONDS.toSeconds(stopTime - startTime);
		this.queue[0] = Integer.toString(n);
		this.queue[1] = Float.toString(t);
	}
}