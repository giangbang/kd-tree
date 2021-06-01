//Java implementation of kd tree
import java.util.List;
import java.util.ArrayList;
import java.util.PriorityQueue;
import java.lang.Double;

public class kdTree {
	
	private class treeNode {
		public double[] val;
		public treeNode left, right;
		
		public treeNode(double[] val) {
			this.val = val;
			this.left = null;
			this.right = null;
		}
		
		// Compare the data point of a tree node to a given data point
		public boolean equals(double[] p) {
			for (int i = 0; i < p.length; ++i) {
				if (p[i] != this.val[i]) return false;
			}
			return true;
		}
		
		public boolean isLeaf() {
			return this.left == null && this.right == null;
		}
	}
	
	private int k = 0;
	private treeNode root = null;
	private int size = 0;
	
	public kdTree(int k) {
		this.k = k;
	}
	
	
	// ==========================================================
	// add the point to the set 
	// ==========================================================
	public void insert(double[] p) {
		check(p);
		++this.size;
		this.root = this.insertHelper(root, p, 0);
	}
	
	private treeNode insertHelper(treeNode node, double[] p, int depth) {
		if (node == null) node = new treeNode(p);
		else {
			int d = depth % this.k;
			if (node.val[d] >= p[d]) 
				node.left = insertHelper(node.left, p, depth+1);
			else 
				node.right = insertHelper(node.right, p, depth+1);
		}
		return node;
	}
	
	
	// ==========================================================
	// does the set contain given point ? 
	// ==========================================================
	public boolean contains(double[] queryPoint) {
		check(queryPoint);
		return this.containsHelper(root, queryPoint, 0);
	}
	
	private boolean containsHelper(treeNode node, double[] p, int depth) {
		if (node == null) return false;
		int d = depth % this.k;
		if (node.equals(p)) return true;
		if (node.val[d] >= p[d]) 
			return containsHelper(node.left, p, depth+1);
		else
			return containsHelper(node.right, p, depth+1);
	}
	
	
	// ==========================================================
	// range search all data points inside min and max
	// ==========================================================
	public List<double[]> rangeSearch(double[] min, double[] max) {
		check(min); check(max);
		
		List<double[]> res = new ArrayList<>();
		rangeSearchHelper(res, root, min, max, 0);
		return res;
	}
	
	private void rangeSearchHelper(List res, treeNode node, double[] min, 
			double[] max, int depth) {
		if (node == null) return;
		int d = depth % this.k;
		if (node.val[d] < max[d])   // queried range overlaps with right subtree
			rangeSearchHelper(res, node.right, min, max, depth+1);
		if (node.val[d] >= min[d])  // queried range overlaps with left subtree
			rangeSearchHelper(res, node.left, min, max, depth+1);
		if (inside(node.val, min, max)) // current node is inside the queried range
			res.add(node.val); 
	}
	
	
	// ==========================================================
	// k nearest neighbors
	// ==========================================================
	public List<double[]> nearestNeighbors(int k, double[] p) {
		check(p);
		assert(k > 0);
		
		PriorityQueue<qNode> top_nearest = new PriorityQueue<>();
		nearestHelper(top_nearest, root, p, k, 0);
		
		List<double[]> ret = new ArrayList<>();
		for (qNode node : top_nearest) {
			ret.add(node.data);
		}
		return ret;
	}
	
	private void nearestHelper(PriorityQueue<qNode> q, treeNode node, 
			double[] p, int k, int depth) {
		if (node == null) return;
		int d = depth % this.k;
		boolean go_left = true;
		if (node.val[d] >= p[d]) 
			nearestHelper(q, node.left, p, k, depth+1); 
		else {
			nearestHelper(q, node.right, p, k, depth+1);
			go_left = false;
		}
		
		double dis = squaredDistance(node.val, p);
		if (q.size() < k) q.add(new qNode(node.val, dis));
		else if (dis < q.peek().sdis) {
			q.poll();
			q.add(new qNode(node.val, dis));
		}
		
		if (q.size() < k || squaredDistance(node.val[d], p[d]) < q.peek().sdis) {
			if (go_left) 
				nearestHelper(q, node.right, p, k, depth+1);
			else
				nearestHelper(q, node.left, p, k, depth+1);
		}
	}
	
	// placeholder in priority queue
	private static class qNode implements Comparable<qNode>{
		double[] data;
		double sdis = 0; // squared distance
		
		public qNode(double[] p, double dis) {
			this.data = p;
			this.sdis = dis;
		}
		
		public int compareTo(qNode that) {
			return Double.compare(that.sdis, this.sdis);
		}
	}
	
	
	// ==========================================================
	// other utility functions
	// ==========================================================
	
	public int size() {
		return this.size;
	}
		
	public boolean isEmpty() {
		return root == null;
	}
	
	private static double squaredDistance(double x, double y) {
		double dif = x - y;
		return dif * dif;
	}
	
	private static double squaredDistance(double[] x, double[] y) {
		double res = 0.;
		for (int i = 0; i < x.length; ++i) {
			double dif = x[i] - y[i];
			res += dif * dif;
		}
		return res;
	}
	
	// check if the data point p is inside the range(min, max)
	private static boolean inside(double[] p, double[] min, double[] max) {
		for (int i = 0; i < p.length; ++i) {
			if (p[i] < min[i] || p[i] > max[i]) 
				return false;
		}
		return true;
	}
	
	// sanity check for the number of dimension
	private void check(double[] p) {
		if (p.length != this.k) 
			throw new 
			IllegalArgumentException("Number of dimension does not match with k");
	}
	
	// unit testing of the methods
	public static void main(String[] args) {
		
		kdTree tree = new kdTree(2);
		// insert and contains
		tree.insert(new double[]{1.2, 1.4});
		tree.insert(new double[]{1.4, 2.5});
		System.out.println(tree.contains(new double[]{1.4, 2.3}));
		System.out.println(tree.contains(new double[]{1.4, 2.5}));
		// >> false
		// >> true
		
		
		tree.insert(new double[]{1.3, 1.8});
		tree.insert(new double[]{1.6, 1.});
		// range search
		System.out.println();
		List<double[]> range = tree.rangeSearch(new double[]{1,1}, new double[]{2,2});
		System.out.println(range.size());
		for (double[] i : range) {
			for (double j : i) System.out.print(j + " ");
			System.out.println();
		}
		// >> 3
		// >> 1.6 1.0
		// >> 1.3 1.8
		// >> 1.2 1.4
		
		
		// k nearest neighbors
		tree.insert(new double[]{1.3, 1.});
		List<double[]> nei = tree.nearestNeighbors(2, new double[]{1.5,2});
		System.out.println(nei.size());
		for (double[] i : nei) {
			for (double j : i) System.out.print(j + " ");
			System.out.println();
		}
		// >> 2
		// >> 1.4 2.5
		// >> 1.3 1.8
	}
}