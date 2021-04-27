/*
 * TestMatrix.java
 * 
 */

public class TestMatrix {
	
	public static void main (String[] args) {
		double[][] consistent = {{4,3,1},{-5,1,4},{1,-2,1},{2,9,5}};  // (2, 1, -1)
		Matrix m1 = new Matrix(consistent);
		double[][] inconsistent = {{-3,1,0},{2,0,2},{-1,-2,-7},{0,-1,3}};
		Matrix m2 = new Matrix(inconsistent);
		double[][] dependent = {{2,-1,1},{3,-4,-2},{5,-10,-8},{0,0,0}};
		Matrix m3 = new Matrix(dependent);
				
		System.out.println(m1.reducedRowEchelonForm().toString() );
		System.out.println(m2.reducedRowEchelonForm().toString() );
		System.out.println(m3.reducedRowEchelonForm().toString() );
	}
}

