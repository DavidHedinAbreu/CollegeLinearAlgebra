/*
 * TestMatrix.java
 * 
 */

public class TestMatrix {
	
	public static void main (String[] args) {
		double[][] consistent = {{4,3,1},{-5,1,4},{1,-2,1}};  //     ,{2,9,5} (2, 1, -1)
		Matrix m1 = new Matrix(consistent);
		double[][] inconsistent = {{-3,1,0},{2,0,2},{-1,-2,-7},{0,-1,3}};  
		Matrix m2 = new Matrix(inconsistent);
		double[][] dependent = {{2,-1,1},{3,-4,-2},{5,-10,-8},{0,0,0}};
		Matrix m3 = new Matrix(dependent);
		double[][] twoBy2 = {{4,3},{-5,1}};  //     (2, 1, -1)
		Matrix m4 = new Matrix(twoBy2);

				
		System.out.println(m1.toString() );
		System.out.println(m2.toString() );
		Matrix m5 = m1.combineColumnwise(m2);
		System.out.println(m5.toString() );
		Matrix m6 = m5.subMatrixColumnwise(0,6);
		System.out.println(m6.toString() );

	}
}

