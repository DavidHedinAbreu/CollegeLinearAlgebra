/*
 * TestMatrix.java
 * 
 */

public class TestMatrix {
	
	public static void main (String[] args) {
		double[][] a = {{1,6,11,16,21},{2,7,12,17,22},{3,8,13,18,23},{4,9,14,19,24},{5,10,15,20,25}};
		Matrix m1 = new Matrix(a);
		double[][] b = {{1,2},{3,4}};
		Matrix m2 = new Matrix(b);
		System.out.println("\n" + m1.toString());
		System.out.println("\n" + m2.toString());
		
		System.out.println(m2.replace(0.5,1,1));
		
		// Matrix m3 = m1.multiply(m2);
		// System.out.println(m3.toString());
		//System.out.println(m1.magnitude() );
		//System.out.println("Determinant is " + m2.determinant() );
	}
}

