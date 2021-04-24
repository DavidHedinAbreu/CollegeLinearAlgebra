/*
 * Matrix
 * Author: David Hedin-Abreu
 * 
 * Systems of linear equations
Row reduction and echelon forms
Matrix operations, including inverses
Block matrices
Linear dependence and independence
Subspaces and bases and dimensions
Orthogonal bases and orthogonal projections
Gram-Schmidt process
Linear models and least-squares problems
Determinants and their properties
Cramer's Rule
Eigenvalues and eigenvectors
Diagonalization of a matrix
Symmetric matrices
Positive definite matrices
Similar matrices
Linear transformations
Singular Value Decomposition
 */
 
import java.util.Vector;

public class Matrix{
	
	double[][] matrixArray;
	
	public Matrix (double[][] a) 
	{
		matrixArray = a;
	}
	
	public double magnitude()  // complete
	{
		double[][] a = this.matrixArray;
		double output=0;
		if( a.length != 2 )
		{
			System.out.println("Could not find magnitude.");
		}
		else
		{
			for(int term = 0; term < a[0].length; term++)
			{
				output += Math.pow(a[0][term] - a[1][term],2);
			}
			output = Math.sqrt(output);
		}
		return output;
	}
	
	public double angleBetween(Matrix m2) // complete
	{
		double[][] a = this.matrixArray;
		double[][] b = m2.matrixArray;
		double output=0;
		if( !(a.length == 2 && b.length == 2) )
		{
			System.out.println("Could not find magnitude.");
		}
		else
		{
			output = Math.acos( this.dotProduct(m2) / ( this.magnitude() * m2.magnitude() ) );
		}
		return output;
	}
	
	public double dotProduct(Matrix m2)  // complete
	{
		double[][] a = this.matrixArray;
		double[][] b = m2.matrixArray;
		double output = 0;
		if( !(a.length == 2 && a[0].length == b[0].length) )
		{
			System.out.println("Could not find dot-product.");
		}
		else
		{
			for(int term = 0; term < a[0].length; term++)
			{
				output += (a[1][term] - a[0][term])*(b[1][term] - b[0][term]);
			}
		}
		return output;
	}
	
	public Matrix interpretRREF()  // must handle inconsistent systems
	{
		// given row starts at 0.  apply swapScaleReplace. advance row by 1, repeat, until at bottom row of matrix.
		// https://textbooks.math.gatech.edu/ila/parametric-form.html
	}

	public Matrix rowReduceToRREF()   // completed: does not interpret RowReducedEchelonForm to provide solution(s).
	{
		double[][] output = this.matrixArray.clone();
		double scaleFactor;
		double temp;
		for(int topRow = 0, leftColumn = 0; topRow < output[0].length; topRow++)
		{
			// If the leftmost column in the given row is 0, find the row with the the left-most entry of all rows at or lower than the row given. 
			// If there are more than one, find the first. Swap that row for the row given, so the leftmost entry is in the given row.
			// After this routine the pivot for the row is at output[leftColumn][row].
			if(output[leftColumn][topRow] == 0)
			{
				for(int row = topRow; row < output[0].length; row++)
				{
					if(output[leftColumn][row] != 0)  
					{
						for(int col = 0; col < output.length; col++)
						{
							temp = output[col][topRow];
							output[col][topRow] = output[col][row];
							output[col][row] = temp;
						}

					}
				}
			}
			// Scale the row given so the pivot (leftmost entry) is 1.  
			scaleFactor = 1.0/output[leftColumn][topRow];
			for(int col = 0; col < output.length; col++)
			{
				output[col][topRow] = scaleFactor*output[col][topRow];
			}
			// Row-replace any lower row with a leftmost entry in the same column, so that entry becomes 0.  
			for(int row = topRow + 1; row < output[0].length; row++)
			{
				if(output[leftColumn][row] != 0)
				{
					for (int col = leftColumn; col < output.length; col++)
					{
						scaleFactor = -1*output[leftColumn][row];
						output[col][row] = scaleFactor*output[col][row];
					}
				}
			}
		}
		return new Matrix(output);
	}

	public Matrix swap(int rowA, int rowB)   // completed: NOT a helpter for "rowReduce"
	{
		double[][] a = this.matrixArray;
		double[][] output = a.clone();
		double temp;
		for(int column = 0; column < a.length; column++)
		{
			temp = output[column][rowA];
			output[column][rowA] = a[column][rowB];
			output[column][rowB] = temp;
		}
		return new Matrix(output);
	}

	public Matrix scale(double scalar, int row)   // completed: NOT a helper method to "rowReduce"
	{
		double[][] output = this.matrixArray.clone();
		for(int column = 0; column < output.length; column++)
		{
			output[column][row] = scalar*output[column][row];
		}
		return new Matrix(output);
	}

	public Matrix replace(double scalar, int rowScaled, int rowReplaced)   // completed: NOT a helper method to "rowReduce". Note if rowScaled=rowReplaced, this is just a scale operation.
	{
		double[][] output = this.matrixArray.clone();
		for(int column = 0; column < output.length; column++)
		{
			output[column][rowReplaced] = output[column][rowReplaced] + scalar*output[column][rowScaled];
		}
		return new Matrix(output);
	}

	public Matrix invert()
	{
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}

	public Matrix diagonalize()
	{
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}

	public Matrix crossProductOf(Matrix a, Matrix b)
	{
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}
	
	public Matrix add(Matrix a, Matrix b)
	{
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}
	
	public Matrix multiply(Matrix m2)
	{
		double[][] a = this.matrixArray;
		double[][] b = m2.matrixArray;
		double[][] output;
		double value;
		if( !(a.length == b[0].length && a[0].length == b.length) )
		{
			System.out.println("Could not multiply matrices.");
			output = null;
		}
		else
		{
			output = new double[a[0].length][b.length];
			for(int row = 0; row < a[0].length; row++)
			{
				for(int col = 0; col < b.length; col++)
				{
					value = 0; 
					for(int pos = 0; pos < a.length; pos++) 
						value += a[pos][row] * b[col][pos]; 
					output[col][row] = value; 
				}
			}
		}
		return new Matrix(output);
	}
	
	public Matrix add(double a, Matrix b)
	{
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}

	public Matrix multiply(double a, Matrix b)
	{
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}

	public Matrix divide(double a)
	{
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}
		
	public static Matrix identity(int dimension)
	{
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}
	
	public CoFactor cofactor(int column)  // complete (assistant method to "determinant")
	{
		double s;
		Matrix minorMatrix;
		double[][] a = this.matrixArray;
		if( !(a.length >= 3 && a[0].length >= 3) )
		{
			System.out.println("Could not make CoFactor.");
			s = 0;
			double[][] bogus = {};
			minorMatrix = new Matrix(bogus);
		}
		else
		{
			double[][] b = new double[a.length-1][a[0].length-1];
			s = a[column][0];
			for(int colA = 0, colB = 0; colA < a.length; colA++)
			{
				if(colA != column)
				{
					for(int row = 1; row < a[0].length; row++)
					{
						b[colB][row-1] = a[colA][row];
					}
					colB++;
				}
			}
			minorMatrix = new Matrix(b);
		}

		return new CoFactor(s, minorMatrix);
	}
	
	public double determinant()  // complete
	{
		double[][] a = this.matrixArray;
		double output = 0;
		CoFactor parentCoFactor;
		int determinantTermSign = 1;
		if( a.length == 2 && a[0].length == 2 )
		{
			output += determinantTermSign*a[0][0]*a[1][1] - a[1][0]*a[0][1];
		}
		else
		{
			for(int column = 0; column < a.length; column++)
			{
				parentCoFactor = this.cofactor(column);
				output += determinantTermSign*parentCoFactor.scalar*parentCoFactor.minorMatrix.determinant();
				determinantTermSign *= -1;
			}
		}
		return output;
	}
	
	public Matrix inverse()
	{
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}
		
	public String toString()  // complete
	{
		double[][] a = this.matrixArray;
		String output = "";
		if(a == null)
			output = "null";
		else
		{
			for(int row = 0; row < a[0].length; row++)
			{
				for(int column = 0; column < a.length; column++)
				{
					if(a[column][row]%1 != 0)
						output += a[column][row] + " ";
					else 
						if(a[column][row] < 10)
							output += "   " + (int)a[column][row] + " ";
						else if(a[column][row] < 100)
							output += "  " + (int)a[column][row] + " ";
						else if(a[column][row] < 1000)
							output += " " + (int)a[column][row] + " ";
				}
				output += "\n";
			}
		}
		return output;
	}
	
	public class CoFactor  // complete
	{
		double scalar;
		Matrix minorMatrix;
		public CoFactor(double s, Matrix m)
		{
			this.scalar = s;
			this.minorMatrix = m;
		}
		
		public String toString()
		{
			double[][] a = this.minorMatrix.matrixArray;
			String output = "";
			if(a == null)
				output = "null";
			else
			{
				for(int row = 0; row < a[0].length; row++)
				{
					for(int column = 0; column < a.length; column++)
					{
						if(a[column][row]%1 != 0)
							output += a[column][row] + " ";
						else 
							if(a[column][row] < 10)
								output += "   " + (int)a[column][row] + " ";
							else if(a[column][row] < 100)
								output += "  " + (int)a[column][row] + " ";
							else if(a[column][row] < 1000)
								output += " " + (int)a[column][row] + " ";
					}
					output += "\n";
				}
			}
			return this.scalar + ", \n" + output;
		}
	}	
}

