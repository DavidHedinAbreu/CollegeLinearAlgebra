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
		double[][] a = this.matrixArray.clone();
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
		double[][] a = this.matrixArray.clone();
		double[][] b = m2.matrixArray.clone();
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
		double[][] a = this.matrixArray.clone();
		double[][] b = m2.matrixArray.clone();
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
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}

	public Matrix reducedRowEchelonForm()   // completed: does not interpret RowReducedEchelonForm to provide solution(s).
	{
		Matrix outputMatrix = this; // makes a clone
		for(int topRow = 0, leftColumn = 0; topRow < outputMatrix.matrixArray[0].length && leftColumn < outputMatrix.matrixArray.length; leftColumn++) 
		{
			outputMatrix.swap(leftColumn, topRow);
			
			if( outputMatrix.scale(leftColumn, topRow) )  // skip rowReduceDown if there were no non-zero values in leftColumn.
			{
				outputMatrix.rowReduceDown(leftColumn, topRow);
				topRow++;   // move topRow down if and only if the scale and rowReduceDown steps took place. Otherwise stay in the same row.
			}
		}
		outputMatrix.rowReduceUp();
		return outputMatrix;
	}
	
	// If the leftmost column in the given row is 0, find the row with the the left-most entry of all rows at or lower than the row given. 
	// If there are more than one, find the first. Swap that row for the row given, so the leftmost entry is in the given row.
	// After this routine the pivot for the row is at output[leftColumn][row].
	public void swap(int leftColumn, int topRow)   // completed: a helper for "rowReduceDown" 
	{
		double[][] output = this.matrixArray.clone();
		double temp;
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
		this.matrixArray = output;
	}

	public boolean scale(int leftColumn, int topRow)   // completed:  a helper method to "rowReduceDown"
	{
		double[][] output = this.matrixArray.clone();
		boolean doRowReduceDown = true;
		// if no row in leftColumn is has nonzero entry, don't scale row because scalefactor would be Infinite/NaN, pass flag so don't rowReduceDown either.
		if(output[leftColumn][topRow] == 0)  
		{
			doRowReduceDown = false;  
		}
		else
		{
			double scaleFactor = 1.0/output[leftColumn][topRow];
			for(int col = 0; col < output.length; col++)
			{
				output[col][topRow] = scaleFactor*output[col][topRow];
			}

			this.matrixArray = output;
		}
		return doRowReduceDown;
	}

	public void rowReduceDown(int leftColumn, int topRow)  // completed: helper method for RREF
	{
		double[][] output = this.matrixArray.clone();
		double scaleFactor;
		for(int row = topRow + 1; row < output[0].length; row++)
		{
			if(output[leftColumn][row] != 0)
			{
				scaleFactor = -1*output[leftColumn][row];
				for (int col = leftColumn; col < output.length; col++)
				{
					output[col][row] += scaleFactor*output[col][topRow];
				}
			}
		}
		
		this.matrixArray = output;
	}
	
	public void rowReduceUp()  // completed: helper method for RREF
	{
		// Start at bottomRow = length - 1, leftColumn = 0. If the value at that position is zero, advance leftColumn, checking each time whether the value there is zero.
		// If all the positions in the row are 0, decrease bottomRow by 1, reset leftColumn to 0, and perform zero checks on the row above, left to right.
		// When a nonzero element is found, it should be 1.  Row-replace up using that row, so all positions above bottomRow at leftColumn have value 0, then begin zero 
		// check on the row above.
		double[][] output = this.matrixArray.clone();
		double scaleFactor;
		for(int bottomRow = output[0].length - 1; bottomRow >= 0; bottomRow--)
		{
			for(int leftColumn = 0; leftColumn < output.length; leftColumn++)
			{
				if(output[leftColumn][bottomRow] !=0)  // if we found a pivot
				{
					for(int row = bottomRow - 1; row >= 0; row--)  // row-replace up
					{
						scaleFactor = -1*output[leftColumn][row];

						for(int col = 0; col < output.length; col++)
						{
							output[col][row] += scaleFactor*output[col][bottomRow];
						}
					} 
					leftColumn = output.length;  // stop searching for pivot
				}
			}
		}
		this.matrixArray = output;
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
		double[][] a = this.matrixArray.clone();
		double[][] b = m2.matrixArray.clone();
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
		double[][] a = this.matrixArray.clone();
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
		double[][] a = this.matrixArray.clone();
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
		double[][] a = this.matrixArray.clone();
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
			double[][] a = this.minorMatrix.matrixArray.clone();
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

		/* Set default starting point for row-replacement up, for the case of a consistent system.
		int bottomRow = output[0].length - 1, rightColumn = output.length - 2;  
		// Check to see whether there is no solution / the system is inconsistent: every entry on the bottom row zero, except for the last which is non-zero pivot.
		// If so, scale the last entry to 1, then row-replace up in the last column.
		if(output[output.length-1][output[0].length-1] == 1 && output[output.length-2][output[0].length-1] == 0) // If inconsistent
		{	
			output[output.length-1][output[0].length-1] = 1;
			for(int row = output[0].length-2; row >= 0; row--)
			{
				output[output.length-1][row] = 0;
			}
		}
		// In the case of a dependent or inconsistent system, set starting place for row-replacement to the first non-zero entry "1" on the left, in the second to last row.  
		// Begin row replacement from that pivot point instead of the default.
		if(output[output.length-2][output[0].length-1] == 0)  
		{
			bottomRow = output[0].length - 2;
			for(rightColumn = 0; rightColumn < output.length - 1 && output[rightColumn][bottomRow] == 0; rightColumn++)
			{
				// no code block, purpose is to advance rightColumn to the right place.
			}
		}
		// Start at bottom row, rightmost column (left of solution vector).  Row-replace any higher row with a rightmost column entry, so that entry becomes zero.  
		// If the bottom has zero as bottom-rightmost entry, move up until the row does have a nonzero rightmost column entry.
		for(; bottomRow >= 0; bottomRow--, rightColumn--)
		{
			while(output[rightColumn][bottomRow] == 0)
			{
				bottomRow--;
			}
			for(int row = bottomRow - 1; row >= 0; row--)
			{
				if(output[rightColumn][row] != 0)
				{
					scaleFactor = -1*output[rightColumn][row];
					for (int col = rightColumn; col <= output.length - 1; col++)
					{
						output[col][row] += scaleFactor*output[col][bottomRow];
					}
				}
			}
			whatTF(output);
		}
		return new Matrix(output);
	} */
