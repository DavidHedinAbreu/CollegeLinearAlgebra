/*
 * Matrix
 * Author: David Hedin-Abreu
 * 
 * Systems of linear equations
Row reduction and echelon forms  DONE
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
	
	private double[][] matrixArray;
	
	public Matrix (double[][] a) 
	{
		matrixArray = a;
	}
	
	public Matrix combineColumnwise(Matrix m2)  // complete
	{
		double[][] a = this.matrixArray;
		double[][] b = m2.matrixArray;
		Matrix m;
		if(!(a[0].length == b[0].length))
		{
			System.out.println("Matrices or vectors could not be combined.  They must have the same number of rows.");
			m = new Matrix(null);
		}
		else
		{
			double[][] output = new double[a.length + b.length][a[0].length];
			int numRows = a[0].length;
			for(int col = 0; col < a.length; col++)
			{
				for(int row = 0; row < numRows; row++)
				{
					output[col][row] = a[col][row];
				}
			}
			for(int col = a.length; col < a.length + b.length; col++)
			{
				for(int row = 0; row < numRows; row++)
				{
					output[col][row] = b[col - a.length][row];
				}
			}
			m = new Matrix(output);
		}
		return m;
	}
	
	public Matrix subMatrixColumnwise(int colA, int colB)  // complete
	{
		double[][] a = this.matrixArray;
		Matrix m;
		if(colA > colB || colA < 0 || colB >= a.length)
		{
			System.out.println("Could not find submatrix taken columnwise. Illegal bounds.");
			m = new Matrix(null);
		}
		else 
		{
			double[][] output = new double[colB - colA+1][a[0].length];
			int numRows = a[0].length;
			for(int col = colA, pos = 0; col <= colB && pos < colB - colA + 1; col++, pos++)
			{
				for(int row = 0; row < numRows; row++)
				{
					output[pos][row] = a[col][row];
				}
			}
			m = new Matrix(output);
		}
		return m;
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
	
	public double angle(Matrix m2) // complete
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
	
	public String interpretRREF()  // must handle inconsistent systems
	{
		// give rank, classification (consistent with solution, dependent with parametric solution, or inconsistent)
		// https://textbooks.math.gatech.edu/ila/parametric-form.html
		double[][] interpret = this.matrixArray.clone();
		String interpretation = "";
		int row = interpret[0].length - 1, pivot = -1;
		for(int col = interpret.length - 2; col >= 0; col--)
		{
			if(interpret[col][row] == 1)
				pivot = row;
		}
		if(pivot == -1 && interpret[interpret.length - 1][row] == 1)
			interpretation = "Inconsistent";
		else if(pivot == -1 && interpret[interpret.length - 1][row] == 0)
			interpretation = "Dependent";
		else
			interpretation = "Consistent";
		
		if(interpretation.equals("Dependent"))
		{
			interpretation += "\n";
			for(row = 0; row < interpret[0].length - 1; row++)
			{
				for(int col = 0; col < interpret.length - 1; col++)
				{
					if(interpret[col][row] != 0)
					{
						if( ! interpretation.substring(interpretation.length() - 1).equals("\n"))
							interpretation += " +";
						if(correctBinary(interpret[col][row]) != 1)
							interpretation += " " + correctBinary(interpret[col][row]);
						interpretation += "X" + col;
					}
				}
				interpretation += " = " + correctBinary(interpret[interpret.length - 1][row]) + "\n";
			}
		}
		if(interpretation.equals("Consistent"))
		{
			interpretation += "\n";
			for(row = 0; row < interpret[0].length ; row++)
			{
				for(int col = 0; col < interpret.length - 1; col++)
				{
					if(interpret[col][row] != 0)
					{
						if( ! interpretation.substring(interpretation.length() - 1).equals("\n"))
							interpretation += " +";
						if(correctBinary(interpret[col][row]) != 1)
							interpretation += " " + correctBinary(interpret[col][row]);
						interpretation += "X" + col;
					}
				}
				interpretation += " = " + correctBinary(interpret[interpret.length - 1][row]) + "\n";
			}
		}
		
		return interpretation;  // new Matrix(matrixContents);
	}

	public Matrix rREF()   // completed: does not interpret RowReducedEchelonForm to provide solution(s).
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

	public void swap(int leftColumn, int topRow)   // completed: a helper for "rREF" 
	{
		// If the leftmost column in the given row is 0, find the row with the the left-most entry of all rows at or lower than the row given. 
		// If there are more than one, find the first. Swap that row for the row given, so the leftmost entry is in the given row.
		// After this routine the pivot for the row is at output[leftColumn][row].
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

	public boolean scale(int leftColumn, int topRow)   // completed:  a helper method to "rREF"
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

	public void rowReduceDown(int leftColumn, int topRow)  // completed: helper method for "rREF"
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
	
	public void rowReduceUp()  // completed: helper method for "rREF"
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
	
	public Matrix diagonalize()
	{
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}

	public Matrix crossProduct(Matrix m2)
	{
		double[][] matrixContents = {{1, 0}, {0, 1}};
		return new Matrix(matrixContents);
	}
	
	public Matrix add(Matrix m2)  // completed
	{
		double[][] a = this.matrixArray.clone();
		double[][] b = m2.matrixArray.clone();
		double[][] output;
		if( !(a.length == b.length && a[0].length == b[0].length) )
		{
			System.out.println("Could not add matrices.");
			output = null;
		}
		else
		{
			int numCols = a.length, numRows = a[0].length;
			output = new double[numCols][numRows];
			for(int col = 0; col < numCols; col++)
			{
				for(int row = 0; row < numRows; row++)
				{
					output[col][row] = a[col][row] + b[col][row]; 
				}
			}
		}
		return new Matrix(output);
	}
	
	public Matrix multiply(Matrix m2)  // complete
	{
		double[][] a = this.matrixArray.clone();
		double[][] b = m2.matrixArray.clone();
		double[][] output;
		if( !(a.length == b[0].length) )
		{
			System.out.println("Could not multiply matrices, number columns of A does not equal number rows of B.");
			output = null;
		}
		else
		{
			output = new double[b.length][a[0].length];
			for(int row = 0; row < a[0].length; row++)
			{
				for(int col = 0; col < b.length; col++)
				{
					for(int pos = 0; pos < a.length; pos++) 
					{
						output[col][row] += a[pos][row] * b[col][pos]; 
					}
				}
			}
			// If values that should be zero are non-zero because of error representing the numbers in binary, fix that.
			for(int col = 0; col < output.length; col++)
			{
				for(int row = 0; row < output[0].length; row++)
				{
					output[col][row] = correctBinary(output[col][row]);
				}
			}
		}
		return new Matrix(output);
	}
	
	public Matrix add(double scalar)  // complete
	{
		double[][] output = this.matrixArray.clone();
		for(int col = 0; col < output.length; col++)
		{
			for(int row = 0; row < output[0].length; row++)
			{
				output[col][row] = output[col][row] + scalar;
			}
		}
		return new Matrix(output);	}

	public Matrix multiply(double scalar)  // complete
	{
		double[][] output = this.matrixArray.clone();
		for(int col = 0; col < output.length; col++)
		{
			for(int row = 0; row < output[0].length; row++)
			{
				output[col][row] = output[col][row] * scalar;
			}
		}
		// If values that should be zero are non-zero because of error representing the numbers in binary, fix that.
		for(int col = 0; col < output.length; col++)
		{
			for(int row = 0; row < output[0].length; row++)
			{
				output[col][row] = correctBinary(output[col][row]);
			}
		}

		return new Matrix(output);
	}

	public static Matrix identity(int dimension)  // complete
	{
		double[][] matrixContents = new double[dimension][dimension];
		for(int col = 0, row = 0; col < dimension; col++, row++)
		{
			matrixContents[col][row] = 1;
		}
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
			minorMatrix = new Matrix(null);
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
	
	public Matrix minorMatrix(int skippedColumn, int skippedRow)  //  complete (assistant method to "inverse")
	{
		Matrix minorMatrix;
		double[][] majorArray = this.matrixArray.clone();
		int numCols = majorArray.length, numRows = majorArray[0].length;
		if( !(numCols >= 3 && numRows >= 3) )
		{
			System.out.println("Could not make minor matrix.");
			minorMatrix = new Matrix(null);
		}
		else
		{
			double[][] minorArray = new double[numCols-1][numRows-1];
			int colMinor = 0, rowMinor;
			for(int colMajor = 0; colMajor < numCols; colMajor++)
			{
				rowMinor = 0;
				for(int rowMajor = 0; rowMajor < numRows; rowMajor++)
				{
					// System.out.println("colMinor " + colMinor + " rowMinor " + rowMinor + ", colMajor " + colMajor + " rowMajor " + rowMajor);
					if( !(rowMajor == skippedRow || colMajor == skippedColumn) )
					{
						minorArray[colMinor][rowMinor] = majorArray[colMajor][rowMajor];
						// System.out.println("   minorArray[colMinor][rowMinor] " + minorArray[colMinor][rowMinor]) ;
						rowMinor++;
					}
				}
				if( !(colMajor == skippedColumn) )
				{
					colMinor++;
				}			}
			minorMatrix = new Matrix(minorArray);
		}

		return minorMatrix;
	}
	
	public Matrix transpose()  // complete
	{
		double[][] a = this.matrixArray.clone();
		int numColsA = a.length, numRowsA = a[0].length;
		double[][] transposeA = new double[numRowsA][numColsA];  // number of rows in transposed matrix = number of columns in original and vice versa.
		for(int col = 0; col < numColsA; col++)
			for(int row = 0; row < numRowsA; row++)
				transposeA[row][col] = a[col][row];
		return new Matrix(transposeA);
	}
		
	public Matrix inverse()  // complete
	{
		double[][] a = this.matrixArray.clone();
		int numColsA = a.length, numRowsA = a[0].length;
		double[][] inverse = new double[numColsA][numRowsA];
		Matrix m;
		
		if(!(numColsA == numRowsA) )
		{
			System.out.println("Cannot find inverse of non-square matrix.");
			m = new Matrix(null);
		}
		else if(this.determinant() == 0)
		{
			System.out.println("The matrix is singular (non-invertible), determinant 0.");
			m = new Matrix(null);
		}
		else if(numColsA == 2 && numRowsA == 2)
		{
			inverse[0][0] = a[1][1];
			inverse[1][1] = a[0][0];
			inverse[1][0] = -1*a[1][0];
			inverse[0][1] = -1*a[0][1];
			m = new Matrix(inverse);
			inverse = m.multiply(1/this.determinant() ).matrixArray;
			m = new Matrix(inverse);
		}
		else {
			// calculate matrix of determinants of minors
			for(int col = 0; col < numColsA; col++)
				for(int row = 0; row < numRowsA; row++)
				{
					inverse[col][row] = this.minorMatrix(col, row).determinant();
				}
			
			// apply "checkerboard" or minuses
			for(int col = 0; col < numColsA; col++)
				for(int row = 0; row < numRowsA; row++)
				{
					if(col%2 == 1 ^ row%2 == 1)  // if exactly one of the coordinates is odd
						inverse[col][row] = -1*inverse[col][row];
				}

			// transpose checkerboarded matrix of minors
			m = new Matrix(inverse);
			inverse = m.transpose().matrixArray;

			// multiply matrix by 1/determinant of A.
			m = new Matrix(inverse);
			inverse = m.multiply(1/this.determinant() ).matrixArray;

			m = new Matrix(inverse);
		}
		return m;
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
						output += correctBinary(a[column][row]) + " ";
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
	
	// Corrects binary misrepresentation rounding error for integer values.
	public static double correctBinary(double d)   // complete
	{
		double corrected;
		double tolerance = 1E-15;
		if(Math.abs(d)%1 > 1-tolerance || Math.abs(d)%1 < tolerance)
			corrected = Math.round(d);
		else
			corrected = d;
		return corrected;
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

