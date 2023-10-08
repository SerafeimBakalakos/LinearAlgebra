namespace MGroup.LinearAlgebra.MatrixFormatConversions
{
	using System;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Text;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Matrices;

	public class CsrToCscConverter : IMatrixFormatConverter<CsrMatrix, CscMatrix>
	{
		/// <summary>
		/// Maps between the values arrays, such that: cscValues[<see cref="_mapCsrToCsc"/>[k]] = csrValues[k]
		/// </summary>
		private int[] _mapCsrToCsc;

		public CscMatrix ConvertAndStorePatternConversion(CsrMatrix original)
		{
			// This method allocates and stores an extra int[nnz], where nnz is the number of non-zero entries of the
			// original matrix.

			int m = original.NumRows;
			int n = original.NumColumns;
			int nnz = original.RawValues.Length;
			var cscValues = new double[nnz];
			var cscRowIndices = new int[nnz];
			var cscColOffsets = new int[n + 1];

			_mapCsrToCsc = new int[nnz];

			CsrToCscAndPattern(m, n, original.RawRowOffsets, original.RawColIndices, original.RawValues,
				cscColOffsets, cscRowIndices, cscValues, _mapCsrToCsc);
			return CscMatrix.CreateFromArrays(m, n, cscValues, cscRowIndices, cscColOffsets, false);
		}

		public CscMatrix ConvertOnce(CsrMatrix original)
		{
			int m = original.NumRows;
			int n = original.NumColumns;
			int nnz = original.RawValues.Length;
			var cscValues = new double[nnz];
			var cscRowIndices = new int[nnz];
			var cscColOffsets = new int[n + 1];

			CsrToCsc(m, n, original.RawRowOffsets, original.RawColIndices, original.RawValues,
				cscColOffsets, cscRowIndices, cscValues);

			return CscMatrix.CreateFromArrays(m, n, cscValues, cscRowIndices, cscColOffsets, false);
		}

		public void ConvertUsingStoredPatternConversion(CsrMatrix original, CscMatrix result) 
		{
			if (_mapCsrToCsc == null || result == null)
			{
				throw new InvalidOperationException("Must call the method ConvertAndSavePatternConversion(...), " +
					"before calling ConvertUsingStoredPatternConversion(...)");
			}

			double[] csrValues = original.RawValues;
			double[] cscValues = result.RawValues;
			
			if (_mapCsrToCsc.Length != csrValues.Length || _mapCsrToCsc.Length != cscValues.Length) 
			{
				// Does not ensure that the patterns are identical, but at least it is better than throwing crytpic OutOfRangeException.
				throw new InvalidOperationException("The method ConvertAndSavePatternConversion(...) was called previously" +
					"for a matrix with different sparsity pattern than this one.");
			}

			for (int k = 0; k < csrValues.Length; k++) 
			{
				cscValues[_mapCsrToCsc[k]] = csrValues[k];
			}
		}
		

		/// <summary>
		/// Compute B = A for CSR matrix A, CSC matrix B.
		/// Also, with the appropriate arguments can also be used to:
		///   - compute B = A^T for CSR matrix A, CSR matrix B
		///   - compute B = A^T for CSC matrix A, CSC matrix B
		///   - convert CSC->CSR
		/// Complexity: Linear. Specifically O(nnz(A) + max(numRows,numCols))
		/// </summary>
		/// <param name="numRows">Number of rows in A.</param>
		/// <param name="numCols">Number of columns in A.</param>
		/// <param name="Ap">Row pointers. Size = numRows+1.</param>
		/// <param name="Aj">Column indices. Size = nnz(A). They are not assumed to be in sorted order.</param>
		/// <param name="Ax">Non-zero values. Size = nnz(A).</param>
		/// <param name="Bp">Preallocated ouput argument. Column pointers. Size = numCols+1</param>
		/// <param name="Bi">Preallocated ouput argument. Row indices. Size = nnz(A). They will be in sorted order.</param>
		/// <param name="Bx">Preallocated ouput argument. Non-zero values. Size = nnz(A).</param>
		internal static void CsrToCsc(int numRows, int numCols, int[] Ap, int[] Aj, double[] Ax, int[] Bp, int[] Bi, double[] Bx)
		{
			int nnz = Ap[numRows];
			FindCscColumnOffsets(numCols, nnz, Aj, Bp);

			// Find the arrays of row indices Bi[] and values Bx[]
			for (int row = 0; row < numRows; row++)
			{
				for (int k = Ap[row]; k < Ap[row + 1]; k++)
				{
					int col = Aj[k];

					// First time a specific "col" is met, "dest" is the offset of its start.
					// Corrupt "Bp[col]" by incrementing it, so the next time "col" is met,
					// "dest" will be the offset of the next entry of this column
					int dest = Bp[col];
					Bp[col]++;

					Bi[dest] = row;
					Bx[dest] = Ax[k];
				}
			}

			RestoreCscColumnOffsets(numCols, Bp);
		}

		/// <summary>
		/// Also, with the appropriate arguments can also be used to convert CSC->CSR
		/// </summary>
		/// <param name="numRows">Number of rows in A.</param>
		/// <param name="numCols">Number of columns in A.</param>
		/// <param name="Ap">Row pointers. Size = numRows+1.</param>
		/// <param name="Aj">Column indices. Size = nnz(A). They are not assumed to be in sorted order.</param>
		/// <param name="Ax">Non-zero values. Size = nnz(A).</param>
		/// <param name="Bp">Preallocated ouput argument. Column pointers. Size = numCols+1</param>
		/// <param name="Bi">Preallocated ouput argument. Row indices. Size = nnz(A). They will be in sorted order.</param>
		/// <param name="Bx">Preallocated ouput argument. Non-zero values. Size = nnz(A).</param>
		/// <param name="mapCsrToCsc">
		/// Maps between the values arrays, such that: 
		/// <paramref name="Bx"/>[<paramref name="mapCsrToCsc"/>[k]] = <paramref name="Ax"/>[k]
		/// </param>
		internal static void CsrToCscAndPattern(int numRows, int numCols, int[] Ap, int[] Aj, double[] Ax, int[] Bp, int[] Bi, 
			double[] Bx, int[] mapCsrToCsc)
		{
			int nnz = Ap[numRows];
			FindCscColumnOffsets(numCols, nnz, Aj, Bp);

			// Find the arrays of row indices Bi[] and values Bx[]
			for (int row = 0; row < numRows; row++)
			{
				for (int k = Ap[row]; k < Ap[row + 1]; k++)
				{
					int col = Aj[k];

					// First time a specific "col" is met, "dest" is the offset of its start.
					// Corrupt "Bp[col]" by incrementing it, so the next time "col" is met,
					// "dest" will be the offset of the next entry of this column
					int dest = Bp[col];
					Bp[col]++;

					Bi[dest] = row;
					Bx[dest] = Ax[k];
					mapCsrToCsc[k] = dest;
				}
			}

			RestoreCscColumnOffsets(numCols, Bp);
		}

		private static void FindCscColumnOffsets(int numCols, int nnz, int[] Aj, int[] Bp)
		{
			// Count the non-zero entries per column of A 
			for (int k = 0; k < nnz; k++)
			{
				Bp[Aj[k]]++;
			}

			// Sum the nnz per column to get the correct column offsets Bp[]
			for (int col = 0, sum = 0; col < numCols; col++)
			{
				int temp = Bp[col];
				Bp[col] = sum;
				sum += temp;
			}
			Bp[numCols] = nnz;
		}

		private static void RestoreCscColumnOffsets(int numCols, int[] Bp)
		{
			// At this point each entry Bp[col] has been incremented once for each entry of the corresponding column col.
			// Thus Bp[col] = offset of 1st entry of column col+1. Restore it to its correct state.
			for (int col = 0, last = 0; col <= numCols; col++)
			{
				int temp = Bp[col];
				Bp[col] = last;
				last = temp;
			}
		}
	}
}
