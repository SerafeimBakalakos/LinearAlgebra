namespace MGroup.LinearAlgebra.MatrixFormatConversions
{
	using System;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Text;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Matrices;

	public class CsrToCscConverter
	{
		/// <summary>
		/// Maps between the values arrays, such that:  cscValues[<see cref="_mapCsrToCsc"/>[k]] = csrValues[k]
		/// </summary>
		private int[] _mapCsrToCsc;


		/// <summary>
		/// Converts the matrix <paramref name="original"/>, without reading or overwritting the stored sparsity pattern 
		/// conversion. Use this when you only need to convert between the two matrix storage formats once. Otherwise, it is 
		/// faster to use <see cref="ConvertAndSavePatternConversion(CsrMatrix)"/> once and 
		/// <see cref="ConvertUsingStoredPatternConversion(CsrMatrix, CscMatrix)"/> multiple times.
		/// </summary>
		/// <param name="original"></param>
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

		/// <summary>
		/// Converts the matrix <paramref name="original"/> and also stores the conversion between the two formats for this exact
		/// sparsity pattern. Therefore subsequent conversions can be done much faster by calling 
		/// <see cref="ConvertUsingStoredPatternConversion(CsrMatrix, CscMatrix)"/>, provided those input and output matrices 
		/// have the exact same sparsity pattern as the one stored by the last call to 
		/// <see cref="ConvertAndSavePatternConversion(CsrMatrix)"/>. This method allocates and stores an extra int[nnz], where 
		/// nnz is the number of non-zero entries in <paramref name="original"/>. If the conversion needs to be done only once,
		/// use <see cref="ConvertOnce(CsrMatrix)"/> instead, since that method does not allocate any extra memory and does less 
		/// work.
		/// </summary>
		/// <param name="original"></param>
		public CscMatrix ConvertAndSavePatternConversion(CsrMatrix original)
		{
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

		/// <summary>
		/// Converts the matrix <paramref name="original"/> and the writes into <paramref name="result"/>, overwriting
		/// its non-zero values array. Both <paramref name="original"/> and <paramref name="result"/> must match the same 
		/// sparsity pattern, created and saved by the last call to <see cref="ConvertAndSavePatternConversion(CsrMatrix)"/>. 
		/// Use this method together with <see cref="ConvertUsingStoredPatternConversion(CsrMatrix, CscMatrix)"/>, when the input
		/// and output matrix have the same pattern over multiple conversions between the 2 storage formats. Otherwise 
		/// <see cref="ConvertOnce(CsrMatrix)"/> is faster.
		/// </summary>
		/// <param name="original">
		/// The original matrix to convert. Must have the same sparsity pattern as the input matrix of the last call to
		/// <see cref="ConvertAndSavePatternConversion(CsrMatrix)"/>.
		/// </param>
		/// <param name="result">
		/// It's non-zero values will be overwritten with the corresponding entries of <paramref name="original"/>, but not its 
		/// pattern arrays. It must have the same sparsity pattern as the returned matrix of the last call to 
		/// <see cref="ConvertAndSavePatternConversion(CsrMatrix)"/>. It does not need to be cleared first.
		/// </param>
		public void ConvertUsingStoredPatternConversion(CsrMatrix original, CscMatrix result) 
		{
			double[] csrValues = original.RawValues;
			double[] cscValues = result.RawValues;
			Debug.Assert(cscValues.Length == csrValues.Length);
			Debug.Assert(_mapCsrToCsc.Length == csrValues.Length);

			for (int k = 0; k < csrValues.Length; k++) 
			{
				cscValues[_mapCsrToCsc[k]] = csrValues[k];
			}
		}
		

		/// <summary>
		/// Compute B = A for CSR matrix A, CSC matrix B.
		/// Also, with the appropriate arguments can also be used to:
		///   - compute B = A ^ t for CSR matrix A, CSR matrix B
		///   - compute B = A ^ t for CSC matrix A, CSC matrix B
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
