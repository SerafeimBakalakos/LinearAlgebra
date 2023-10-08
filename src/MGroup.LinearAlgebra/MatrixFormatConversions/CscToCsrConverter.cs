namespace MGroup.LinearAlgebra.MatrixFormatConversions
{
	using System;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Text;
	using MGroup.LinearAlgebra.Matrices;

	public class CscToCsrConverter : IMatrixFormatConverter<CscMatrix, CsrMatrix>
	{
		/// <summary>
		/// Maps between the values arrays, such that: cscValues[<see cref="_mapCscToCsr"/>[k]] = csrValues[k]
		/// </summary>
		private int[] _mapCscToCsr;

		public CsrMatrix ConvertAndStorePatternConversion(CscMatrix original)
		{
			// This method allocates and stores an extra int[nnz], where nnz is the number of non-zero entries of the
			// original matrix.

			int m = original.NumRows;
			int n = original.NumColumns;
			int nnz = original.RawValues.Length;
			var csrValues = new double[nnz];
			var csrColIndices = new int[nnz];
			var csrRowOffsets = new int[m + 1];

			_mapCscToCsr = new int[nnz];

			CsrToCscConverter.CsrToCscAndPattern(n, m, original.RawColOffsets, original.RawRowIndices, original.RawValues,
				csrRowOffsets, csrColIndices, csrValues, _mapCscToCsr);
			return CsrMatrix.CreateFromArrays(m, n, csrValues, csrColIndices, csrRowOffsets, false);
		}

		public CsrMatrix ConvertOnce(CscMatrix original)
		{
			int m = original.NumRows;
			int n = original.NumColumns;
			int nnz = original.RawValues.Length;
			var csrValues = new double[nnz];
			var csrColIndices = new int[nnz];
			var csrRowOffsets = new int[m + 1];

			CsrToCscConverter.CsrToCsc(n, m, original.RawColOffsets, original.RawRowIndices, original.RawValues,
				csrRowOffsets, csrColIndices, csrValues);

			return CsrMatrix.CreateFromArrays(m, n, csrValues, csrColIndices, csrRowOffsets, false);
		}

		public void ConvertUsingStoredPatternConversion(CscMatrix original, CsrMatrix result)
		{
			if (_mapCscToCsr == null || result == null)
			{
				throw new InvalidOperationException("Must call the method ConvertAndSavePatternConversion(...), " +
					"before calling ConvertUsingStoredPatternConversion(...)");
			}

			double[] cscValues = original.RawValues;
			double[] csrValues = result.RawValues;

			if (_mapCscToCsr.Length != cscValues.Length || _mapCscToCsr.Length != csrValues.Length)
			{
				// Does not ensure that the patterns are identical, but at least it is better than throwing crytpic OutOfRangeException.
				throw new InvalidOperationException("The method ConvertAndSavePatternConversion(...) was called previously" +
					"for a matrix with different sparsity pattern than this one.");
			}

			for (int k = 0; k < cscValues.Length; k++)
			{
				csrValues[_mapCscToCsr[k]] = cscValues[k];
			}
		}
	}
}
