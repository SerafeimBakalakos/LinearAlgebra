namespace MGroup.LinearAlgebra.Tests.MatrixFormatConversions
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.MatrixFormatConversions;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;

	public class CscToCsrConverterTests : MatrixFormatConverterTestSuite<CscMatrix, CsrMatrix>
	{
		protected override void AssertEqualRawArrays(CsrMatrix expected, CsrMatrix computed, MatrixComparer comparer)
		{
			comparer.AssertEqual(expected.RawValues, computed.RawValues);
			comparer.AssertEqual(expected.RawRowOffsets, computed.RawRowOffsets);
			comparer.AssertEqual(expected.RawColIndices, computed.RawColIndices);
		}

		protected override IMatrixFormatConverter<CscMatrix, CsrMatrix> CreateConverter() => new CscToCsrConverter();

		protected override IIndexable2D CreateDefault(Pattern pattern, double multiplier)
		{
			if (pattern == Pattern.One)
			{
				double[,] values = ArrayUtilities.DoToAll(SparseRectangular10by5.Matrix, x => x * multiplier);
				return Matrix.CreateFromArray(values);
			}
			else
			{
				double[,] values = ArrayUtilities.DoToAll(SquareInvertible10by10.Matrix, x => x * multiplier);
				return Matrix.CreateFromArray(values);
			}
		}

		protected override CscMatrix CreateOriginal(Pattern pattern, double multiplier)
		{
			if (pattern == Pattern.One)
			{
				double[] values = ArrayUtilities.DoToAll(SparseRectangular10by5.CscValues, x => x * multiplier);
				return CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
					values, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
					true);
			}
			else
			{
				double[] values = ArrayUtilities.DoToAll(SquareInvertible10by10.CscValues, x => x * multiplier);
				return CscMatrix.CreateFromArrays(SquareInvertible10by10.Order, SquareInvertible10by10.Order,
					values, SquareInvertible10by10.CscRowIndices, SquareInvertible10by10.CscColOffsets,
					true);
			}
		}

		protected override CsrMatrix CreateResult(Pattern pattern, double multiplier)
		{
			if (pattern == Pattern.One)
			{
				double[] values = ArrayUtilities.DoToAll(SparseRectangular10by5.CsrValues, x => x * multiplier);
				return CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
					values, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
					true);
			}
			else
			{
				double[] values = ArrayUtilities.DoToAll(SquareInvertible10by10.CsrValues, x => x * multiplier);
				return CsrMatrix.CreateFromArrays(SquareInvertible10by10.Order, SquareInvertible10by10.Order,
					values, SquareInvertible10by10.CsrColIndices, SquareInvertible10by10.CsrRowOffsets,
					true);
			}
		}
	}
}
