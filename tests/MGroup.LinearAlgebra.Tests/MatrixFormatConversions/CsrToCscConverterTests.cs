namespace MGroup.LinearAlgebra.Tests.MatrixFormatConversions
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.MatrixFormatConversions;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;

	public class CsrToCscConverterTests : MatrixFormatConverterTestSuite<CsrMatrix, CscMatrix>
	{
		protected override void AssertEqualRawArrays(CscMatrix expected, CscMatrix computed, MatrixComparer comparer)
		{
			comparer.AssertEqual(expected.RawValues, computed.RawValues);
			comparer.AssertEqual(expected.RawColOffsets, computed.RawColOffsets);
			comparer.AssertEqual(expected.RawRowIndices, computed.RawRowIndices);
		}

		protected override IMatrixFormatConverter<CsrMatrix, CscMatrix> CreateConverter() => new CsrToCscConverter();

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

		protected override CsrMatrix CreateOriginal(Pattern pattern, double multiplier)
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

		protected override CscMatrix CreateResult(Pattern pattern, double multiplier)
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
	}
}
