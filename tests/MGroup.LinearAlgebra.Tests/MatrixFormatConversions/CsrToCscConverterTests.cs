namespace MGroup.LinearAlgebra.Tests.MatrixFormatConversions
{
	using System;
	using System.Collections.Generic;
	using System.ComponentModel.DataAnnotations;
	using System.Text;
	using CSparse;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.MatrixFormatConversions;
	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;

	using Xunit;
	using Xunit.Sdk;

	public static class CsrToCscConverterTests
	{
		private enum Pattern
		{
			One, Two
		}

		private static readonly MatrixComparer comparer = new MatrixComparer(1E-15);

		// [DONE] ConvertAndStore() pattern 1
		// [DONE] ConvertAndStore() pattern 1, then ConvertOnce() pattern 1
		// [DONE] ConvertAndStore() pattern 1, then ConvertOnce() pattern 2
		// [DONE] ConvertAndStore() pattern 1, then ConvertOnce() pattern 2, then ConvertUsingPattern pattern 1

		// [DONE] ConvertAndStore() pattern 1, then ConvertUsingPattern() pattern 1
		// [DONE] ConvertAndStore() pattern 1, then multiple ConvertUsingPattern() pattern 1

		// [DONE] ConvertOnce() pattern 1
		// [DONE] ConvertOnce() pattern 1, then ConvertUsingPattern() pattern 1 (without ConvertAndStore(): converter must throw exception here.

		// [DONE] ConvertUsingPattern() pattern 1 (without ConvertAndStore()): converter must throw exception here.


		// [DONE] ConvertAndStore() pattern 1, then ConvertUsingPattern() pattern 1, then ConvertOnce() pattern 1
		// [DONE] ConvertAndStore() pattern 1, then ConvertUsingPattern() pattern 1, then ConvertOnce() pattern 2
		// [DONE] ConvertAndStore() pattern 1, then ConvertUsingPattern() pattern 1, then ConvertOnce() pattern 2, then ConvertUsingPattern() pattern 1
		// [DONE] ConvertAndStore() pattern 1, then ConvertUsingPattern() pattern 1, then ConvertOnce() pattern 2, then ConvertAndStore() pattern 1



		// [DONE] ConvertAndStore() pattern 1, then ConvertUsingPattern() pattern 2 -> fail
		// [DONE] ConvertAndStore() pattern 1, then ConvertUsingPattern() pattern 1, then ConvertUsingPattern() pattern 2 -> fail
		// [DONE] ConvertAndStore() pattern 1, then ConvertUsingPattern() pattern 1, then ConvertAndStore() pattern 2
		// ConvertAndStore() pattern 1, then ConvertUsingPattern() pattern 1, then ConvertAndStore() pattern 2, then ConvertUsingPattern() pattern 1 -> fail
		// [DONE] ConvertAndStore() pattern 1, then ConvertUsingPattern() pattern 1, then ConvertAndStore() pattern 2, then ConvertUsingPattern() pattern 2 

		[Fact]
		public static void TestConvertAndStorePattern()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertOncePattern1()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			// ConvertOnce() should work after calls to ConvertAndStorePatternConversion()
			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertOncePattern2()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			// ConvertOnce() should work after calls to ConvertAndStorePatternConversion() even with different pattern
			pattern = Pattern.Two;
			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertOncePattern2_ConvertUsingPattern1()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted1;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted1 = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted1, pattern, multiplier);

			pattern = Pattern.Two;
			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			CscMatrix converted2 = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted2, pattern, multiplier);

			// ConvertOnce() should not have affected the stored pattern conversion
			pattern = Pattern.One;
			multiplier = 0.001;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted1);
			CheckConvertedMatrix(converted1, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertUsingPattern1()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertAndStorePattern2()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted);
			CheckConvertedMatrix(converted, pattern, multiplier);

			pattern = Pattern.Two;
			multiplier = 0.001;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertAndStorePattern2_ConvertUsingPattern1()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted);
			CheckConvertedMatrix(converted, pattern, multiplier);

			pattern = Pattern.Two;
			multiplier = 0.001;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			pattern = Pattern.One;
			multiplier = -0.305;
			AssertFailureOfConvertingWithDifferentStoredPattern(converter, converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertAndStorePattern2_ConvertUsingPattern2()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted);
			CheckConvertedMatrix(converted, pattern, multiplier);

			pattern = Pattern.Two;
			multiplier = 0.001;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			pattern = Pattern.Two;
			multiplier = -0.305;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertUsingPattern2()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted);
			CheckConvertedMatrix(converted, pattern, multiplier);

			pattern = Pattern.Two;
			multiplier = 0.001;
			AssertFailureOfConvertingWithDifferentStoredPattern(converter, converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertOncePattern1()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted);
			CheckConvertedMatrix(converted, pattern, multiplier);

			multiplier = 0.001;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertOncePattern2()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted);
			CheckConvertedMatrix(converted, pattern, multiplier);

			pattern = Pattern.Two;
			multiplier = 0.001;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertOncePattern2_ConvertAndStorePattern1()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted);
			CheckConvertedMatrix(converted, pattern, multiplier);

			pattern = Pattern.Two;
			multiplier = 0.001;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			pattern = Pattern.One;
			multiplier = -0.305;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertOncePattern2_ConvertUsingPattern1()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted1;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted1 = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted1, pattern, multiplier);

			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted1);
			CheckConvertedMatrix(converted1, pattern, multiplier);

			pattern = Pattern.Two;
			multiplier = 0.001;
			original = CreateOriginal(pattern, multiplier);
			CscMatrix converted2 = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted2, pattern, multiplier);

			pattern = Pattern.One;
			multiplier = -0.305;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted1);
			CheckConvertedMatrix(converted1, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_ConvertUsingPattern2()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			pattern = Pattern.Two;
			multiplier = -1.5;
			AssertFailureOfConvertingWithDifferentStoredPattern(converter, converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertAndStorePattern1_MultipleConvertUsingPattern1()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			for (int t = 0; t < 10; ++t)
			{
				multiplier = 0.01 * (t + 1);
				original = CreateOriginal(pattern, multiplier);
				converter.ConvertUsingStoredPatternConversion(original, converted);
				CheckConvertedMatrix(converted, pattern, multiplier);
			}
		}

		[Fact]
		public static void TestConvertOnce()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public static void TestConvertOncePattern1_ConvertUsingPattern1()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted, pattern, multiplier);

			// Calling ConvertOnce(), but not ConvertAndStorePatternConversion(), before ConvertUsingStoredPatternConversion()
			// must fail
			Assert.Throws<InvalidOperationException>(() => converter.ConvertUsingStoredPatternConversion(original, converted));
		}

		[Fact]
		public static void TestConvertUsingPattern()
		{
			var converter = new CsrToCscConverter();
			Pattern pattern;
			double multiplier;
			CsrMatrix original;
			CscMatrix converted = null;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);

			// Calling ConvertUsingStoredPatternConversion(), without calling ConvertAndStorePatternConversion() first must fail
			Assert.Throws<InvalidOperationException>(() => converter.ConvertUsingStoredPatternConversion(original, converted));
		}


		#region utilities
		private static void CheckConvertedMatrix(CscMatrix converted, Pattern pattern, double multiplier)
		{
			Matrix defaultExpected = CreateDefault(pattern, multiplier);
			CscMatrix resultExpected = CreateResult(pattern, multiplier);
			AssertEqualEntries(defaultExpected, converted);
			AssertEqualRawArrays(resultExpected, converted);
		}

		private static void AssertFailureOfConvertingWithDifferentStoredPattern(
			IMatrixFormatConverter<CsrMatrix, CscMatrix> converter, CscMatrix converted, Pattern pattern, double multiplier)
		{
			CsrMatrix original = CreateOriginal(pattern, multiplier);

			// If the sparsity pattern is different, an InvalidOperationException() may be thrown, but that is not guaranteed.
			try
			{
				converter.ConvertUsingStoredPatternConversion(original, converted);
			}
			catch (InvalidOperationException ex)
			{
				return;
			}

			// If no such exception is thrown, the converted matrix must be incorrect.
			Assert.Throws<TrueException>(() => CheckConvertedMatrix(converted, pattern, multiplier));
		}

		private static Matrix CreateDefault(Pattern pattern, double multiplier)
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

		private static CsrMatrix CreateOriginal(Pattern pattern, double multiplier)
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

		private static CscMatrix CreateResult(Pattern pattern, double multiplier)
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

		//private static CsrMatrix DoToNonZeros(CsrMatrix original, Func<double, double> func)
		//{
		//	double[] modifiedValues = ArrayUtilities.DoToAll(original.RawValues, func);
		//	return CsrMatrix.CreateFromArrays(original.NumRows, original.NumColumns, modifiedValues, 
		//		original.RawColIndices.Copy(), original.RawRowOffsets.Copy(), true);
		//}

		//private static CscMatrix DoToNonZeros(CscMatrix original, Func<double, double> func)
		//{
		//	double[] modifiedValues = ArrayUtilities.DoToAll(original.RawValues, func);
		//	return CscMatrix.CreateFromArrays(original.NumRows, original.NumColumns, modifiedValues,
		//		original.RawRowIndices.Copy(), original.RawColOffsets.Copy(), true);
		//}


		private static void AssertEqualEntries(Matrix expected, CscMatrix computed)
		{
			comparer.AssertEqual(expected, computed);
		}

		private static void AssertEqualRawArrays(CscMatrix expected, CscMatrix computed)
		{
			comparer.AssertEqual(expected.RawValues, computed.RawValues);
			comparer.AssertEqual(expected.RawColOffsets, computed.RawColOffsets);
			comparer.AssertEqual(expected.RawRowIndices, computed.RawRowIndices);
		}
		#endregion
	}
}
