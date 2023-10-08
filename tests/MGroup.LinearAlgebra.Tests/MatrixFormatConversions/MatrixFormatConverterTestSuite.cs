namespace MGroup.LinearAlgebra.Tests.MatrixFormatConversions
{
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.MatrixFormatConversions;

	using MGroup.LinearAlgebra.Tests.TestData;
	using MGroup.LinearAlgebra.Tests.Utilities;

	using Xunit.Sdk;

	using Xunit;
	using System;

	public abstract class MatrixFormatConverterTestSuite<TMatrixIn, TMatrixOut> 
		where TMatrixIn : IMatrixView 
		where TMatrixOut : IMatrixView
	{
		public enum Pattern
		{
			One, Two
		}

		protected readonly MatrixComparer _comparer = new MatrixComparer(1E-15);

		#region test methods
		[Fact]
		public void TestConvertAndStorePattern()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public void TestConvertAndStorePattern1_ConvertOncePattern1()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertAndStorePattern1_ConvertOncePattern2()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertAndStorePattern1_ConvertOncePattern2_ConvertUsingPattern1()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted1;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted1 = converter.ConvertAndStorePatternConversion(original);
			CheckConvertedMatrix(converted1, pattern, multiplier);

			pattern = Pattern.Two;
			multiplier = -1.5;
			original = CreateOriginal(pattern, multiplier);
			TMatrixOut converted2 = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted2, pattern, multiplier);

			// ConvertOnce() should not have affected the stored pattern conversion
			pattern = Pattern.One;
			multiplier = 0.001;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted1);
			CheckConvertedMatrix(converted1, pattern, multiplier);
		}

		[Fact]
		public void TestConvertAndStorePattern1_ConvertUsingPattern1()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertAndStorePattern2()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertAndStorePattern2_ConvertUsingPattern1()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertAndStorePattern2_ConvertUsingPattern2()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertUsingPattern2()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertOncePattern1()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertOncePattern2()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertOncePattern2_ConvertAndStorePattern1()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertAndStorePattern1_ConvertUsingPattern1_ConvertOncePattern2_ConvertUsingPattern1()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted1;

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
			TMatrixOut converted2 = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted2, pattern, multiplier);

			pattern = Pattern.One;
			multiplier = -0.305;
			original = CreateOriginal(pattern, multiplier);
			converter.ConvertUsingStoredPatternConversion(original, converted1);
			CheckConvertedMatrix(converted1, pattern, multiplier);
		}

		[Fact]
		public void TestConvertAndStorePattern1_ConvertUsingPattern2()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertAndStorePattern1_MultipleConvertUsingPattern1()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertOnce()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);
			converted = converter.ConvertOnce(original);
			CheckConvertedMatrix(converted, pattern, multiplier);
		}

		[Fact]
		public void TestConvertOncePattern1_ConvertUsingPattern1()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted;

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
		public void TestConvertUsingPattern()
		{
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter = CreateConverter();
			Pattern pattern;
			double multiplier;
			TMatrixIn original;
			TMatrixOut converted = default;

			pattern = Pattern.One;
			multiplier = 1.0;
			original = CreateOriginal(pattern, multiplier);

			// Calling ConvertUsingStoredPatternConversion(), without calling ConvertAndStorePatternConversion() first must fail
			Assert.Throws<InvalidOperationException>(() => converter.ConvertUsingStoredPatternConversion(original, converted));
		}
		#endregion

		#region utilities methods
		private void CheckConvertedMatrix(TMatrixOut converted, Pattern pattern, double multiplier)
		{
			IIndexable2D defaultExpected = CreateDefault(pattern, multiplier);
			TMatrixOut resultExpected = CreateResult(pattern, multiplier);
			AssertEqualEntries(defaultExpected, converted);
			AssertEqualRawArrays(resultExpected, converted, _comparer);
		}

		private void AssertEqualEntries(IIndexable2D expected, TMatrixOut computed)
		{
			_comparer.AssertEqual(expected, computed);
		}

		private void AssertFailureOfConvertingWithDifferentStoredPattern(
			IMatrixFormatConverter<TMatrixIn, TMatrixOut> converter, TMatrixOut converted, Pattern pattern, double multiplier)
		{
			TMatrixIn original = CreateOriginal(pattern, multiplier);

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
		#endregion

		#region abstract methods
		protected abstract void AssertEqualRawArrays(TMatrixOut expected, TMatrixOut computed, MatrixComparer comparer);

		protected abstract IMatrixFormatConverter<TMatrixIn, TMatrixOut> CreateConverter();

		protected abstract IIndexable2D CreateDefault(Pattern pattern, double multiplier);

		protected abstract TMatrixIn CreateOriginal(Pattern pattern, double multiplier);

		protected abstract TMatrixOut CreateResult(Pattern pattern, double multiplier);
		#endregion
	}
}
