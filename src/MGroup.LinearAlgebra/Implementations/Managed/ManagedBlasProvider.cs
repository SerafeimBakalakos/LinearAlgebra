namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	using DotNumerics.LinearAlgebra.CSLapack;

	/// <summary>
	/// Provides managed C# implementations of the linear algebra operations defined by <see cref="IBlasProvider"/>.
	/// </summary>
	public partial class ManagedBlasProvider : IBlasProvider
	{
		//TODO: perhaps these should not be static.
		private static readonly DAXPY daxpy = new DAXPY();
		private static readonly DDOT ddot = new DDOT();
		private static readonly DGEMM dgemm = new DGEMM();
		private static readonly DGEMV dgemv = new DGEMV();
		private static readonly DNRM2 dnrm2 = new DNRM2();
		private static readonly DSCAL dscal = new DSCAL();
		private static readonly DTRSV dtrsv = new DTRSV();

		public static ManagedBlasProvider UniqueInstance { get; } = new ManagedBlasProvider();

		private ManagedBlasProvider() { } // private constructor for singleton pattern

		#region BLAS Level 2

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dgemv.aspx
		/// </summary>
		public void Dgemv(TransposeMatrix transA, int m, int n,
			double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX,
			double beta, double[] y, int offsetY, int incY)
			=> dgemv.Run(transA.Translate(), m, n, alpha, a, offsetA, ldA, x, offsetX, incX, beta, ref y, offsetY, incY);

		public void DgemvRowMajor(TransposeMatrix transA, int m, int n, double[] a, double[] x, double[] y)
		{
			if (transA == TransposeMatrix.NoTranspose)
			{
				for (var i = 0; i < m; ++i)
				{
					var rowStart = i * n;
					double sum = 0;
					for (var j = 0; j < n; ++j)
					{
						sum += a[rowStart + j] * x[j];
					}
					y[i] = sum;
				}
			}
			else
			{
				for (var j = 0; j < n; ++j)
				{
					double sum = 0;
					for (var i = 0; i < m; ++i)
					{
						sum += a[i * n + j] * x[i];
					}
					y[j] = sum;
				}
			}
		}

		public void Dspmv(StoredTriangle uplo, int n,
			double alpha, double[] a, int offsetA, double[] x, int offsetX, int incX,
			double beta, double[] y, int offsetY, int incY)
		{
			// y = alpha * L * x + beta * y 
			CblasLevel2Implementations.LowerTimesVectorPackedRowMajor(
				CblasLevel2Implementations.Diagonal.Regular, n, alpha, a, offsetA, x, offsetX, incX, beta, y, offsetY, incY);

			// y = alpha * U * x + y, where U has 0 diagonal
			CblasLevel2Implementations.UpperTimesVectorPackedColMajor(
				CblasLevel2Implementations.Diagonal.Zero, n, alpha, a, offsetA, x, offsetX, incX, 1.0, y, offsetY, incY);
		}

		public void Dtpmv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
			double[] a, int offsetA, double[] x, int offsetX, int incX)
		{
			// The copy may be avoidable in trangular operations, if we start the dot products from the bottom
			var input = new double[x.Length];
			Array.Copy(x, input, x.Length);

			CblasLevel2Implementations.Diagonal managedDiag = (diag == DiagonalValues.NonUnit) ?
				CblasLevel2Implementations.Diagonal.Regular : CblasLevel2Implementations.Diagonal.Unit;
			if (UseUpperImplementation(uplo, transA))
			{
				CblasLevel2Implementations.UpperTimesVectorPackedColMajor(
					managedDiag, n, 1.0, a, offsetA, input, offsetX, incX, 0.0, x, offsetX, incX);
			}
			else
			{
				CblasLevel2Implementations.LowerTimesVectorPackedRowMajor(
					managedDiag, n, 1.0, a, offsetA, input, offsetX, incX, 0.0, x, offsetX, incX);
			}
		}

		public void Dtpsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
			double[] a, int offsetA, double[] x, int offsetX, int incX)
		{
			bool unit = (diag == DiagonalValues.Unit) ? true : false;
			if (UseUpperImplementation(uplo, transA))
			{
				CblasLevel2Implementations.BackSubstitutionPackedColMajor(unit, n, a, offsetA, x, offsetX, incX);
			}
			else CblasLevel2Implementations.ForwardSubstitutionPackedRowMajor(unit, n, a, offsetA, x, offsetX, incX);
		}

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dtrsv.aspx
		/// </summary>
		public void Dtrsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
			double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
			=> dtrsv.Run(uplo.Translate(), transA.Translate(), diag.Translate(), n, a, offsetA, ldA, ref x, offsetX, incX);
		#endregion

		#region BLAS Level 3

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dgemm.aspx
		/// </summary>
		public void Dgemm(TransposeMatrix transA, TransposeMatrix transB, int m, int n, int k, double alpha,
			double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB, double beta, double[] c, int offsetC, int ldC)
			=> dgemm.Run(transA.Translate(), transB.Translate(), m, n, k, alpha, a, offsetA, ldA, b, offsetB, ldB,
				beta, ref c, offsetC, ldC);
		#endregion

		[Conditional("DEBUG")]
		private static void AssertNotNull(double[] x)
		{
			Debug.Assert(x != null, "Vector x is null");
		}

		[Conditional("DEBUG")]
		private static void AssertSameLength(double[] x, double[] y)
		{
			Debug.Assert(x != null, "Vector x is null");
			Debug.Assert(y != null, "Vector y is null");
			Debug.Assert(x.Length == y.Length, "Vectors have different length");
		}

		[Conditional("DEBUG")]
		private static void AssertSameLength(double[] x, double[] y, double[] z)
		{
			Debug.Assert(x != null, "Vector x is null");
			Debug.Assert(y != null, "Vector y is null");
			Debug.Assert(z != null, "Vector z is null");
			Debug.Assert(x.Length == y.Length, "Vectors x,y have different length");
			Debug.Assert(y.Length == z.Length, "Vectors y,z have different length");
		}

		[Conditional("DEBUG")]
		private static void AssertVector(double[] x, int n, int offset, int inc)
		{
			Debug.Assert(x != null, "Vector is null");
			Debug.Assert(offset >= 0, "offset >= 0");
			Debug.Assert(inc > 0, "inc > 0");
			Debug.Assert(n >= 0, "n >= 0");

			if (n > 0)
			{
				Debug.Assert(offset + (n - 1) * inc < x.Length, "Out of bounds");
			}
		}

		private static bool UseUpperImplementation(StoredTriangle uplo, TransposeMatrix transA)
		{
			//if (transA == TransposeMatrix.ConjugateTranspose)
			//    throw new ArgumentException("Cannot use conjugate transpose operations for double matrices and vectors.");

			if (uplo == StoredTriangle.Upper && transA == TransposeMatrix.NoTranspose) return true;
			if (uplo == StoredTriangle.Upper && transA == TransposeMatrix.Transpose) return false;
			if (uplo == StoredTriangle.Lower && transA == TransposeMatrix.NoTranspose) return true;
			if (uplo == StoredTriangle.Lower && transA == TransposeMatrix.Transpose) return false;
			throw new Exception("This code should not have been reached");
		}
	}
}
