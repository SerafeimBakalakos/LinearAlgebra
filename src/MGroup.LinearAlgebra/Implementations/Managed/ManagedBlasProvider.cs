#pragma warning disable SA1407 // Arithmetic expressions should declare precedence
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
		private static void AssertMatrixColMajor(double[] a, int m, int n, int offset, int leadDim)
		{
			Debug.Assert(a != null, "Matrix cannot be null");
			Debug.Assert(m >= 0, "Required: m >= 0");
			Debug.Assert(n >= 0, "Required: n >= 0");
			Debug.Assert(leadDim >= 1, "Required: leadDim > 1");
			Debug.Assert(leadDim >= m, "Required: leadDim >= m");
			Debug.Assert(offset >= 0, "Required: offset >= 0");
			Debug.Assert(offset < a.Length, "Required: offset < a.Length");

			if (m > 0 && n > 0)
			{
				int lastIndex = offset + (n - 1) * leadDim + (m - 1);
				Debug.Assert(lastIndex < a.Length, "Not enough space in matrix array");
			}
		}

		[Conditional("DEBUG")]
		private static void AssertMatrixRowMajor(double[] a, int m, int n, int offset, int leadDim)
		{
			Debug.Assert(a != null, "Matrix cannot be null");
			Debug.Assert(m >= 0, "Required: m >= 0");
			Debug.Assert(n >= 0, "Required: n >= 0");
			Debug.Assert(leadDim >= 1, "Required: leadDim > 1");
			Debug.Assert(leadDim >= n, "Required: leadDim >= n");
			Debug.Assert(offset >= 0, "Required: offset >= 0");
			Debug.Assert(offset < a.Length, "Required: offset < a.Length");

			if (m > 0 && n > 0)
			{
				int lastIndex = offset + (m - 1) * leadDim + (n - 1);
				Debug.Assert(lastIndex < a.Length, "Not enough space in matrix array");
			}
		}

		[Conditional("DEBUG")]
		private static void AssertMatrixTriangular(double[] a, int n, int offset, int ldA)
		{
			Debug.Assert(a != null, "Matrix cannot be null");
			Debug.Assert(n >= 0, "Required: n >= 0");
			Debug.Assert(ldA >= Math.Max(1, n), "Required: ldA >= n");
			Debug.Assert(offset >= 0, "Required: offset >= 0");
			Debug.Assert(offset < a.Length, "Required: offset < a.Length");

			if (n > 0)
			{
				// Last used element is bottom of last column
				int lastIndex = offset + (n - 1) * ldA + (n - 1);
				Debug.Assert(lastIndex < a.Length, "Not enough space in matrix array");
			}
		}

		[Conditional("DEBUG")]
		private static void AssertMvmDimensions(int m, int n, double[] a, double[] x, double[] y)
		{
			Debug.Assert(a != null, "Matrix A cannot be null");
			Debug.Assert(x != null, "LHS vector cannot be null");
			Debug.Assert(y != null, "RHS vector cannot be null");
			Debug.Assert(m >= 1, "Required: m >= 1");
			Debug.Assert(n >= 1, "Required: n >= 1");
			Debug.Assert(a.Length == m * n, "Not enough space in matrix array");
			Debug.Assert(x.Length == n, "Not enough space in LHS vector array");
			Debug.Assert(y.Length == m, "Not enough space in RHS vector array");
		}

		[Conditional("DEBUG")]
		private static void AssertMvmTransposeDimensions(int m, int n, double[] a, double[] x, double[] y)
		{
			Debug.Assert(a != null, "Matrix A cannot be null");
			Debug.Assert(x != null, "LHS vector cannot be null");
			Debug.Assert(y != null, "RHS vector cannot be null");
			Debug.Assert(m >= 1, "Required: m >= 1");
			Debug.Assert(n >= 1, "Required: n >= 1");
			Debug.Assert(a.Length == m * n, "Not enough space in matrix array");
			Debug.Assert(x.Length == m, "Not enough space in LHS vector array");
			Debug.Assert(y.Length == n, "Not enough space in RHS vector array");
		}

		[Conditional("DEBUG")]
		private static void AssertVector(double[] x, int n, int offset, int inc)
		{
			Debug.Assert(x != null, "Vector cannot be null");
			Debug.Assert(n >= 0, "Required: n >= 0");
			Debug.Assert(inc > 0, "Required: inc > 0");
			Debug.Assert(offset >= 0, "Required: offset >= 0");
			Debug.Assert(offset < x.Length, "Required: offset x.Length");

			if (n > 0)
			{
				int lastIndex = offset + (n - 1) * inc;
				Debug.Assert(lastIndex < x.Length, "Not enough space in vector array");
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
#pragma warning restore SA1407 // Arithmetic expressions should declare precedence
