using System;

namespace MGroup.LinearAlgebra.Implementations.Managed
{
	/// <summary>
	/// Provides managed C# implementations of the linear algebra operations defined by <see cref="IBlasProvider"/>. Uses loop unrolling by a factor of 4. Slight performance benefit for large vectors and matrices.
	/// </summary>
	public class ManagedUnrolledBlasProvider : IBlasProvider
	{
		public static ManagedUnrolledBlasProvider UniqueInstance { get; } = new ManagedUnrolledBlasProvider();

		private ManagedUnrolledBlasProvider() { } // private constructor for singleton pattern

		#region BLAS Level 1
		public void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
			=> throw new NotImplementedException();

		/// <summary>
		/// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/daxpy.aspx
		/// </summary>
		public void Daxpy(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
			=> throw new NotImplementedException();

		public double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
			=> throw new NotImplementedException();

		public double Dnrm2(int n, double[] x, int offsetX, int incX)
			=> throw new NotImplementedException();

		public void Dscal(int n, double alpha, double[] x, int offsetX, int incX)
		{
			if (incX == 1) // Contiguous access as a special case, to allow SIMD optimizations.
			{
				int i = 0;
				int limit = n - n % 4; // unroll by 4
				for (; i < limit; i += 4)
				{
					x[offsetX + i] *= alpha;
					x[offsetX + i + 1] *= alpha;
					x[offsetX + i + 2] *= alpha;
					x[offsetX + i + 3] *= alpha;
				}

				for (; i < n; i++) // remainder entries
				{
					x[offsetX + i] *= alpha;
				}
			}
			else
			{
				for (int i = 0; i < n; i++)
				{
					x[offsetX + i * incX] *= alpha;
				}
			}
		}
		#endregion

		#region BLAS Level 2

		public void Dgemv(TransposeMatrix transA, int m, int n, double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
			=> throw new NotImplementedException();

		public void DgemvRowMajor(TransposeMatrix transA, int m, int n, double[] a, double[] x, double[] y)
			=> throw new NotImplementedException();

		public void Dspmv(StoredTriangle uplo, int n, double alpha, double[] a, int offsetA, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
			=> throw new NotImplementedException();

		public void Dtpmv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n, double[] a, int offsetA, double[] x, int offsetX, int incX)
			=> throw new NotImplementedException();

		public void Dtpsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n, double[] a, int offsetA, double[] x, int offsetX, int incX)
			=> throw new NotImplementedException();

		public void Dtrsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
			=> throw new NotImplementedException();
		#endregion

		#region BLAS Level 3

		public void Dgemm(TransposeMatrix transA, TransposeMatrix transB, int m, int n, int k, double alpha, double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB, double beta, double[] c, int offsetC, int ldC)
			=> throw new NotImplementedException();
		#endregion
	}
}
