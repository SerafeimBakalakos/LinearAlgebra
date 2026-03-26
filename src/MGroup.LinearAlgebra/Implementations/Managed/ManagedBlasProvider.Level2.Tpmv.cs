#pragma warning disable SA1407 // Arithmetic expressions should declare precedence
namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
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

	}
}
#pragma warning restore SA1407 // Arithmetic expressions should declare precedence
