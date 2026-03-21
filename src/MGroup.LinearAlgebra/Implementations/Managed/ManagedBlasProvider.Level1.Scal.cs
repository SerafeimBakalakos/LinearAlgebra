#pragma warning disable SA1407 // Arithmetic expressions should declare precedence
namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Dscal(int n, double alpha, double[] x, int offsetX, int incX)
		{
			AssertVector(x, n, offsetX, incX);
			if (incX == 1) // Contiguous access as a special case, to allow SIMD optimizations.
			{
				for (int i = 0; i < n; i++)
				{
					x[offsetX + i] *= alpha;
				}
			}
			else
			{
				int ix = offsetX;
				for (int i = 0; i < n; i++)
				{
					x[ix] *= alpha;
					ix += incX;
				}
			}
		}

		public void Dscal(double alpha, double[] x)
		{
			AssertNotNull(x);
			int n = x.Length;
			for (int i = 0; i < n; i++)
			{
				x[i] *= alpha;
			}
		}

		public void DscalTo(int n, double alpha, double[] x, int offsetX, int incX, double[] z, int offsetZ, int incZ)
		{
			AssertVector(x, n, offsetX, incX);
			AssertVector(z, n, offsetZ, incZ);
			if (incX == 1 && incZ == 1) // Contiguous access as a special case, to allow SIMD optimizations.
			{
				for (int i = 0; i < n; i++)
				{
					z[offsetZ + i] = alpha * x[offsetX + i];
				}
			}
			else
			{
				int ix = offsetX;
				int iz = offsetZ;
				for (int i = 0; i < n; i++)
				{
					z[iz] = alpha * x[ix];
					ix += incX;
					iz += incZ;
				}
			}
		}

		public void DscalTo(double alpha, double[] x, double[] z)
		{
			AssertSameLength(x, z);
			int n = x.Length;
			for (int i = 0; i < n; i++)
			{
				z[i] = alpha * x[i];
			}
		}
	}
}
#pragma warning restore SA1407 // Arithmetic expressions should declare precedence
