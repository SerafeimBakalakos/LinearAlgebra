#pragma warning disable SA1407 // Arithmetic expressions should declare precedence
namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Daxpy(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
		{
			AssertVector(x, n, offsetX, incX);
			AssertVector(y, n, offsetY, incY);

			if (alpha == 0)
			{
				return;
			}

			if (incX == 1 && incY == 1) // Contiguous access as a special case, to allow SIMD optimizations.
			{
				for (int i = 0; i < n; i++)
				{
					y[offsetY + i] += alpha * x[offsetX + i];
				}
			}
			else
			{
				int ix = offsetX;
				int iy = offsetY;
				for (int i = 0; i < n; i++)
				{
					y[iy] += alpha * x[ix];
					ix += incX;
					iy += incY;
				}
			}
		}

		public void Daxpy(double alpha, double[] x, double[] y)
		{
			AssertSameLength(x, y);

			if (alpha == 0)
			{
				return;
			}

			int n = x.Length;
			for (int i = 0; i < n; i++)
			{
				y[i] += alpha * x[i];
			}
		}

		public void DaxpyTo(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY, double[] z, int offsetZ, int incZ)
		{
			AssertVector(x, n, offsetX, incX);
			AssertVector(y, n, offsetY, incY);
			AssertVector(z, n, offsetZ, incZ);

			if (alpha == 0)
			{
				Dcopy(n, y, offsetY, incY, z, offsetZ, incZ);
				return;
			}

			if (incX == 1 && incY == 1 && incZ == 1) // Contiguous access as a special case, to allow SIMD optimizations.
			{
				for (int i = 0; i < n; i++)
				{
					z[offsetZ + i] = alpha * x[offsetX + i] + y[offsetY + i];
				}
			}
			else
			{
				int ix = offsetX;
				int iy = offsetY;
				int iz = offsetZ;
				for (int i = 0; i < n; i++)
				{
					z[iz] = alpha * x[ix] + y[iy];
					ix += incX;
					iy += incY;
					iz += incZ;
				}
			}
		}

		public void DaxpyTo(double alpha, double[] x, double[] y, double[] z)
		{
			AssertSameLength(x, y, z);
			int n = x.Length;
			if (alpha == 0)
			{
				Array.Copy(y, z, n);
				return;
			}

			for (int i = 0; i < n; i++)
			{
				z[i] = alpha * x[i] + y[i];
			}
		}
	}
}
#pragma warning restore SA1407 // Arithmetic expressions should declare precedence
