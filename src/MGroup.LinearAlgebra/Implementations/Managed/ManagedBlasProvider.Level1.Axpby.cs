#pragma warning disable SA1407 // Arithmetic expressions should declare precedence
namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
		{
			AssertVector(x, n, offsetX, incX);
			AssertVector(y, n, offsetY, incY);

			if (alpha == 0)
			{
				Dscal(n, beta, y, offsetY, incY);
				return;
			}

			if (beta == 0)
			{
				DscalTo(n, alpha, x, offsetX, incX, y, offsetY, incY);
				return;
			}

			if (incX == 1 && incY == 1) // Contiguous access as a special case, to allow SIMD optimizations.
			{
				for (int i = 0; i < n; i++)
				{
					int iy = offsetY + i;
					y[iy] = alpha * x[offsetX + i] + beta * y[iy];
				}
			}
			else
			{
				int ix = offsetX;
				int iy = offsetY;
				for (int i = 0; i < n; i++)
				{
					y[iy] = alpha * x[ix] + beta * y[iy];
					ix += incX;
					iy += incY;
				}
			}
		}

		public void Daxpby(double alpha, double[] x, double beta, double[] y)
		{
			AssertSameLength(x, y);

			if (alpha == 0)
			{
				Dscal(beta, y);
				return;
			}

			if (beta == 0)
			{
				DscalTo(alpha, x, y);
			}

			int n = x.Length;
			for (int i = 0; i < n; i++)
			{
				y[i] = alpha * x[i] + beta * y[i];
			}
		}

		public void DaxpbyTo(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY, double[] z, int offsetZ, int incZ)
		{
			AssertVector(x, n, offsetX, incX);
			AssertVector(y, n, offsetY, incY);
			AssertVector(z, n, offsetZ, incZ);

			if (alpha == 0)
			{
				DscalTo(n, beta, y, offsetY, incY, z, offsetZ, incZ);
				return;
			}

			if (beta == 0)
			{
				DscalTo(n, alpha, x, offsetX, incX, z, offsetZ, incZ);
				return;
			}

			if (incX == 1 && incY == 1 && incZ == 1) // Contiguous access as a special case, to allow SIMD optimizations.
			{
				for (int i = 0; i < n; i++)
				{
					z[offsetZ + i] = alpha * x[offsetX + i] + beta * y[offsetY + i];
				}
			}
			else
			{
				int ix = offsetX;
				int iy = offsetY;
				int iz = offsetZ;
				for (int i = 0; i < n; i++)
				{
					z[iz] = alpha * x[ix] + beta * y[iy];
					ix += incX;
					iy += incY;
					iz += incZ;
				}
			}
		}

		public void DaxpbyTo(double alpha, double[] x, double beta, double[] y, double[] z)
		{
			AssertSameLength(x, y, z);

			if (alpha == 0)
			{
				DscalTo(beta, y, z);
				return;
			}

			if (beta == 0)
			{
				DscalTo(alpha, x, z);
				return;
			}

			int n = x.Length;
			for (int i = 0; i < n; i++)
			{
				z[i] = alpha * x[i] + beta * y[i];
			}
		}
	}
}
#pragma warning restore SA1407 // Arithmetic expressions should declare precedence
