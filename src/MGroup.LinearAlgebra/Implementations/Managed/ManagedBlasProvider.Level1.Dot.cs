#pragma warning disable SA1407 // Arithmetic expressions should declare precedence
namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
		{
			AssertVector(x, n, offsetX, incX);
			AssertVector(y, n, offsetY, incY);
			double sum = 0;
			if (incX == 1 && incY == 1) // Contiguous access as a special case, to allow SIMD optimizations.
			{
				for (int i = 0; i < n; i++)
				{
					sum += x[offsetX + i] * y[offsetY + i];
				}
			}
			else
			{
				int ix = offsetX;
				int iy = offsetY;
				for (int i = 0; i < n; i++)
				{
					sum += x[ix] * y[iy];
					ix += incX;
					iy += incY;
				}
			}

			return sum;
		}

		public double Ddot(double[] x, double[] y)
		{
			AssertSameLength(x, y);
			double sum = 0;
			int n = x.Length;
			for (int i = 0; i < n; i++)
			{
				sum += x[i] * y[i];
			}

			return sum;
		}
	}
}
#pragma warning restore SA1407 // Arithmetic expressions should declare precedence
