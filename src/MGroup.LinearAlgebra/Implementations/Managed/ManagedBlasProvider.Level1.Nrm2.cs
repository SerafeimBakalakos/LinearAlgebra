#pragma warning disable SA1407 // Arithmetic expressions should declare precedence
namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public double Dnrm2(int n, double[] x, int offsetX, int incX)
		{
			AssertVector(x, n, offsetX, incX);

			if (n == 0)
			{
				return 0.0;
			}

			// NRM2 is guaranteed to avoid numerical overflow and underflow.
			// We must use a scaling strategy that enforces at every step:
			// scale = max(|x₁|, ..., |xᵢ|), sum = Σ(xⱼ / scale)²   for j ≤ i
			// This is done by updating the scale and sum, whenever a bigger value appears.
			// At the end we multiply the sqrt(sum) with scale to get rid of it.
			double scale = 0.0;
			double sum = 1.0;

			int ix = offsetX;
			for (int i = 0; i < n; i++)
			{
				double xi = x[ix];

				if (xi != 0.0)
				{
					double absXi = Math.Abs(xi);

					if (absXi > scale) 
					{
						double r = scale / absXi;
						sum = 1.0 + sum * r * r;
						scale = absXi;
					}
					else
					{
						double r = absXi / scale;
						sum += r * r;
					}
				}

				ix += incX;
			}

			return scale * Math.Sqrt(sum);
		}

		public double Dnrm2(double[] x)
		{
			// The algorithm is dependency-bound. No SIMD optimizations are possible. This overload provides consistency only.
			return Dnrm2(x.Length, x, 0, 1);
		}
	}
}
#pragma warning restore SA1407 // Arithmetic expressions should declare precedence
