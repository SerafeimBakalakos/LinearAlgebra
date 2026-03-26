#pragma warning disable SA1407 // Arithmetic expressions should declare precedence
namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Dclear(int n, double[] x, int offsetX, int incX)
		{
			AssertVector(x, n, offsetX, incX);
			if (incX == 1)
			{
				Array.Clear(x, offsetX, n);
			}
			else
			{
				int ix = offsetX;
				for (int i = 0; i < n; i++)
				{
					x[ix] = 0.0;
					ix += incX;
				}
			}
		}
	}
}
#pragma warning restore SA1407 // Arithmetic expressions should declare precedence
