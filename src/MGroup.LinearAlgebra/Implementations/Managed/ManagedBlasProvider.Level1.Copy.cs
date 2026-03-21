#pragma warning disable SA1407 // Arithmetic expressions should declare precedence
namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Dcopy(int n, double[] x, int offsetX, int incX, double[] z, int offsetZ, int incZ)
		{
			AssertVector(x, n, offsetX, incX);
			AssertVector(z, n, offsetZ, incZ);
			if (incX == 1 && incZ == 1)
			{
				Array.Copy(x, offsetX, z, offsetZ, n);
			}
			else
			{
				int ix = offsetX;
				int iz = offsetZ;
				for (int i = 0; i < n; i++)
				{
					z[iz] = x[ix];
					ix += incX;
					iz += incZ;
				}
			}
		}
	}
}
#pragma warning restore SA1407 // Arithmetic expressions should declare precedence
