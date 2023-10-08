using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Tests.Utilities
{
	public static class ArrayUtilities
	{
		public static bool AreEqual(int[] expected, int[] computed)
		{
			if (expected.Length != computed.Length)
			{
				return false;
			}
			for (int i = 0; i < expected.Length; ++i)
			{
				if (expected[i] != computed[i])
				{
					return false;
				}
			}
			return true;
		}

		public static bool AreEqual(int[,] expected, int[,] computed)
		{
			if ((expected.GetLength(0) != computed.GetLength(0)) || (expected.GetLength(1) != computed.GetLength(1)))
			{
				return false;
			}
			for (int i = 0; i < expected.GetLength(0); ++i)
			{
				for (int j = 0; j < expected.GetLength(1); ++j)
				{
					if (expected[i, j] != computed[i, j])
					{
						return false;
					}
				}
			}
			return true;
		}

		public static double[] DoToAll(double[] array, Func<double, double> func)
		{
			int n = array.Length;
			var result = new double[n];
			for (int i = 0; i < n; ++i) 
			{
				result[i] = func(array[i]);
			}
			return result;
		}

		public static double[,] DoToAll(double[,] array, Func<double, double> func)
		{
			int m = array.GetLength(0);
			int n = array.GetLength(1);
			var result = new double[m, n];
			for (int i = 0; i < m; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					result[i, j] = func(array[i, j]);
				}
			}
			return result;
		}
	}
}
