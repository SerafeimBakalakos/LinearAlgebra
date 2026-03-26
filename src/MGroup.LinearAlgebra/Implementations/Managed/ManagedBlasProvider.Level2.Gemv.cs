#pragma warning disable SA1407 // Arithmetic expressions should declare precedence
namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Dgemv(TransposeMatrix transA, int m, int n,
			double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX,
			double beta, double[] y, int offsetY, int incY)
		{
			if (transA == TransposeMatrix.NoTranspose)
			{
				DgemvNoTranspose(m, n, alpha, a, offsetA, ldA, x, offsetX, incX, beta, y, offsetY, incY);
			}
			else
			{
				DgemvColMajorTranspose(m, n, alpha, a, offsetA, ldA, x, offsetX, incX, beta, y, offsetY, incY);
			}
		}

		public void DgemvNoTranspose(int m, int n,
			double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX,
			double beta, double[] y, int offsetY, int incY)
		{
			AssertMatrixColMajor(a, m, n, offsetA, ldA);
			AssertVector(x, n, offsetX, incX);
			AssertVector(y, m, offsetY, incY);

			if (m == 0 || n == 0)
			{
				return;
			}

			if (alpha == 0.0)
			{
				OptimizationForAlpha0(m, beta, y, offsetY, incY);
				return;
			}

			if (incX == 1 && incY == 1) // Special fast path that is SIMD friendly
			{
				if (beta == 0.0)
				{
					Array.Clear(y, offsetY, m);
				}
				else if (beta != 1.0)
				{
					// Inline Dscal(m, beta, y, offsetY, incY);
					for (int i = 0; i < m; i++)
					{
						y[offsetY + i] *= beta;
					}
				}

				int colA = offsetA;
				for (int j = 0; j < n; j++)
				{
					double xj = x[offsetX + j];
					if (xj == 0.0)
					{
						colA += ldA;
						continue;
					}

					double alphaTimesX = alpha * xj;
					for (int i = 0; i < m; i++)
					{
						y[offsetY + i] += alphaTimesX * a[colA + i];
					}

					colA += ldA; // a[colA + i] is faster than a[offsetA + j * ldA + i]
				}
			}
			else
			{
				if (beta == 0.0)
				{
					// Inline Dclear(m, y, offsetY, incY);
					int iy = offsetY;
					for (int i = 0; i < m; i++)
					{
						y[iy] = 0.0;
						iy += incY;
					}
				}
				else if (beta != 1.0)
				{
					// Inline Dscal(m, beta, y, offsetY, incY);
					int iy = offsetY;
					for (int i = 0; i < m; i++)
					{
						y[iy] *= beta;
						iy += incY;
					}
				}

				int colA = offsetA;
				int indexX = offsetX;
				for (int j = 0; j < n; j++)
				{
					double xj = x[indexX];
					indexX += incX;
					if (xj == 0.0)
					{
						colA += ldA;
						continue;
					}

					double alphaTimesX = alpha * xj;
					int iy = offsetY;
					for (int i = 0; i < m; i++)
					{
						y[iy] += alphaTimesX * a[colA + i];
						iy += incY;
					}

					colA += ldA; // a[colA + i] is faster than a[offsetA + j * ldA + i]
				}
			}
		}

		public void DgemvColMajorNoTranspose(int m, int n, double[] a, double[] x, double[] y)
		{
			AssertMvmDimensions(m, n, a, x, y);
			Array.Clear(y, 0, m);
			
			int colA = 0;
			for (int j = 0; j < n; j++)
			{
				double xj = x[j];
				if (xj == 0.0)
				{
					colA += m;
					continue;
				}

				for (int i = 0; i < m; i++)
				{
					y[i] += xj * a[colA + i];
				}

				colA += m;
			}
		}

		public void DgemvColMajorTranspose(int m, int n,
			double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX,
			double beta, double[] y, int offsetY, int incY)
		{
			AssertMatrixColMajor(a, m, n, offsetA, ldA);
			AssertVector(x, m, offsetX, incX);
			AssertVector(y, n, offsetY, incY);

			if (m == 0 || n == 0)
			{
				return;
			}

			if (alpha == 0.0)
			{
				OptimizationForAlpha0(n, beta, y, offsetY, incY);
				return;
			}

			if (incX == 1 && incY == 1) // Special fast path that is SIMD friendly
			{
				int colA = offsetA;
				for (int i = 0; i < n; i++)
				{
					double sum = 0;
					for (int j = 0; j < m; j++)
					{
						sum += a[colA + j] * x[offsetX + j];
					}

					int iy = offsetY + i;
					if (beta == 0.0)
					{
						y[iy] = alpha * sum;
					}
					else if (beta == 1.0)
					{
						y[iy] += alpha * sum;
					}
					else
					{
						y[iy] = beta * y[iy] + alpha * sum;
					}

					colA += ldA; // a[colA + j] is faster than a[offsetA + i * ldA + j]
				}
			}
			else
			{
				int iy = offsetY;
				int colA = offsetA;
				for (int i = 0; i < n; i++)
				{
					double sum = 0;
					int jx = offsetX;
					for (int j = 0; j < m; j++)
					{
						// sum += A^T[i,j] * x[j] <=> sum += A[j,i] * x[j]
						sum += a[colA + j] * x[jx];
						jx += incX;
					}

					if (beta == 0.0)
					{
						y[iy] = alpha * sum;
					}
					else if (beta == 1.0)
					{
						y[iy] += alpha * sum;
					}
					else
					{
						y[iy] = beta * y[iy] + alpha * sum;
					}

					iy += incY;
					colA += ldA; // a[colA + j] is faster than a[offsetA + i * ldA + j]
				}
			}
		}

		public void DgemvColMajorTranspose(int m, int n, double[] a, double[] x, double[] y)
		{
			AssertMvmTransposeDimensions(m, n, a, x, y);

			int colA = 0;
			for (int i = 0; i < n; i++)
			{
				double sum = 0;
				for (int j = 0; j < m; j++)
				{
					sum += a[colA + j] * x[j];
				}

				y[i] = sum;
				colA += m;
			}
		}

		public void DgemvRowMajorNoTranspose(int m, int n,
			double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX,
			double beta, double[] y, int offsetY, int incY)
		{
			AssertMatrixRowMajor(a, m, n, offsetA, ldA);
			AssertVector(x, n, offsetX, incX);
			AssertVector(y, m, offsetY, incY);

			if (m == 0 || n == 0)
			{
				return;
			}

			if (alpha == 0.0)
			{
				OptimizationForAlpha0(m, beta, y, offsetY, incY);
				return;
			}

			if (incX == 1 && incY == 1) // Special fast path that is SIMD friendly
			{
				int rowA = offsetA;
				for (int i = 0; i < m; i++)
				{
					double sum = 0;
					for (int j = 0; j < n; j++)
					{
						sum += a[rowA + j] * x[offsetX + j];
					}

					int iy = offsetY + i;
					if (beta == 0.0)
					{
						y[iy] = alpha * sum;
					}
					else if (beta == 1.0)
					{
						y[iy] += alpha * sum;
					}
					else
					{
						y[iy] = beta * y[iy] + alpha * sum;
					}

					rowA += ldA; // a[rowA + j] is faster than a[offsetA + i * ldA + j]
				}
			}
			else
			{
				int iy = offsetY;
				int rowA = offsetA;
				for (int i = 0; i < m; i++)
				{
					double sum = 0;
					int jx = offsetX;
					for (int j = 0; j < n; j++)
					{
						sum += a[rowA + j] * x[jx];
						jx += incX;
					}

					if (beta == 0.0)
					{
						y[iy] = alpha * sum;
					}
					else if (beta == 1.0)
					{
						y[iy] += alpha * sum;
					}
					else
					{
						y[iy] = beta * y[iy] + alpha * sum;
					}

					iy += incY;
					rowA += ldA; // a[rowA + j] is faster than a[offsetA + i * ldA + j]
				}
			}
		}

		public void DgemvRowMajorNoTranspose(int m, int n, double[] a, double[] x, double[] y)
		{
			AssertMvmDimensions(m, n, a, x, y);

			int rowA = 0;
			for (int i = 0; i < m; i++)
			{
				double sum = 0;
				for (int j = 0; j < n; j++)
				{
					sum += a[rowA + j] * x[j];
				}

				y[i] = sum;
				rowA += n;
			}
		}

		public void DgemvRowMajorTranspose(int m, int n,
			double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX,
			double beta, double[] y, int offsetY, int incY)
		{
			AssertMatrixRowMajor(a, m, n, offsetA, ldA);
			AssertVector(x, m, offsetX, incX);
			AssertVector(y, n, offsetY, incY);

			if (m == 0 || n == 0)
			{
				return;
			}

			if (alpha == 0.0)
			{
				OptimizationForAlpha0(n, beta, y, offsetY, incY);
				return;
			}

			if (incX == 1 && incY == 1) // Special fast path that is SIMD friendly
			{
				if (beta == 0.0)
				{
					Array.Clear(y, offsetY, n);
				}
				else if (beta != 1.0)
				{
					// Inline Dscal(m, beta, y, offsetY, incY);
					for (int i = 0; i < n; i++)
					{
						y[offsetY + i] *= beta;
					}
				}

				int rowA = offsetA;
				for (int j = 0; j < m; j++)
				{
					double xj = x[offsetX + j];
					if (xj == 0.0)
					{
						rowA += ldA;
						continue;
					}

					double alphaTimesX = alpha * xj;
					for (int i = 0; i < n; i++)
					{
						y[offsetY + i] += alphaTimesX * a[rowA + i];
					}

					rowA += ldA;
				}
			}
			else
			{
				if (beta == 0.0)
				{
					// Inline Dclear(n, y, offsetY, incY);
					int iy = offsetY;
					for (int i = 0; i < n; i++)
					{
						y[iy] = 0.0;
						iy += incY;
					}
				}
				else if (beta != 1.0)
				{
					// Inline Dscal(n, beta, y, offsetY, incY);
					int iy = offsetY;
					for (int i = 0; i < n; i++)
					{
						y[iy] *= beta;
						iy += incY;
					}
				}

				int rowA = offsetA;
				int indexX = offsetX;
				for (int j = 0; j < m; j++)
				{
					double xj = x[indexX];
					indexX += incX;
					if (xj == 0.0)
					{
						rowA += ldA;
						continue;
					}

					double alphaTimesX = alpha * xj;
					int iy = offsetY;
					for (int i = 0; i < n; i++)
					{
						// y[i] += alpha * A^T[i, j] * x[j] => y[i] += alpha * A[j, i] * x[j]
						y[iy] += alphaTimesX * a[rowA + i];
						iy += incY;
					}

					rowA += ldA; // a[rowA + i] is faster than a[offsetA + j * ldA + i]
				}
			}
		}

		public void DgemvRowMajorTranspose(int m, int n, double[] a, double[] x, double[] y)
		{
			AssertMvmTransposeDimensions(m, n, a, x, y);
			Array.Clear(y, 0, n);

			int rowA = 0;
			for (int j = 0; j < m; j++)
			{
				double xj = x[j];
				if (xj == 0.0)
				{
					rowA += n;
					continue;
				}

				for (int i = 0; i < n; i++)
				{
					y[i] += xj * a[rowA + i];
				}

				rowA += n;
			}
		}

		private void OptimizationForAlpha0(int lengthY, double beta, double[] y, int offsetY, int incY)
		{
			if (beta == 0.0)
			{
				Dclear(lengthY, y, offsetY, incY);
			}
			else if (beta == 1.0) // Do nothing
			{
			}
			else
			{
				Dscal(lengthY, beta, y, offsetY, incY);
			}
		}
	}
}
#pragma warning restore SA1407 // Arithmetic expressions should declare precedence
