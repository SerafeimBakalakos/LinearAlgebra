#pragma warning disable SA1407 // Arithmetic expressions should declare precedence
namespace MGroup.LinearAlgebra.Implementations.Managed
{
	using System;
	using System.Diagnostics;

	public partial class ManagedBlasProvider : IBlasProvider
	{
		public void Dtrmv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
		{
			if (uplo == StoredTriangle.Lower)
			{
				if (transA == TransposeMatrix.NoTranspose)
				{
					DtrmvColMajorLowerNoTranspose(diag, n, a, offsetA, ldA, x, offsetX, incX);
				}
				else
				{
					DtrmvColMajorLowerTranspose(diag, n, a, offsetA, ldA, x, offsetX, incX);
				}
			}
			else
			{
				if (transA == TransposeMatrix.NoTranspose)
				{
					DtrmvColMajorUpperNoTranspose(diag, n, a, offsetA, ldA, x, offsetX, incX);
				}
				else
				{
					DtrmvColMajorUpperTranspose(diag, n, a, offsetA, ldA, x, offsetX, incX);
				}
			}
		}

		public void DtrmvRowMajor(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
		{
			if (uplo == StoredTriangle.Lower)
			{
				if (transA == TransposeMatrix.NoTranspose)
				{
					DtrmvColMajorUpperTranspose(diag, n, a, offsetA, ldA, x, offsetX, incX);
				}
				else
				{
					DtrmvColMajorUpperNoTranspose(diag, n, a, offsetA, ldA, x, offsetX, incX);
				}
			}
			else
			{
				if (transA == TransposeMatrix.NoTranspose)
				{
					DtrmvColMajorLowerTranspose(diag, n, a, offsetA, ldA, x, offsetX, incX);
				}
				else
				{
					DtrmvColMajorLowerNoTranspose(diag, n, a, offsetA, ldA, x, offsetX, incX);
				}
			}
		}

		public void DtrmvColMajorLowerNoTranspose(DiagonalValues diag, int n, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
		{
			AssertMatrixTriangular(a, n, offsetA, ldA);
			AssertVector(x, n, offsetX, incX);
			bool unitDiag = diag == DiagonalValues.Unit;

			if (n == 0)
			{
				return;
			}

			if (incX == 1) // Fast path
			{
				int colA = offsetA;

				for (int j = 0; j < n; j++)
				{
					double xj = x[offsetX + j];

					if (xj != 0.0)
					{
						double temp = xj;

						// Diagonal
						if (!unitDiag)
						{
							temp *= a[colA + j];
						}

						// Below diagonal: i = j+1 → n-1
						for (int i = j + 1; i < n; i++)
						{
							x[offsetX + i] += temp * a[colA + i];
						}

						x[offsetX + j] = temp;
					}

					colA += ldA;
				}
			}
			else
			{
				int jx = offsetX;
				int colA = offsetA;

				for (int j = 0; j < n; j++)
				{
					double xj = x[jx];

					if (xj != 0.0)
					{
						double temp = xj;

						if (!unitDiag)
						{
							temp *= a[colA + j];
						}

						int ix = jx + incX;

						for (int i = j + 1; i < n; i++)
						{
							x[ix] += temp * a[colA + i];
							ix += incX;
						}

						x[jx] = temp;
					}

					jx += incX;
					colA += ldA;
				}
			}
		}

		public void DtrmvColMajorLowerTranspose(DiagonalValues diag, int n, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
		{
			AssertMatrixTriangular(a, n, offsetA, ldA);
			AssertVector(x, n, offsetX, incX);
			bool unitDiag = diag == DiagonalValues.Unit;

			if (n == 0) return;

			if (incX == 1)
			{
				for (int i = n - 1; i >= 0; i--)
				{
					int colA = offsetA + i * ldA;

					double sum = x[offsetX + i];

					if (!unitDiag)
					{
						sum *= a[colA + i];
					}

					for (int j = i + 1; j < n; j++)
					{
						sum += a[offsetA + j * ldA + i] * x[offsetX + j];
					}

					x[offsetX + i] = sum;
				}
			}
			else
			{
				int ix = offsetX + (n - 1) * incX;

				for (int i = n - 1; i >= 0; i--)
				{
					double sum = x[ix];

					if (!unitDiag)
					{
						sum *= a[offsetA + i * ldA + i];
					}

					int jx = ix + incX;
					for (int j = i + 1; j < n; j++)
					{
						sum += a[offsetA + j * ldA + i] * x[jx];
						jx += incX;
					}

					x[ix] = sum;
					ix -= incX;
				}
			}
		}

		public void DtrmvColMajorUpperNoTranspose(DiagonalValues diag, int n, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
		{
			AssertMatrixTriangular(a, n, offsetA, ldA);
			AssertVector(x, n, offsetX, incX);
			bool unitDiag = diag == DiagonalValues.Unit;

			if (n == 0) return;

			if (incX == 1)
			{
				int colA = offsetA + (n - 1) * ldA;

				for (int j = n - 1; j >= 0; j--)
				{
					double xj = x[offsetX + j];

					if (xj != 0.0)
					{
						double temp = xj;

						if (!unitDiag)
						{
							temp *= a[colA + j];
						}

						for (int i = 0; i < j; i++)
						{
							x[offsetX + i] += temp * a[colA + i];
						}

						x[offsetX + j] = temp;
					}

					colA -= ldA;
				}
			}
			else
			{
				int jx = offsetX + (n - 1) * incX;
				int colA = offsetA + (n - 1) * ldA;

				for (int j = n - 1; j >= 0; j--)
				{
					double xj = x[jx];

					if (xj != 0.0)
					{
						double temp = xj;

						if (!unitDiag)
						{
							temp *= a[colA + j];
						}

						int ix = offsetX;
						for (int i = 0; i < j; i++)
						{
							x[ix] += temp * a[colA + i];
							ix += incX;
						}

						x[jx] = temp;
					}

					jx -= incX;
					colA -= ldA;
				}
			}
		}

		public void DtrmvColMajorUpperTranspose(DiagonalValues diag, int n, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
		{
			AssertMatrixTriangular(a, n, offsetA, ldA);
			AssertVector(x, n, offsetX, incX);
			bool unitDiag = diag == DiagonalValues.Unit;

			if (n == 0) return;

			if (incX == 1)
			{
				for (int i = 0; i < n; i++)
				{
					int colA = offsetA + i * ldA;

					double sum = x[offsetX + i];

					if (!unitDiag)
					{
						sum *= a[colA + i];
					}

					for (int j = 0; j < i; j++)
					{
						sum += a[offsetA + j * ldA + i] * x[offsetX + j];
					}

					x[offsetX + i] = sum;
				}
			}
			else
			{
				int ix = offsetX;

				for (int i = 0; i < n; i++)
				{
					double sum = x[ix];

					if (!unitDiag)
					{
						sum *= a[offsetA + i * ldA + i];
					}

					int jx = offsetX;
					for (int j = 0; j < i; j++)
					{
						sum += a[offsetA + j * ldA + i] * x[jx];
						jx += incX;
					}

					x[ix] = sum;
					ix += incX;
				}
			}
		}
	}
}
#pragma warning restore SA1407 // Arithmetic expressions should declare precedence
