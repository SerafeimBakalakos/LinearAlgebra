namespace MGroup.LinearAlgebra.Iterative.Preconditioning.IncompleteFactorizations
{
	using System;
	using System.Collections.Generic;
	using System.Runtime.CompilerServices;
	using System.Text;

	using MGroup.LinearAlgebra.Exceptions;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	// Adapted from https://docs.octave.org/doxygen/stable/dd/ddc/____ichol_____8cc_source.html
	public class IncompleteCholesky : IPreconditioner
	{
		private void Ichol0(CscMatrix matrix, bool michol)
		{
			// Input matrix internal arrays
			int n = matrix.NumColumns;
			int[] cidx = matrix.RawColOffsets;
			int[] ridx = matrix.RawRowIndices;
			double[] data = matrix.RawValues;

			// Allocate work arrays
			//TODO: These should be allocated once per IChol0 object and then cleared or rest to -1
			int[] Lfirst = new int[n];
			int[] Llist = new int[n];
			int[] iw = new int[n];
			double[] dropsums = new double[n];

			// Initialize work arrays (dropsums is initialized to 0 during allocation)
			for (int i = 0; i < n; i++)
			{
				Lfirst[i] = -1;
				Llist[i] = -1;
				iw[i] = -1;
			}

			// Loop over all columns
			for (int k = 0; k < n; k++)
			{
				int j1 = cidx[k];
				int j2 = cidx[k + 1];
				for (int j = j1; j < j2; j++)
				{
					iw[ridx[j]] = j;
				}

				int jrow = Llist[k];

				// Iterate over each nonzero element in the actual row.
				while (jrow != -1)
				{
					int jjrow = Lfirst[jrow];
					int jend = cidx[jrow + 1];
					for (int jj = jjrow; jj < jend; jj++)
					{
						int r = ridx[jj];
						int jw = iw[r];
						double tl = data[jj] * data[jjrow];
						if (jw != -1)
						{
							data[jw] -= tl;
						}
						else
						{
							// Because of the symmetry of the matrix, we know
							// the drops in the column r are also in the column k.
							if (michol == true)
							{
								dropsums[r] -= tl;
								dropsums[k] -= tl;
							}
						}
					}

					// Update the linked list and the first entry of the actual column.
					if ((jjrow + 1) < jend)
					{
						Lfirst[jrow]++;
						int j = jrow;
						jrow = Llist[jrow];
						Llist[j] = Llist[ridx[Lfirst[j]]];
						Llist[ridx[Lfirst[j]]] = j;
					}
					else
					{
						jrow = Llist[jrow];
					}
				}

				if (michol == true)
				{
					data[j1] += dropsums[k];
				}

				// Test for j1 == j2 must be first to avoid invalid ridx[j1] access
				if ((j1 == j2) || (ridx[j1] != k))
				{
					throw new InvalidPivotException("ichol: encountered a pivot equal to 0");
				}

				if (!ichol_checkpivot(data[j1]))
				{
					break;
				}

				data[cidx[k]] = Math.Sqrt(data[j1]);

				// Update Llist and Lfirst with the k-column information.
				// Also, scale the column elements by the pivot and reset the working array iw.
				if (k < (n - 1))
				{
					iw[ridx[j1]] = -1;
					for (int i = j1 + 1; i < j2; i++)
					{
						iw[ridx[i]] = -1;
						data[i] /= data[j1];
					}

					Lfirst[k] = j1;
					if ((Lfirst[k] + 1) < j2)
					{
						Lfirst[k]++;
						int jjrow = ridx[Lfirst[k]];
						Llist[k] = Llist[jjrow];
						Llist[jjrow] = k;
					}
				}
			}
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private static void CheckPivotReal(double pivot)
		{
			if (pivot < 0)
			{
				throw new InvalidPivotException("ichol: negative pivot encountered");
			}
		}

		public IPreconditioner CopyWithInitialSettings() => throw new NotImplementedException();
		public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector) => throw new NotImplementedException();
		public void UpdateMatrix(IMatrixView matrix, bool isPatternModified) => throw new NotImplementedException();
	}
}
