namespace MGroup.LinearAlgebra.MatrixFormatConversions
{
	using System;
	using System.Collections.Generic;
	using System.Text;
	using MGroup.LinearAlgebra.Matrices;

	public interface IMatrixFormatConverter<TMatrixIn, TMatrixOut>
	{
		/// <summary>
		/// Converts the matrix <paramref name="original"/> and also stores the conversion between the two formats for this exact
		/// sparsity pattern. Therefore subsequent conversions can be done much faster by calling 
		/// <see cref="ConvertUsingStoredPatternConversion(TMatrixIn, TMatrixOut)"/>, provided those input and output matrices 
		/// have the exact same sparsity pattern as the one stored by the last call to 
		/// <see cref="ConvertAndStorePatternConversion(TMatrixIn)"/>. This method may allocate and store extra memory to 
		/// faciliate subsequent conversions. If the conversion needs to be done only once, use
		/// <see cref="ConvertOnce(TMatrixIn)"/> instead, since that method does not allocate any extra memory and does less work.
		/// </summary>
		/// <param name="original">The matrix in the original storage format before conversion.</param>
		TMatrixOut ConvertAndStorePatternConversion(TMatrixIn original);

		/// <summary>
		/// Converts the matrix <paramref name="original"/>, without reading or overwritting the stored sparsity pattern 
		/// conversion. Use this when you only need to convert between the two matrix storage formats once. Otherwise, it is 
		/// faster to use <see cref="ConvertAndStorePatternConversion(TMatrixIn)"/> once and 
		/// <see cref="ConvertUsingStoredPatternConversion(TMatrixIn, TMatrixOut)"/> multiple times.
		/// </summary>
		/// <param name="original">The matrix in the original storage format before conversion.</param>
		TMatrixOut ConvertOnce(TMatrixIn original);

		/// <summary>
		/// Converts the matrix <paramref name="original"/> and the writes into <paramref name="result"/>, overwriting
		/// its non-zero values array. Both <paramref name="original"/> and <paramref name="result"/> must have the same 
		/// sparsity pattern, created and saved by the last call to <see cref="ConvertAndStorePatternConversion(TMatrixIn)"/>. 
		/// Use this method together with <see cref="ConvertUsingStoredPatternConversion(TMatrixIn, TMatrixOut)"/>, when the 
		/// input and output matrix have the same pattern over multiple conversions between the 2 storage formats. Otherwise 
		/// <see cref="ConvertOnce(TMatrixIn)"/> is faster.
		/// </summary>
		/// <param name="original">
		/// The original matrix to convert. Must have the same sparsity pattern as the input matrix of the last call to
		/// <see cref="ConvertAndStorePatternConversion(TMatrixIn)"/>.
		/// </param>
		/// <param name="result">
		/// It's non-zero values will be overwritten with the corresponding entries of <paramref name="original"/>, but not its 
		/// pattern arrays. It must have the same sparsity pattern as the returned matrix of the last call to 
		/// <see cref="ConvertAndStorePatternConversion(TMatrixIn)"/>. It does not need to be cleared first.
		/// </param>
		/// <exception cref="InvalidOperationException">
		/// Thrown if this method is called without calling <see cref="ConvertAndStorePatternConversion(TMatrixIn)"/> first.
		/// It may also be thrown if this method is called for an input matrix with different sparsity pattern than the sparsity
		/// pattern of the input matrix of the previous call to <see cref="ConvertAndStorePatternConversion(TMatrixIn)"/>. 
		/// However, there is no guarantee that this method can always check this case, since it is inefficient/impossible, 
		/// unless the two sparsity patterns have a different number on non-zero entries. 
		/// </exception>
		void ConvertUsingStoredPatternConversion(TMatrixIn original, TMatrixOut result);
	}
}
