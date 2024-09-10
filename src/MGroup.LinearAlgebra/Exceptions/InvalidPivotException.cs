namespace MGroup.LinearAlgebra.Exceptions
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	/// <summary>
	/// The exception that is thrown when the value of a pivot element is illegal for the specific linear algebra algorithm.
	/// That usually means the pivot is zero or negative, although that depends on the algorithm itself. This exception does not
	/// necessarily mean that the matrix is singular or indefinite, although this is the most common reason for invalid pivot
	/// values. Instead, it means that a specific algorithm (e.g. incomplete holesky) will breakdown for a specific matrix.
	/// </summary>
	public class InvalidPivotException : Exception
	{
		/// <summary>
		/// Initializes a new instance of the <see cref="InvalidPivotException"/> class.
		/// </summary>
		public InvalidPivotException()
		{ }

		/// <summary>
		/// Initializes a new instance of the <see cref="InvalidPivotException"/> class with a specified error message.
		/// </summary>
		/// <param name="message">The error message that explains the reason for the exception.</param>
		public InvalidPivotException(string message)
			: base(message)
		{ }

		/// <summary>
		/// Initializes a new instance of the <see cref="InvalidPivotException"/> class with a specified error message 
		/// and a reference to the inner exception that is the cause of this exception.
		/// </summary>
		/// <param name="message">The error message that explains the reason for the exception.</param>
		/// <param name="inner">
		/// The exception that is the cause of the current exception. If the innerException parameter is not a null reference,
		/// the current exception is raised in a catch block that handles the inner exception.
		/// </param>
		public InvalidPivotException(string message, Exception inner)
			: base(message, inner)
		{ }
	}
}
