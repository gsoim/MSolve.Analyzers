﻿namespace MGroup.Analyzers.Interfaces
{
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.MSolve.Discretization.Interfaces;

	public interface INonLinearProvider : IAnalyzerProvider
	{
		double CalculateRhsNorm(IVectorView rhs);

		//TODO: Very generic name. There is also a similar method in IImplictIntegrationProvider.
		void ProcessInternalRhs(ISubdomain subdomain, IVectorView solution, IVector rhs);
	}
}
