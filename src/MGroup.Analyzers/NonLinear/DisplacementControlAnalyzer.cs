using System;
using System.Collections.Generic;
using System.Diagnostics;
using MGroup.Analyzers.Interfaces;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization.Interfaces;
using MGroup.MSolve.Logging;
using MGroup.Solvers;
using MGroup.Solvers.LinearSystems;

namespace MGroup.Analyzers.NonLinear
{
	public class DisplacementControlAnalyzer : NonLinearAnalyzerBase
	{
		/// <summary>
		/// This class solves the linearized geometrically nonlinear system of equations according to displacement control incremental-iterative method.
		/// This only works if there are no nodal loads or any loading condition other than prescribed displacements.
		/// </summary>
		/// <param name="model">Instance of the model that will be solved</param>
		/// <param name="solver">Instance of the solver that will solve the linear system of equations</param>
		/// <param name="provider">Instance of the problem type to be solved</param>
		/// <param name="subdomainUpdaters">Instance that updates constraints, right-hand-side vector, updates and resets state</param>
		/// <param name="numIncrements">Number of total load increments</param>
		/// <param name="maxIterationsPerIncrement">Number of maximum iterations within a load increment</param>
		/// <param name="numIterationsForMatrixRebuild">Number of iterations for the rebuild of the siffness matrix within a load increment</param>
		/// <param name="residualTolerance">Tolerance for the convergence criterion of the residual forces</param>
		private DisplacementControlAnalyzer(IModel model, ISolver solver, INonLinearProvider provider,
			IReadOnlyDictionary<int, INonLinearSubdomainUpdater> subdomainUpdaters,
			int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance) :
			base(model, solver, provider, subdomainUpdaters, numIncrements, maxIterationsPerIncrement,
				numIterationsForMatrixRebuild, residualTolerance)
		{ }

		/// <summary>
		/// Solves the nonlinear equations and calculates the displacements vector
		/// </summary>
		public override void Solve()
		{
			InitializeLogs();

			DateTime start = DateTime.Now;
			UpdateInternalVectors();
			for (int increment = 0; increment < numIncrements; increment++)
			{
				double errorNorm = 0;
				ClearIncrementalSolutionVector();
				UpdateRhs(increment);
				ScaleSubdomainConstraints(increment);

				double firstError = 0;
				int iteration = 0;
				for (iteration = 0; iteration < maxIterationsPerIncrement; iteration++)
				{
					AddEquivalentNodalLoadsToRHS(increment, iteration);
					solver.Solve();

					Dictionary<int, IVector> internalRhsVectors = CalculateInternalRhs(increment, iteration);
					errorNorm = UpdateResidualForcesAndNorm(increment, internalRhsVectors); // This also sets the rhs vectors in linear systems.
																							//Console.WriteLine($"Increment {increment}, iteration {iteration}: norm2(error) = {errorNorm}");

					if (iteration == 0) firstError = errorNorm;

					if (TotalDisplacementsPerIterationLog != null) TotalDisplacementsPerIterationLog.StoreDisplacements(uPlusdu);

					if (errorNorm < residualTolerance)
					{
						foreach (var subdomainLogPair in IncrementalLogs)
						{
							int subdomainID = subdomainLogPair.Key;
							TotalLoadsDisplacementsPerIncrementLog log = subdomainLogPair.Value;
							log.LogTotalDataForIncrement(increment, iteration, errorNorm,
								uPlusdu[subdomainID], internalRhsVectors[subdomainID]);
						}
						break;
					}

					SplitResidualForcesToSubdomains();
					if ((iteration + 1) % numIterationsForMatrixRebuild == 0)
					{
						provider.Reset();
						BuildMatrices();
					}
				}
				Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", iteration, firstError, errorNorm);
				SaveMaterialStateAndUpdateSolution();
			}

			// TODO: Logging should be done at each iteration. And it should be done using pull observers
			DateTime end = DateTime.Now;
			StoreLogResults(start, end);
		}

		/// <summary>
		/// Initializes internal vectors
		/// </summary>
		protected override void InitializeInternalVectors()
		{
			base.InitializeInternalVectors();
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				subdomainUpdaters[linearSystem.Subdomain.ID].ScaleConstraints(1 / (double)numIncrements);
			}
		}

		/// <summary>
		/// Adds equivalent nodal loads created by the prescribed DOFs to the right-hand-side vector.
		/// </summary>
		private void AddEquivalentNodalLoadsToRHS(int currentIncrement, int iteration)
		{
			if (iteration != 0)
				return;

			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;

				double scalingFactor = 1; //((double)currentIncrement + 2) / (currentIncrement + 1); //2; //
				IVector equivalentNodalLoads = provider.DirichletLoadsAssembler.GetEquivalentNodalLoads(linearSystem.Subdomain,
					u[id], scalingFactor);
				linearSystem.RhsVector.SubtractIntoThis(equivalentNodalLoads);

				model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
			}
		}

		/// <summary>
		/// Scales the subdomain constraints.
		/// This does nothing at all, as it is written right now.
		/// </summary>
		private void ScaleSubdomainConstraints(int currentIncrement)
		{
			if (currentIncrement == 0)
				return;

			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				//int idx = FindSubdomainIdx(linearSystems, linearSystem);
				double scalingFactor = 1; // ((double)currentIncrement + 2) / (currentIncrement + 1);
				subdomainUpdaters[linearSystem.Subdomain.ID].ScaleConstraints(scalingFactor);
			}
		}

		public class Builder : NonLinearAnalyzerBuilderBase
		{
			public Builder(IModel model, ISolver solver, INonLinearProvider provider, int numIncrements) :
				base(model, solver, provider, numIncrements)
			{
				MaxIterationsPerIncrement = 1000;
				NumIterationsForMatrixRebuild = 1;
				ResidualTolerance = 1E-3;
			}

			public DisplacementControlAnalyzer Build()
			{
				return new DisplacementControlAnalyzer(model, solver, provider, SubdomainUpdaters,
					numIncrements, maxIterationsPerIncrement, numIterationsForMatrixRebuild, residualTolerance);
			}
		}
	}
}
