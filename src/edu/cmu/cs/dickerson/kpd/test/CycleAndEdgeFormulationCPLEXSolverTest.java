package edu.cmu.cs.dickerson.kpd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.util.List;
import java.util.Random;

import org.junit.Test;

import edu.cmu.cs.dickerson.kpd.solver.CycleAndEdgeFormulationCPLEXSolver;
import edu.cmu.cs.dickerson.kpd.solver.CycleFormulationCPLEXSolver;
import edu.cmu.cs.dickerson.kpd.solver.exception.SolverException;
import edu.cmu.cs.dickerson.kpd.solver.solution.Solution;
import edu.cmu.cs.dickerson.kpd.structure.Cycle;
import edu.cmu.cs.dickerson.kpd.structure.Pool;
import edu.cmu.cs.dickerson.kpd.structure.alg.CycleGenerator;
import edu.cmu.cs.dickerson.kpd.structure.alg.CycleMembership;
import edu.cmu.cs.dickerson.kpd.structure.alg.FailureProbabilityUtil;
import edu.cmu.cs.dickerson.kpd.structure.generator.SaidmanPoolGenerator;

public class CycleAndEdgeFormulationCPLEXSolverTest {

	@Test
	public void testNewSolver() {
		int numPairs = 50;
		int numAlts = numPairs/10;

		final int MAX_CHAIN = 4;
		
		// For a few instances, check that the new IP model gives the same objective result
		// as the cycle formulation
		for(int repeatIdx=0; repeatIdx<3; repeatIdx++) {
			Pool pool = new SaidmanPoolGenerator(new Random(repeatIdx)).generate(numPairs, numAlts);

			FailureProbabilityUtil.ProbabilityDistribution failDist = FailureProbabilityUtil.ProbabilityDistribution.CONSTANT;
			FailureProbabilityUtil.setFailureProbability(pool, failDist, new Random(1), 0.7);			
			
			CycleGenerator cg = new CycleGenerator(pool);
			List<Cycle> cycles = cg.generateCyclesAndChains(3, MAX_CHAIN, true, false, 1);

			CycleMembership membership = new CycleMembership(pool, cycles);

			try {
				Solution newSol = new CycleAndEdgeFormulationCPLEXSolver(pool, cycles).solve(MAX_CHAIN);
				Solution cSol = new CycleFormulationCPLEXSolver(pool, cycles, membership).solve();
				assertEquals(newSol.getObjectiveValue(), cSol.getObjectiveValue(), 1e-11);

			} catch(SolverException e) {
				fail(e.getMessage());
			}
		}
		
	}
}
