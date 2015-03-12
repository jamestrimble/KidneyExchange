package edu.cmu.cs.dickerson.kpd.james;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import org.jgrapht.generate.GraphGenerator;
import org.jgrapht.generate.RandomGraphGenerator;
import org.junit.Test;

import edu.cmu.cs.dickerson.kpd.fairness.solver.FairnessCPLEXSolver;
import edu.cmu.cs.dickerson.kpd.fairness.solver.FairnessGurobiSolver;
import edu.cmu.cs.dickerson.kpd.solver.CycleFormulationCPLEXSolver;
import edu.cmu.cs.dickerson.kpd.solver.CycleFormulationGurobiSolver;
import edu.cmu.cs.dickerson.kpd.solver.exception.SolverException;
import edu.cmu.cs.dickerson.kpd.solver.solution.Solution;
import edu.cmu.cs.dickerson.kpd.structure.Cycle;
import edu.cmu.cs.dickerson.kpd.structure.Edge;
import edu.cmu.cs.dickerson.kpd.structure.Pool;
import edu.cmu.cs.dickerson.kpd.structure.Vertex;
import edu.cmu.cs.dickerson.kpd.structure.alg.CycleGenerator;
import edu.cmu.cs.dickerson.kpd.structure.alg.CycleMembership;
import edu.cmu.cs.dickerson.kpd.structure.alg.FailureProbabilityUtil;
import edu.cmu.cs.dickerson.kpd.structure.generator.SaidmanPoolGenerator;
import edu.cmu.cs.dickerson.kpd.structure.generator.factories.AllMatchVertexPairFactory;

public class FailureAwareExperiment {

	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
		int numPairs = 50;
		int numAlts = numPairs/10;

		PrintWriter writer = new PrintWriter("results/results.txt", "UTF-8");
		
		// Test 10 random graphs, make sure the fairness solver (with no fairness constraint) and the base solver return the same value of answer
		for(int repeatIdx=0; repeatIdx<10; repeatIdx++) {
			Pool pool = new SaidmanPoolGenerator(new Random()).generate(numPairs, numAlts);

			FailureProbabilityUtil.ProbabilityDistribution failDist = FailureProbabilityUtil.ProbabilityDistribution.CONSTANT;
			FailureProbabilityUtil.setFailureProbability(pool, failDist, new Random(1), 0.7);			
			
			pool.writeToWmdFile("instances/instance" + repeatIdx);
			
			CycleGenerator cg = new CycleGenerator(pool);
			List<Cycle> cycles = cg.generateCyclesAndChains(3, 4, true, false, 1);

			CycleMembership membership = new CycleMembership(pool, cycles);

			try {
				Solution cSol = new CycleFormulationGurobiSolver(pool, cycles, membership).solve();
				writer.println("instance" + repeatIdx + "," + cSol.getObjectiveValue());

			} catch(SolverException e) {
				fail(e.getMessage());
			}
		}
		
		writer.close();
	}
}
