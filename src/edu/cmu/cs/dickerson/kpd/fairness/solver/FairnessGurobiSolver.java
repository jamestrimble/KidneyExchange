package edu.cmu.cs.dickerson.kpd.fairness.solver;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.cmu.cs.dickerson.kpd.helper.IOUtil;
import edu.cmu.cs.dickerson.kpd.solver.CPLEXSolver;
import edu.cmu.cs.dickerson.kpd.solver.GurobiEnvironment;
import edu.cmu.cs.dickerson.kpd.solver.exception.SolverException;
import edu.cmu.cs.dickerson.kpd.solver.solution.Solution;
import edu.cmu.cs.dickerson.kpd.structure.Cycle;
import edu.cmu.cs.dickerson.kpd.structure.Edge;
import edu.cmu.cs.dickerson.kpd.structure.Pool;
import edu.cmu.cs.dickerson.kpd.structure.Vertex;
import edu.cmu.cs.dickerson.kpd.structure.alg.CycleMembership;
import edu.cmu.cs.dickerson.kpd.structure.alg.FailureProbabilityUtil;
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

public class FairnessGurobiSolver {

	private double[] altWeights;
	private Set<Vertex> specialV;
	private CycleMembership membership;
	protected List<Cycle> cycles;
	private Pool pool;
	
	public FairnessGurobiSolver(Pool pool, List<Cycle> cycles, CycleMembership membership, Set<Vertex> specialV) {
		// Use the normal fairness solver without any failure probabilities
		this(pool, cycles, membership, specialV, false);
	}
	
	public FairnessGurobiSolver(Pool pool, List<Cycle> cycles, CycleMembership membership, Set<Vertex> specialV, boolean usingFailureProbabilities) {
		this.pool = pool;
		this.cycles = cycles;
		this.membership = membership;
		this.specialV = specialV;
		
		// Calculate highly-sensitized weights for each cycle
		if(usingFailureProbabilities) {
			calcSpecialWeightsWithFailureProbs();
		} else {
			calcSpecialWeightsWithoutFailureProbs();
		}	
		
	}

	/**
	 * Computes h(c), the non-failure-aware weight of a cycle where only transplants to 
	 * a vertex in the "special vertex" (e.g., highly-sensitized) set count
	 */
	private void calcSpecialWeightsWithoutFailureProbs() {
		
		// Each cycle will have a new, adjusted weight
		altWeights = new double[cycles.size()];
		int cycleIdx = 0;

		// For each cycle, create a new weight that only takes special vertices' successful transplants into account
		for(Cycle c : cycles) {
			double altWeight = 0.0;
			for(Edge e : c.getEdges()) {
				Vertex recipient = pool.getEdgeTarget(e);
				if(specialV.contains(recipient)) {
					altWeight += pool.getEdgeWeight(e);
				}
			}
			altWeights[cycleIdx++] = altWeight;
		}
	}
	
	
	/**
	 * Computes h(c), the discounted (failure-aware) utility of a cycle where only transplants to 
	 * a vertex in the "special vertex" (e.g., highly-sensitized) set count
	 */
	private void calcSpecialWeightsWithFailureProbs() {
		
		// Each cycle will have a new, adjusted weight
		altWeights = new double[cycles.size()];
		int cycleIdx = 0;

		// For each cycle, create a new utility that only takes special vertices' successful transplants into account
		for(Cycle c : cycles) {
			
			boolean isChain = Cycle.isAChain(c, pool);
			
			if(isChain) {
				altWeights[cycleIdx++] = FailureProbabilityUtil.calculateDiscountedChainUtility(c, pool, specialV);
			} else {
				altWeights[cycleIdx++] = FailureProbabilityUtil.calculateDiscountedCycleUtility(c, pool, specialV);
			}
			
		}
	}
	
	
	
	
	public Solution solve(double alpha) throws SolverException {
		
		IOUtil.dPrintln(getClass().getSimpleName(), "Solving main fairness IP with a = " + alpha);
		
		try {
			Solution sol = new Solution();
			
			GRBEnv env = GurobiEnvironment.getGrbEnv();
			GRBModel model = new GRBModel(env);
			env.set(GRB.DoubleParam.MIPGap, 0);
			
			// One decision variable per cycle
			GRBVar[] x = new GRBVar[cycles.size()];
			for (int i=0; i<cycles.size(); i++) {
				x[i] = model.addVar(0, 1, 0, GRB.BINARY, "x" + i);
			}
			
			// Decision variables multiplied by weight of corresponding cycle
			double[] weights = new double[x.length];
			int cycleIdx = 0;
			for(Cycle c : cycles) {
				weights[cycleIdx++] = c.getWeight();
			}

			
			// Objective:
			// Maximize \sum_{all cycles c} altWeight_c * decVar_c
			GRBLinExpr scalProdExpr = new GRBLinExpr();
			for (int i=0; i<cycles.size(); i++) {
				scalProdExpr.addTerm(weights[i], x[i]);
			}
			
			model.update();
			
			model.setObjective(scalProdExpr, GRB.MAXIMIZE);					
			
			model.update();
			
			
			// Subject to: 
			// \sum_{cycles c containing v} decVar_c <=1   \forall v
			for(Vertex v : pool.vertexSet()) {
				
				Set<Integer> cycleColIDs = membership.getMembershipSet(v);
				if(null == cycleColIDs || cycleColIDs.isEmpty()) {
					continue;
				}
				
				GRBLinExpr sum = new GRBLinExpr();
				for(Integer cycleColID : cycleColIDs) {
					sum.addTerm(1.0, x[cycleColID]);
				}
				model.addConstr(sum, GRB.LESS_EQUAL, 1.0, "");
			}
			
			
			// \sum_c altWeight_c * decVar_c >= alpha*|special|
			GRBLinExpr linExpr = new GRBLinExpr();
			for (int i=0; i<cycles.size(); i++) {
				linExpr.addTerm(altWeights[i], x[i]);
			}
			model.addConstr(linExpr, GRB.GREATER_EQUAL, alpha * specialV.size(), "");
			

			// Solve the model, get base statistics (solve time, objective value, etc)
			model.update();
			model.optimize();
			double objVal = model.get(GRB.DoubleAttr.ObjVal);
			sol.setObjectiveValue(objVal);
			sol.setSolveTime(0);
		
			// Figure out which cycles were included in the final solution
			double[] vals = new double[x.length];
			for (int i=0; i<x.length; i++) {
				vals[i] = x[i].get(GRB.DoubleAttr.X);
			}
					
			int nCols = x.length; // Possibly wrong? Was cplex.getNcols();
			for(cycleIdx=0; cycleIdx<nCols; cycleIdx++) {
				if(vals[cycleIdx] > 1e-3) {
					sol.addMatchedCycle(cycles.get(cycleIdx));
				}
			}
			
			IOUtil.dPrintln(getClass().getSimpleName(), "Solved IP!  Objective value: " + sol.getObjectiveValue());
			IOUtil.dPrintln(getClass().getSimpleName(), "Number of cycles in matching: " + sol.getMatching().size());
			
			// TODO move to a JUnit test
			// Sanity check to make sure the matching is vertex disjoint
			Set<Vertex> seenVerts = new HashSet<Vertex>();
			for(Cycle c : sol.getMatching()) {
				for(Edge e : c.getEdges()) {
					Vertex v = pool.getEdgeSource(e);
					if(seenVerts.contains(v)) {
						IOUtil.dPrintln(getClass().getSimpleName(), "A vertex (" + v + ") was in more than one matched cycle; aborting.");
					}
					seenVerts.add(v);
				}
			}
			
			
			
			// 
			model.dispose();
			//cplex.end();		
			
			return sol;
			
		} catch (GRBException e) {
			System.err.println("Exception thrown during CPLEX solve: " + e);
			throw new SolverException(e.toString());
		}
	}

	
	public String getID() {
		return "Fairness CPLEX Solver";
	}

	
}
