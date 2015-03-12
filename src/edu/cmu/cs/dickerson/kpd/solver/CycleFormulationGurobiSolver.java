package edu.cmu.cs.dickerson.kpd.solver;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.cmu.cs.dickerson.kpd.helper.IOUtil;
import edu.cmu.cs.dickerson.kpd.solver.exception.SolverException;
import edu.cmu.cs.dickerson.kpd.solver.solution.Solution;
import edu.cmu.cs.dickerson.kpd.structure.Cycle;
import edu.cmu.cs.dickerson.kpd.structure.Edge;
import edu.cmu.cs.dickerson.kpd.structure.Pool;
import edu.cmu.cs.dickerson.kpd.structure.Vertex;
import edu.cmu.cs.dickerson.kpd.structure.alg.CycleMembership;
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

public class CycleFormulationGurobiSolver {
	
	private CycleMembership membership;
	protected List<Cycle> cycles;
	private Pool pool;
	
	public CycleFormulationGurobiSolver(Pool pool, List<Cycle> cycles, CycleMembership membership) {
		this.pool = pool;
		this.cycles = cycles;
		this.membership = membership;
	}

	public Solution solve() throws SolverException {

		IOUtil.dPrintln(getClass().getSimpleName(), "Solving cycle formulation IP.");

		// If no cycles, problem is possibly unbounded; return 0-value empty solution
		if(cycles.size() == 0) {
			return new Solution(0,0,new HashSet<Cycle>());
		}
		
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


			model.dispose();

			return sol;

		} catch (GRBException e) {
			System.err.println("Exception thrown during CPLEX solve: " + e);
			throw new SolverException(e.toString());
		}
	}

	public String getID() {
		return "Cycle Formulation CPLEX Solver";
	}

}
