package edu.cmu.cs.dickerson.kpd.solver;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;

import java.util.ArrayList;
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
import edu.cmu.cs.dickerson.kpd.structure.alg.EdgeMembership;

public class CycleAndEdgeFormulationCPLEXSolver extends CPLEXSolver {
	
	private CycleMembership membership;
	protected List<Cycle> cycles;
	private Pool pool;
	
	private List<Vertex> vertices;
	private List<Edge> edgesFromAltruist;
	private List<Edge> nonAltruistEdges;
	
	private EdgeMembership firstEdgeMembership;	// membership map for edges from altruistic donors
	private EdgeMembership otherEdgeSourceMembership; // for each non-altruistic vertex V, of which edges (not from an altruist) is V the source?
	private EdgeMembership otherEdgeTargetMembership; // for each non-altruistic vertex V, of which edges (not from an altruist) is V the target?
	
	public CycleAndEdgeFormulationCPLEXSolver(Pool pool, List<Cycle> cycles) {
		super(pool);
		this.pool = pool;
		this.cycles = findNonAltruisticCycles(cycles);
		membership = new CycleMembership(pool, this.cycles);
		
		vertices = new ArrayList<Vertex>();
		for (Vertex v : pool.vertexSet()) vertices.add(v);
		
		findEdgesFromAltruist();
		findNonAltruistEdges();
	}
	
	private List<Cycle> findNonAltruisticCycles(List<Cycle> allCycles) {
		List<Cycle> nonAltruisticCycles = new ArrayList<Cycle>();
		for (Cycle c : allCycles) {
			if (!Cycle.isAChain(c, pool)) {
				nonAltruisticCycles.add(c);
			}
		}
		return nonAltruisticCycles;
	}

	private void findEdgesFromAltruist() {
		edgesFromAltruist = new ArrayList<Edge>();
		for (Vertex v : vertices)
			if (v.isAltruist())
				for (Edge e : pool.outgoingEdgesOf(v))
					if (!pool.getEdgeTarget(e).isAltruist())
						edgesFromAltruist.add(e);
		firstEdgeMembership = new EdgeMembership(pool, edgesFromAltruist, true, true, true);
	}

	private void findNonAltruistEdges() {
		nonAltruistEdges = new ArrayList<Edge>();
		for (Vertex v : vertices)
			if (!v.isAltruist())
				for (Edge e : pool.outgoingEdgesOf(v))
					if (!pool.getEdgeTarget(e).isAltruist()) 
						nonAltruistEdges.add(e);
		otherEdgeSourceMembership = new EdgeMembership(pool, nonAltruistEdges, true, false, false);
		otherEdgeTargetMembership = new EdgeMembership(pool, nonAltruistEdges, false, true, false);
	}

	public Solution solve(int maxChain) throws SolverException {

		int maxChainTail = maxChain-1;	// The max number of non-altruistic edges in a chain
		
		IOUtil.dPrintln(getClass().getSimpleName(), "Solving cycle formulation IP.");

		// If no cycles, problem is possibly unbounded; return 0-value empty solution
		if(cycles.size() == 0) {
			return new Solution(0,0,new HashSet<Cycle>());
		}
		
		int n = nonAltruistEdges.size();
		
		try {
			super.initializeCPLEX();

			// One decision variable per cycle
			IloNumVar[] x = cplex.boolVarArray(cycles.size());

			// One decision variable per edge from an altruist
			IloNumVar[] x_a = cplex.boolVarArray(edgesFromAltruist.size());

			// One decision variable per non-altruist edge
			IloNumVar[][] x_n = new IloNumVar[maxChainTail][];
			for (int i=0; i<maxChainTail; i++)
				x_n[i] = cplex.boolVarArray(n);
			
				

			// Decision variables multiplied by weight of corresponding cycle
			double[] weights = new double[x.length];
			int cycleIdx = 0;
			for(Cycle c : cycles) weights[cycleIdx++] = c.getWeight();

			// Decision variables multiplied by weight of corresponding edge from an altruist
			double[] weights_a = new double[x_a.length];
			int edgeIdx = 0;
			for (Edge e : edgesFromAltruist) {
				weights_a[edgeIdx++] = pool.getEdgeWeight(e) * (1-e.getFailureProbability());
			}

			// Decision variables multiplied by weight of corresponding edge between two donor-patient pairs
			double[][] weights_n = new double[maxChainTail][n];
			edgeIdx = 0;
			for (Edge e : nonAltruistEdges) {
				for (int i=0; i<maxChainTail; i++) {
					weights_n[i][edgeIdx] = pool.getEdgeWeight(e) * Math.pow((1-e.getFailureProbability()), i+2);
				}
				edgeIdx++;
			}
			

			IloNumVar[] x_n_flattened = new IloNumVar[maxChainTail*n];
			double[] weights_n_flattened = new double[maxChainTail*n];
			for (int i=0; i<maxChainTail; i++) {
				for (int j=0; j<n; j++) { 
					x_n_flattened[i*n+j] = x_n[i][j];
					weights_n_flattened[i*n+j] = weights_n[i][j];
				}
			}


			// Objective:
			// Maximize sum of chosen cycle weights plus sum of chosen chain-edge weights
			cplex.addMaximize(cplex.sum(
					cplex.scalProd(weights, x),
					cplex.scalProd(weights_a, x_a),
					cplex.scalProd(weights_n_flattened, x_n_flattened)));


			// Subject to: 
			// \sum_{cycles/edges x containing v} decVar_x <=1   \forall v
			for(Vertex v : pool.vertexSet()) {

				IloLinearNumExpr sum = cplex.linearNumExpr(); 
				
				Set<Integer> cycleColIDs = membership.getMembershipSet(v);
				if(null != cycleColIDs)	
					for(Integer cycleColID : cycleColIDs)
						sum.addTerm(1.0, x[cycleColID]);

				Set<Integer> edgeIDs = firstEdgeMembership.getMembershipSet(v);
				if(null != edgeIDs)	
					for(Integer edgeID : edgeIDs) {
						sum.addTerm(1.0, x_a[edgeID]);
					}

				edgeIDs = otherEdgeTargetMembership.getMembershipSet(v);
				if(null != edgeIDs)	
					for(Integer edgeID : edgeIDs)
						for (int i=0; i<maxChainTail; i++)
							sum.addTerm(1.0, x_n[i][edgeID]);

				cplex.addLe(sum, 1.0);
			}
			
			
			// if a non-altruist vertex has an outgoing edge at position i+1 of
			// a chain, then it must have an outgoing edge at position i of a chain
			for(Vertex v : pool.vertexSet()) {
				if (!v.isAltruist()) {
										
					Set<Integer> inEdgeIDsFromAlt = firstEdgeMembership.getMembershipSet(v);
					Set<Integer> inEdgeIDs = otherEdgeTargetMembership.getMembershipSet(v);
					Set<Integer> outEdgeIDs = otherEdgeSourceMembership.getMembershipSet(v);
					
					for (int i=0; i<maxChainTail; i++) {
						IloLinearNumExpr lhs = cplex.linearNumExpr();
						IloLinearNumExpr rhs = cplex.linearNumExpr();
						
						if (i==0) {
							for(int edgeID : inEdgeIDsFromAlt)
								lhs.addTerm(1.0, x_a[edgeID]);
						} else {
							for(int edgeID : inEdgeIDs)
								lhs.addTerm(1.0, x_n[i-1][edgeID]);
						}
						
						for(int edgeID : outEdgeIDs)
							rhs.addTerm(1.0, x_n[i][edgeID]);
						
						cplex.addGe(lhs, rhs);
					}
				}
			}


			// Solve the model, get base statistics (solve time, objective value, etc)
			Solution sol = super.solveCPLEX();

			// Figure out which cycles were included in the final solution
			double[] vals = cplex.getValues(x);
			int nCols = cplex.getNcols();
			for(cycleIdx=0; cycleIdx<x.length; cycleIdx++) {
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
			cplex.clearModel();
			//cplex.end();		

			return sol;

		} catch(IloException e) {
			System.err.println("Exception thrown during CPLEX solve: " + e);
			throw new SolverException(e.toString());
		}
	}

	@Override
	public String getID() {
		return "Cycle and Edge Formulation CPLEX Solver";
	}

}
