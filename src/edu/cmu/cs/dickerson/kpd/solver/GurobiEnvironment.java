package edu.cmu.cs.dickerson.kpd.solver;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;

/**
 * This class provides a static method {@code getGrbEnv()} that returns
 * a unique GRBEnv instance that is used throught the Java program's run
 * @author James Trimble
 *
 */
public class GurobiEnvironment {
	private static GRBEnv grbEnv = null;
	
	/**
	 * Returns a unique Gurobi environment instance.
	 * @return a unique Gurobi environment instance
	 * @throws GRBException
	 */
	public static GRBEnv getGrbEnv() throws GRBException {
		if (grbEnv == null) {
			grbEnv = new GRBEnv();
			//grbEnv.set(GRB.IntParam.Method, 2);
			int cores = Runtime.getRuntime().availableProcessors();
			grbEnv.set(GRB.IntParam.Threads, Math.max(1, cores-1));
		}
		
		return grbEnv;
	}
}
