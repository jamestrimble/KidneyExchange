package ilog.cplex;

import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex.DoubleParam;
import ilog.cplex.IloCplex.IntParam;

public class IloCplex {

	public static enum DoubleParam {
		TiLim}

	public static enum IntParam {
		Threads};

	public void addLe(IloLinearNumExpr sum, double d) {
		// TODO Auto-generated method stub
		
	}

	public IloNumVar[] boolVarArray(int size) {
		// TODO Auto-generated method stub
		return null;
	}

	public void clearModel() {
		// TODO Auto-generated method stub
		
	}

	public void addEq(IloLinearNumExpr sum, double d) {
		// TODO Auto-generated method stub
		
	}

	public void addGe(IloLinearNumExpr sum, int exactRequiredNumPairs) {
		// TODO Auto-generated method stub
		
	}

	public int getNcols() {
		// TODO Auto-generated method stub
		return 0;
	}

	public String getCplexStatus() {
		// TODO Auto-generated method stub
		return null;
	}

	public String getObjValue() {
		// TODO Auto-generated method stub
		return null;
	}

	public Object scalProd(double[] weights, IloNumVar[] x) {
		// TODO Auto-generated method stub
		return null;
	}

	public IloLinearNumExpr linearNumExpr() {
		// TODO Auto-generated method stub
		return null;
	}

	public boolean solve() {
		// TODO Auto-generated method stub
		return false;
	}

	public void addMaximize(Object scalProd) {
		// TODO Auto-generated method stub
		
	}

	public void addGe(Object scalProd, double d) {
		// TODO Auto-generated method stub
		
	}

	public void addMinimize(Object scalProd) {
		// TODO Auto-generated method stub
		
	}

	public Object getObjective() {
		// TODO Auto-generated method stub
		return null;
	}

	public String getStatus() {
		// TODO Auto-generated method stub
		return null;
	}

	public double[] getValues(IloNumVar[] x) {
		// TODO Auto-generated method stub
		return null;
	}

	public IloNumVar[] numVarArray(int size, double d, double e) {
		// TODO Auto-generated method stub
		return null;
	}

	public void setParam(IntParam threads, int maxCPUThreads) {
		// TODO Auto-generated method stub
		
	}

	public void setParam(DoubleParam tilim, double maxSolveSeconds) {
		// TODO Auto-generated method stub
		
	}

}
