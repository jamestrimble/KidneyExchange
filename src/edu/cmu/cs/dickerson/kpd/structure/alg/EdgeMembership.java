package edu.cmu.cs.dickerson.kpd.structure.alg;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.cmu.cs.dickerson.kpd.structure.Edge;
import edu.cmu.cs.dickerson.kpd.structure.Pool;
import edu.cmu.cs.dickerson.kpd.structure.Vertex;

public class EdgeMembership {
	private Map<Vertex, Set<Integer>> membership;
	
	// incSrc should be true if we want to include edges where a given vertex is the source
	// incTgt should be true if we want to include edges where a given vertex is the target
	public EdgeMembership(Pool pool, List<Edge> edges, boolean incSrc, boolean incTgt, boolean incAltruists) {
		// Want to make sure EVERY vertex in the pool is in this map
		membership = new HashMap<Vertex, Set<Integer>>();
		for(Vertex v : pool.vertexSet()) {
			membership.put(v, new HashSet<Integer>());
		}

		int edgeIdx = 0;
		for(Edge e : edges) {
			if (incSrc) {
				Vertex v = pool.getEdgeSource(e);
				if (incAltruists || !v.isAltruist()) membership.get(v).add(edgeIdx);
			}
			if (incTgt) {
				Vertex v = pool.getEdgeTarget(e);
				if (!v.isAltruist()) membership.get(v).add(edgeIdx);
			}
			edgeIdx++;
		}
	}

	public Set<Integer> getMembershipSet(Vertex v) {
		return membership.get(v);
	}
	
	public Set<Vertex> getAllVertices() {
		return membership.keySet();
	}
}