package ca.utoronto.siren.internal;

import java.util.ArrayList;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

import org.cytoscape.model.CyColumn;
import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyRow;
import org.cytoscape.model.CyTable;
import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.ListMultipleSelection;

public class SirenTask extends AbstractTask {
	
	@Tunable(description="Select gene expression attribute(s)")
	public ListMultipleSelection<String> attributeNames;
	private CyNetwork network;
	
	public SirenTask(CyNetwork network) {
		this.network = network;
		
		CyTable table = network.getDefaultNodeTable();
		List<String> columnNames = new ArrayList<String>();
		for (CyColumn column : table.getColumns()) {
			Class<?> type = column.getType();
			if (Double.class.equals(type) || Integer.class.equals(type)) {
				columnNames.add(column.getName());
			}
		}
		Collections.sort(columnNames);
		attributeNames = new ListMultipleSelection<String>(columnNames);
	}
	
	@ProvidesTitle
	public String getTitle() {
		return "Compute SIREN scores for " + network.getRow(network).get(CyNetwork.NAME, String.class);
	}
	
	@Override
	public void run(TaskMonitor taskMonitor) throws Exception {
		List<String> columnNames = attributeNames.getSelectedValues();
		CyColumn[] columns = getColumns(network, columnNames);
		
		List<CyNode> nodes = network.getNodeList();
		List<CyEdge> edges = network.getEdgeList();
		
		int[][] networkMatrix = extractNetworkMatrix(nodes, edges);
		double[][] expressionMatrix = extractExpressionMatrix(network, nodes, columns);
		double[] scores = Siren.computeScores(expressionMatrix, networkMatrix, Siren.DEFAULT_WEIGHT_MATRIX);
		
		CyTable table = network.getDefaultEdgeTable();
		String columnName = "SIREN";
		if (table.getColumn(columnName) == null) {
			table.createColumn(columnName, Double.class, false);
		}
		
		for (int i = 0; i < scores.length; i++) {
			CyEdge edge = edges.get(i);
			network.getRow(edge).set(columnName, scores[i]);
		}
	}

	private static CyColumn[] getColumns(CyNetwork network, List<String> columnNames) {
		CyTable table = network.getDefaultNodeTable();
		CyColumn[] columns = new CyColumn[columnNames.size()];
		int index = 0;
		for (String name : columnNames) {
			columns[index++] = table.getColumn(name);
		}
		return columns;
	}

	@SuppressWarnings("unchecked")
	private static double[][] extractExpressionMatrix(CyNetwork network, List<CyNode> nodes, CyColumn[] columns) {
		double[][] result = new double[nodes.size()][columns.length];
		int nodeIndex = 0;
		for (CyNode node : nodes) {
			for (int i = 0; i < columns.length; i++) {
				CyColumn column = columns[i];
				CyRow row = network.getRow(node);
				Class<Number> type = (Class<Number>) column.getType();
				Number number = row.get(column.getName(), type);
				result[nodeIndex][i] = number == null ? Double.NaN : number.doubleValue(); 
			}
			nodeIndex++;
		}
		return result;
	}

	private static int[][] extractNetworkMatrix(List<CyNode> nodes, List<CyEdge> edges) {
		Map<CyNode, Integer> nodeIndexes = new IdentityHashMap<CyNode, Integer>();
		int nodeIndex = 0;
		for (CyNode node : nodes) {
			nodeIndexes.put(node, nodeIndex++);
		}

		int[][] result = new int[edges.size()][];
		int edgeIndex = 0;
		for (CyEdge edge : edges) {
			result[edgeIndex++] = new int[] {
				nodeIndexes.get(edge.getSource()),
				nodeIndexes.get(edge.getTarget()),
			};
		}
		return result;
	}
}
