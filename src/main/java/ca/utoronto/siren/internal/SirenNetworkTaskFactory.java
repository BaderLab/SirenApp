package ca.utoronto.siren.internal;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.task.NetworkTaskFactory;
import org.cytoscape.work.TaskIterator;

public class SirenNetworkTaskFactory implements NetworkTaskFactory {

	@Override
	public TaskIterator createTaskIterator(CyNetwork network) {
		return new TaskIterator(new SirenTask(network));
	}

	@Override
	public boolean isReady(CyNetwork network) {
		return network != null;
	}
	
}
