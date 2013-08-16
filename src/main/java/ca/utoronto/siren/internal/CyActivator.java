package ca.utoronto.siren.internal;

import java.util.Properties;

import org.cytoscape.service.util.AbstractCyActivator;
import org.cytoscape.task.NetworkTaskFactory;
import org.cytoscape.work.ServiceProperties;
import org.osgi.framework.BundleContext;

public class CyActivator extends AbstractCyActivator {
    public void start(BundleContext context) {
    	Properties properties = new Properties();
    	properties.put(ServiceProperties.PREFERRED_MENU, ServiceProperties.APPS_MENU);
    	properties.put(ServiceProperties.TITLE, "SIREN");
		registerService(context, new SirenNetworkTaskFactory(), NetworkTaskFactory.class, properties);
    }
}
