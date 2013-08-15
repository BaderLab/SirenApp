package ca.utoronto.siren.internal;

import org.junit.Assert;
import org.junit.Test;


public class SirenTest {
	static final double MAX_ERROR = 0.0001;
	
	@Test
	public void testScaleAndCentre() {
		double[] result = Siren.scaleAndCentre(new double[] { 1, 2, 3 });
		Assert.assertArrayEquals(new double[] { -1, 0, 1 }, result, MAX_ERROR);
	}
	
	@Test
	public void testScaleAndCentre2() {
		double[] result = Siren.scaleAndCentre(new double[] { 5, 4, 3, 2, 1, 1, 2, 3, 4, 5 });
		Assert.assertArrayEquals(new double[] { 1.3416408,
				0.6708204,
				0.0000000,
				-0.6708204,
				-1.3416408,
				-1.3416408,
				-0.6708204,
				0.0000000,
				0.6708204,
				1.3416408,
		}, result, MAX_ERROR);
	}
	
	@Test
	public void testComputeKnots() {
		double[] result = Siren.computeKnots(new double[] { 1, 2, 3, 4, 5, 6, 7, 8 }, 8, 2, 1, 8);
		Assert.assertArrayEquals(new double[] { 1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8 }, result, MAX_ERROR);
	}

	@Test
	public void testComputeKnots2() {
		// Overlapping knots
		double[] result = Siren.computeKnots(new double[] { 1, 2, 3, 4, 1, 2, 3, 4 }, 8, 2, 1, 4);
		Assert.assertArrayEquals(new double[] { 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4 }, result, MAX_ERROR);
	}

	@Test
	public void testComputeKnots3() {
		// Interpolated knots
		double[] result = Siren.computeKnots(new double[] { 0, 10 }, 11, 2, 0, 10);
		Assert.assertArrayEquals(new double[] { 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10 }, result, MAX_ERROR);
	}
}
