package ca.utoronto.siren.internal;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Java port of the SIgning Of REgulatory Networks (SIREN) scoring algorithm
 * (originally from R).  This class computes the SIREN scores for a gene-gene
 * network given some expression data.
 */
public class Siren {
	public static double[] computeScores(double[][] expressionMatrix, double[][] networkMatrix, double[][] weightMatrix) {
		int degreesOfFreedom = 10;
		int degree = 2;
		
		double[][][] bMatrix = computeBMatrix(expressionMatrix, degreesOfFreedom, degree);
		double[][] paMatrix = computePaMatrix(bMatrix);
		return computeScores(bMatrix, weightMatrix, paMatrix, networkMatrix);
	}
	
	static double[][][] computeBMatrix(double[][] expressionMatrix, int degreesOfFreedom, int degree) {
		// Assume expressionMatrix is rectangular and has at least 1 row
		int totalConditions = expressionMatrix[0].length;
		int totalGenes = expressionMatrix.length;
		
		double[][][] result = new double[totalGenes][degreesOfFreedom][totalConditions];
		for (int g = 0; g < totalGenes; g++) {
			double[][] basis = computeBSplineBasis(scaleAndCentre(expressionMatrix[g]), degreesOfFreedom, degree);
			for (int c = 0; c < totalConditions; c++) {
				for (int i = 0; i < basis.length; i++) {
					result[g][i][c] = basis[i][c];
				}
			}
		}
		return result;
	}
	
	static double[][] computePaMatrix(double[][][] bMatrix) {
		int totalGenes = bMatrix.length;
		int totalBins = bMatrix[0].length;
		int totalConditions = bMatrix[0][0].length;
		
		double[][] result = new double[totalGenes][totalBins];
		for (int g = 0; g < totalGenes; g++) {
			for (int b = 0; b < totalBins; b++) {
				for (int c = 0; c < totalConditions; c++) {
					result[g][b] += bMatrix[g][b][c];
				}
				result[g][b] /= totalConditions;
			}
		}
		return result;
	}
	
	static double[][] computePabMatrix(double[][] bMatrixA, double[][] bMatrixB) {
		int totalBins = bMatrixA.length;
		int totalConditions = bMatrixA[0].length;
		
		double[][] result = new double[totalBins][totalBins];
		for (int i = 0; i < totalBins; i++) {
			for (int j = 0; j < totalBins; j++) {
				for (int c = 0; c < totalConditions; c++) {
					result[i][j] += bMatrixA[i][c] * bMatrixB[j][c];
				}
				result[i][j] /= totalConditions;
			}
		}
		return result;
	}

	static double[] scaleAndCentre(double[] vector) {
		// Adapted from online_variance() at:
		// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
		
		int count = 0;
		double mean = 0;
		double M2 = 0;
		
		for (int i = 0; i < vector.length; i++) {
			if (Double.isNaN(vector[i])) {
				continue;
			}
			count++;
			double delta = vector[i] - mean;
			mean += delta / count;
			M2 += delta * (vector[i] - mean);
		}
		
		double variance = M2 / (count - 1);
		double sigma = Math.sqrt(variance);
		
		double[] result = new double[vector.length];
		for (int i = 0; i < vector.length; i++) {
			if (Double.isNaN(vector[i])) {
				result[i] = Double.NaN;
				continue;
			}
			result[i] = (vector[i] - mean) / sigma;
		}
		return result;
	}
	
	/**
	 * Returns a vector of length X containing the coefficients of polynomial
	 * j.
	 * 
	 * @param x
	 * @param degree
	 * @param j
	 * @param knots
	 * @return
	 */
	static double[] computeBasis(double[] x, int degreesOfFreedom, int degree, int j, double[] knots) {
		// Adapted from basis() at:
		// http://cran.r-project.org/web/packages/crs/vignettes/spline_primer.pdf
		
		double[] result = new double[x.length];

		// The R implementation has NaNs in the following case.
		double maxX = knots[knots.length - 1];
		boolean hasMultipleOuterBound = maxX == knots[degreesOfFreedom + degree - 2];
		for (int i = 0; i < x.length; i++) {
			if (x[i] == maxX && j > degreesOfFreedom - degree && hasMultipleOuterBound) {
				result[i] = Double.NaN;
			}
		}

		if (degree == 0) {
			for (int i = 0; i < x.length; i++) {
				if (Double.isNaN(result[i])) {
					continue;
				}
				result[i] = (x[i] >= knots[j] && x[i] < knots[j + 1]) ? 1 : 0;
			}
			return result;
		} else {
			double[] basis1 = computeBasis(x, degreesOfFreedom, degree - 1, j, knots);
			double[] basis2 = computeBasis(x, degreesOfFreedom, degree - 1, j + 1, knots);
			
			for (int i = 0; i < x.length; i++) {
				if (Double.isNaN(result[i])) {
					continue;
				}
				double denominator1 = knots[degree + j] - knots[j];
				double denominator2 = knots[j + degree + 1] - knots[j + 1];				
				double alpha1 = denominator1 == 0 ? 0 : (x[i] - knots[j]) / denominator1;
				double alpha2 = denominator2 == 0 ? 0 : (knots[j + degree + 1] - x[i]) / denominator2;
				result[i] = alpha1 * basis1[i] + alpha2 * basis2[i];
			}
			return result;
		}
	}
	
	/**
	 * Returns an X x K B-spline basis matrix, where X is the number of
	 * elements in x, and K = degreesOfFreedom.
	 * @param x
	 * @param degreesOfFreedom
	 * @param degree
	 * @return
	 */
	static double[][] computeBSplineBasis(double[] x, int degreesOfFreedom, int degree) {
		// Adapted from bs() at:
		// http://cran.r-project.org/web/packages/crs/vignettes/spline_primer.pdf

		double minX = Double.MAX_VALUE;
		double maxX = Double.MIN_VALUE;
		
		for (int i = 0; i < x.length; i++) {
			minX = Math.min(minX, x[i]);
			maxX = Math.max(maxX, x[i]);
		}
		
		int k = degreesOfFreedom;
		double[] knots = computeKnots(x, degreesOfFreedom, degree, minX, maxX);
		double[][] result = new double[k][];
		
		// Discard intercept; i.e. process bases j = 1...k-1
		for (int j = 0; j < k; j++) {
			result[j] = computeBasis(x, degreesOfFreedom, degree, j + 1, knots);
		}
		
		boolean hasMultipleOuterBound = maxX == knots[degreesOfFreedom + degree - 2];
		for (int i = 0; i < x.length; i++) {
			if (x[i] == maxX && !hasMultipleOuterBound) {
				result[k - 1][i] = 1; 
			}
		}
		
		return result;
	}
	
	static double[] computeKnots(double sample[], int degreesOfFreedom, int degree, double minX, double maxX) {
		int interiorKnotCount = degreesOfFreedom - degree;
		double[] quantiles = computeQuantiles(sample, interiorKnotCount + 1);
		double[] knots = new double[interiorKnotCount + 2 * (degree + 1)];
		for (int i = 0; i < interiorKnotCount + 1; i++) {
			knots[i + degree] = quantiles[i];
		}
		for (int i = 0; i < degree + 1; i++) {
			knots[i] = minX;
			knots[knots.length - 1 - i] = maxX;
		}
		return knots;
	}

	/**
	 * Returns the kth q-quantiles for q = bins and 0 < k < q.
	 * @param sample sample data
	 * @param bins total number of quantiles
	 * @return the kth q-quantiles for q = bins and 0 < k < q.
	 */
	static double[] computeQuantiles(double[] sample, int bins) {
		double[] x = new double[sample.length];
		System.arraycopy(sample, 0, x, 0, sample.length);
		Arrays.sort(x);
		
		double[] result = new double[bins];
		int q = bins;
		int N = x.length;
		for (int k = 1; k < q; k++) {
			double p = (double) k / q;
			double h = (N - 1) * p + 1;
			int floorH = (int) Math.floor(h);
			result[k] = x[floorH - 1] + (h - floorH) * (x[floorH] - x[floorH - 1]);
		}
		return result;
	}
	
	public static double[][] loadMatrix(String path) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(path));
		int rows = 0;
		int columns = 0;
		
		try {
			String line = reader.readLine();
			while (line != null) {
				try {
					if (columns == 0) {
						columns = line.split("\t").length;
					}
					rows++;
				} finally {
					line = reader.readLine();
				}
			}
		} finally {
			reader.close();
		}

		reader = new BufferedReader(new FileReader(path));
		double[][] result = new double[rows][columns];
		try {
			int rowIndex = 0;
			String line = reader.readLine();
			while (line != null) {
				try {
					int columnIndex = 0;
					for (String value : line.split("\t")) {
							result[rowIndex][columnIndex] = Double.parseDouble(value);
							columnIndex++;
					}
					rowIndex++;
				} finally {
					line = reader.readLine();
				}
			}
		} finally {
			reader.close();
		}
		return result;
	}
	
	static void printVector(double[] vector) {
		for (int i = 0; i < vector.length; i++) {
			System.out.printf("%f\n", vector[i]);
		}
	}

	static void printMatrix(double[][] matrix) {
		int rows = matrix.length;
		int columns = matrix[0].length;
		
		for (int j = 0; j < columns; j++) {
			for (int i = 0; i < rows; i++) {
				System.out.printf("%f\t", matrix[i][j]);
			}
			System.out.println();
		}
		System.out.println();
	}

	static void printMatrix(double[][][] matrix) {
		int layers = matrix.length;
		for (int i = 0; i < layers; i++) {
			System.out.printf("(%d, ,)\n\n", i);
			printMatrix(matrix[i]);
		}
	}
	
	static double[] computeScores(double[][][] bMatrix, double[][] weightMatrix, double[][] paMatrix, double[][] networkMatrix) {
		int totalBins = bMatrix[0].length;
		int totalInteractions = networkMatrix.length;
		
		double[] result = new double[totalInteractions];
		for (int i = 0; i < totalInteractions; i++) {
			int geneA = (int) networkMatrix[i][0];
			int geneB = (int) networkMatrix[i][1];
			
			double[] pA = paMatrix[geneA];
			double[] pB = paMatrix[geneB];
			double[][] pABMatrix = computePabMatrix(bMatrix[geneA], bMatrix[geneB]);
			for (int x = 0; x < totalBins; x++) {
				for (int y = 0; y < totalBins; y++) {
					double xN = Math.log(pABMatrix[x][y] / pA[x] / pB[y]);
					if (xN > 0) {
						result[i] += pABMatrix[x][y] * weightMatrix[x][y] * xN;
					}
				}
			}
		}
		return result;
	}

	static List<Double> parseResults(String filename) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		try {
			List<Double> result = new ArrayList<Double>();
			String line = reader.readLine();
			while (line != null) {
				try {
					String[] parts = line.split("\t");
					if (parts.length == 4) {
						result.add(Double.parseDouble(parts[3]));
					}
				} finally {
					line = reader.readLine();
				}
			}
			return result;
		} finally {
			reader.close();
		}
	}
	
	public static void main(String[] args) throws IOException {
		double[][] weightMatrix = loadMatrix("Weighting_Matrix.txt");
		double[][] expressionMatrix = loadMatrix("Expression_Format.txt");
		double[][] networkMatrix = loadMatrix("Network_Format.txt");
		adjustIndexes(networkMatrix);

		double[] scores = computeScores(expressionMatrix, networkMatrix, weightMatrix);
		List<Double> results = parseResults("Result.txt");
		for (int i = 0; i < scores.length; i++) {
			if (Math.abs(scores[i] - results.get(i)) > 0.0000001) {
				System.out.printf("%d\t%g\n", i, scores[i] - results.get(i));
			}
		}
	}

	static void adjustIndexes(double[][] networkMatrix) {
		// The R implementation of SIREN used networks with 1-based indexes.
		// Since Java uses 0-based indexes, we need to adjust them.
		for (int i = 0; i < networkMatrix[0].length; i++) {
			for (int j = 0; j < networkMatrix.length; j++) {
				networkMatrix[j][i]--;
			}
		}
	}
}
