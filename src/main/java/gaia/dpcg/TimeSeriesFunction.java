/*-----------------------------------------------------------------------------
 *
 *                      Gaia CU7 for Variability Processing
 *                           ObsGe Geneva Observatory
 *
 *        (c) 2005-2020 Gaia Data Processing and Analysis Consortium
 *
 *
 * CU7 variability software is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * CU7 variability software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this CU7 variability software; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301
 * USA
 *
 *-----------------------------------------------------------------------------
 */
package gaia.dpcg;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.analysis.function.HarmonicOscillator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.stat.regression.ModelSpecificationException;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathUtils;



/**
 * Representation of the general univariate time series model of variability.
 * defined by:
 * 
 * y(x) = a_0 + a_1 * x + a_2 * x^2 + ... + a_{n-1} x^(n-1) + sin[2*PI*f*(t-t0)]
 * + cos[2*PI*f*(t-t0)] + .... sin[(l-1) * 2*PI*l*f*(t-t0)] + cos[(l-1) *
 * 2*PI*l*f*(t-t0)]
 * 
 * where:
 * <ul>
 * <li>x are the covariates</li>
 * <li>y are the responses</li>
 * <li>a_i are the regression coefficients (to be determined)</li>
 * <li>f_i(x) are the regressors evaluated in the covariates.</li>
 * </ul>
 * <p>
 * The time series function is represented by a polynomial function together
 * with a list of 0 or more HarmonicOscillators
 * 
 * @since 9.0.0
 * @author Leanne Guy &lt;leanne.guy@unige.ch&gt;
 * @version $Id: TimeSeriesFunction.java 556096 2017-03-31 08:20:31Z lguy $
 * 
 *          Requirements implemented by this class
 * @req CU7-WP711-03000-S-FUN-200
 * @req CU7-WP711-03000-S-FUN-220
 * @req CU7-WP711-03000-S-FUN-240
 * @req CU7-WP711-03000-S-FUN-260
 */
public final class TimeSeriesFunction implements UnivariateDifferentiableFunction, Serializable {

	/** Serialization identifier */
	private static final long serialVersionUID = 910;

	/** The 2 * PI term */
	private static final double TWO_PI = 2.0 * Math.PI;

	/** The polynomial function that defines the model */
	private PolynomialFunction polynomialFunction;

	/** The harmonic functions that define the model */
	private List<HarmonicOscillator> harmonicFunctions;

	/**
	 * Construct a time series function with the given polynomial terms and no
	 * harmonic terms. The first element of the coefficients array is the
	 * constant term. Higher degree coefficients follow in sequence.
	 * 
	 * @param polynomialCoefficients
	 *            the polynomial coefficients
	 */
	public TimeSeriesFunction(final double[] polynomialCoefficients) {
		this(polynomialCoefficients, null);
	}

	/**
	 * Construct a time series function with the given polynomial terms and
	 * harmonic terms.
	 * <p>
	 * The first element of the polynomial coefficients array
	 * {@code polynomialCoefficients} is the constant term. Higher degree
	 * coefficients follow in sequence.
	 * <p>
	 * Each harmonic term {@code harmonicCoefficients} represents a single
	 * harmonic element with the coefficients ordered as:
	 * <ol>
	 * <li>amplitude</li>
	 * <li>frequency</li>
	 * <li>phase</li>
	 * </ol>
	 * 
	 * @param polynomialCoefficients
	 *            the polynomial coefficients array
	 * @param harmonicCoefficients
	 *            the harmonic coefficients 2D array
	 */
	public TimeSeriesFunction(double[] polynomialCoefficients, double[][] harmonicCoefficients) {
		this(polynomialCoefficients, harmonicCoefficients, 0.0);
	}

	/**
	 * Construct a time series function with the given polynomial terms and
	 * harmonic terms.
	 * <p>
	 * The first element of the polynomial coefficients array
	 * {@code polynomialCoefficients} is the constant term. Higher degree
	 * coefficients follow in sequence.
	 * <p>
	 * Each harmonic term {@code harmonicCoefficients} represents a single
	 * harmonic element with the coefficients ordered as:
	 * <ol>
	 * <li>amplitude</li>
	 * <li>frequency</li>
	 * <li>phase</li>
	 * </ol>
	 * 
	 * @param polynomialCoefficients
	 *            the polynomial coefficients array
	 * @param harmonicCoefficients
	 *            the harmonic coefficients 2D array
	 * @param referenceTime the function reference time
	 */
	public TimeSeriesFunction(double[] polynomialCoefficients, double[][] harmonicCoefficients, double referenceTime) {
		this.polynomialFunction = new PolynomialFunction(polynomialCoefficients);
		if (harmonicCoefficients != null) {
			this.harmonicFunctions = new ArrayList<HarmonicOscillator>();
			for (double[] harmonicCoefficient : harmonicCoefficients) {
				harmonicFunctions.add(new HarmonicOscillator(harmonicCoefficient[0], TWO_PI * harmonicCoefficient[1],
						harmonicCoefficient[2]));
			}
		}
	}
	
	/**
	 * Get the degree of the polynomial component
	 * 
	 * @return the polynomialdegree
	 */
	public int getPolynomialDegree() {
		return this.polynomialFunction.degree();
	}

	/**
	 * Get the polynomial component of the time series function
	 * 
	 * @return the polynomialFunction
	 */
	public PolynomialFunction getPolynomialFunction() {
		return this.polynomialFunction;
	}

	/**
	 * Return the list of all harmonic functions
	 * 
	 * @return the harmonicFunctions
	 */
	public List<HarmonicOscillator> getHarmonicFunctions() {
		return harmonicFunctions;
	}

	/**
	 * Add a harmonic function to the list of existing harmonic functions
	 * 
	 * @param a
	 *            teh amplitude
	 * @param frequency
	 *            the frequency
	 * @param phi
	 *            the phase
	 */
	public void addHarmonicFunction(final double a, final double frequency, final double phi) {
		if (harmonicFunctions == null) {
			harmonicFunctions = new ArrayList<HarmonicOscillator>();
		}
		harmonicFunctions.add(new HarmonicOscillator(a, TWO_PI * frequency, phi));
	}

	/**
	 * {@inheritDoc}
	 */
	public double value(double x) {
		double value = 0;
		if (polynomialFunction != null) {
			value += polynomialFunction.value(x);
		}
		if (harmonicFunctions != null) {
			final Iterator<HarmonicOscillator> iter = harmonicFunctions.iterator();
			while (iter.hasNext()) {
				final UnivariateFunction function = iter.next();
				value += function.value(x);
			}
		}
		return value;
	}

	/**
	 * Evaluate the function for an array of input values
	 * 
	 * @param x
	 *            the array of values at which to evaluate the function
	 * @return the array of values of the function
	 */
	public double[] value(double[] x) {
		final double[] values = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			values[i] = this.value(x[i]);
		}
		return values;
	}

	/**
	 * Parametric function where the input array contains the parameters of the
	 * time series function, ordered as follows:
	 * <ul>
	 * <li>Polynomial terms</li>
	 * <li>Harmonic terms (A, f, phi)</li>
	 * </ul>
	 * 
	 * The parameters do not contain the frequency term for harmonics of a
	 * fundamental frequency as these frequencies are not variable parameters
	 * but should always remain integer multiples of the fundamental after any
	 * fitting procedure.
	 */
	public static class Parametric implements ParametricUnivariateFunction {

		/** The degree of the polynomial */
		private int polynomialDegree;

		/** The list of harmonic terms for each fundamental frequency */
		private List<FrequencyHarmonicPair> frequencyHarmonicTerms;

		/**
		 * The total number of parameters (polynomila + harmonic) in the
		 * function
		 */
		private int numModelParameters;

		/**
		 * Construct a parametric time series function with only polynomial
		 * terms
		 * 
		 * @param polynomialDegree
		 *            the degree of the polynomial
		 * @throws MathIllegalArgumentException
		 */
		public Parametric(int polynomialDegree) throws MathIllegalArgumentException {
			if (polynomialDegree < 0) {
				throw new ModelSpecificationException(LocalizedFormats.NON_POSITIVE_POLYNOMIAL_DEGREE,
						polynomialDegree);
			}
			this.polynomialDegree = polynomialDegree;
			numModelParameters = polynomialDegree + 1;
		}

		/**
		 * Construct a parametric time series function with polynomial terms and
		 * one frequency
		 * 
		 * @param polynomialDegree
		 *            the degree of the polynomial
		 * @param frequency
		 *            the fundamental frequency
		 * @throws MathIllegalArgumentException
		 */
		public Parametric(int polynomialDegree, double frequency) throws MathIllegalArgumentException {
			this(polynomialDegree,Stream.of(new FrequencyHarmonicPair(frequency, IntStream.rangeClosed(1, 1).toArray())).collect(Collectors.toList()));
		}

		/**
		 * Construct a parametric time series function with polynomial terms and
		 * one frequency with harmonics
		 * 
		 * @param polynomialDegree
		 *            the degree of the polynomial
		 * @param frequency
		 *            the fundamental frequency
		 * @param numHarmonics
		 *            the number of harmonics
		 * @throws MathIllegalArgumentException
		 */
		public Parametric(int polynomialDegree, double frequency, int numHarmonics)
				throws MathIllegalArgumentException {
			this(polynomialDegree,
					Stream.of(new FrequencyHarmonicPair(frequency, IntStream.rangeClosed(1, numHarmonics).toArray())).collect(Collectors.toList()));
		}

		/**
		 * Construct a parametric time series function with polynomial terms and
		 * one frequency with harmonics
		 * 
		 * @param polynomialDegree
		 *            the degree of the polynomial
		 * @param frequency
		 *            the fundamental frequency
		 * @param harmonics
		 *            the list of harmonic terms
		 * @throws MathIllegalArgumentException
		 */
		public Parametric(int polynomialDegree, double frequency, int[] harmonics) throws MathIllegalArgumentException {
			this(polynomialDegree, Stream.of(new FrequencyHarmonicPair(frequency, harmonics)).collect(Collectors.toList()));
		}

		

//		/**
//		 * Construct a parametric time series function with polynomial terms and
//		 * many frequencies
//		 * 
//		 * @param polynomialDegree
//		 *            the degree of the polynomial
//		 * @param frequencies
//		 *            the fundamental frequencies
//		 * @param numHarmonics
//		 *            the number of harmonics for each fundamental frequency
//		 * @throws MathIllegalArgumentException
//		 */
//		public Parametric(int polynomialDegree, double[] frequencies, int[] numHarmonics)
//				throws MathIllegalArgumentException {
//			this(polynomialDegree, ModellingUtils.getFrequencyHarmonics(frequencies, numHarmonics));
//		}

	
		
	
		
		/**
		 * Construct a TimeSeriesFunction defining the number of polynomial
		 * coefficients, fundamental frequencies and harmonics
		 * 
		 * @param polynomialDegree
		 *            the degree of the polynomial
		 * @param frequencyHarmonics
		 *            the harmonic terms of each fundamental frequency
		 * @throws MathIllegalArgumentException
		 */
		public Parametric(int polynomialDegree, List<FrequencyHarmonicPair> frequencyHarmonics)
				throws MathIllegalArgumentException {
			if (polynomialDegree < 0) {
				throw new ModelSpecificationException(LocalizedFormats.NON_POSITIVE_POLYNOMIAL_DEGREE,
						polynomialDegree);
			}
			this.polynomialDegree = polynomialDegree;
			numModelParameters = polynomialDegree + 1;

			if (frequencyHarmonics != null) {
				// Check that the first harmonic of each frequency is the
				// fundamental, i.e 1
				for (FrequencyHarmonicPair frequencyHarmonic : frequencyHarmonics) {
					int[] harmonics = frequencyHarmonic.getHarmonics();
					if (harmonics.length < 1) {
						throw new RuntimeException(
								"FUNDAMENTAL_FIRST_HARMONIC_NOT_ONE");
					}
					if (harmonics[0] != 1) {
						throw new RuntimeException(
								"FUNDAMENTAL_FIRST_HARMONIC_NOT_ONE");
					}
					numModelParameters += 3 + 2 * (harmonics.length - 1);
				}
				this.frequencyHarmonicTerms = frequencyHarmonics;
			}
		}

		/**
		 * {@inheritDoc}
		 */
		@Override
		public double value(double x, double... params) throws NullArgumentException, DimensionMismatchException {

			validateParameters(params);

			// Evaluate the polynomial component
			int nPolynomialTerms = polynomialDegree + 1;
			double result = params[nPolynomialTerms - 1];
			for (int j = nPolynomialTerms - 2; j >= 0; j--) {
				result = x * result + params[j];
			}

			// Evaluate the harmonic components
			if (frequencyHarmonicTerms == null) {
				return result;
			}

			// Compute the harmonic component. The harmonics of the fundamental
			// must
			// always remain integer multiples of the fundamental
			int idx = polynomialDegree + 1;
			for (FrequencyHarmonicPair frequencyHarmonic : frequencyHarmonicTerms) {
				int[] harmonics = frequencyHarmonic.getHarmonics();
				double fundFreq = Double.NaN;
				for (int h = 0; h < harmonics.length; h++) {
					int harmonic = harmonics[h];
					boolean isFundamental = (harmonic == 1) ? true : false;
					final double amplitude = params[idx];
					fundFreq = (isFundamental) ? params[idx + 1] : fundFreq;
					double phi = (isFundamental) ? params[idx + 2] : params[idx + 1];
					result += amplitude * FastMath.cos(TWO_PI * x * harmonic * fundFreq + phi);
					idx = (isFundamental) ? idx + 3 : idx + 2;
				}
			}
			return result;
		}

		/**
		 * {@inheritDoc}
		 */
		@Override
		public double[] gradient(double x, double... params) {

			MathUtils.checkNotNull(params);

			// The number of terms in the model, taking into account that the
			// frequencies of the harmonics > 1 cannot vary. Minimum of one for
			// the constant
			// terms.
			if (!(numModelParameters == params.length)) {
				
				
				throw new RuntimeException("NUM_MODEL_TERMS_INCORRECT");
			}

			double[] gradient = new double[params.length];

			// Partial derivative w.r.t. polynomial terms
			double xn = 1.0;
			for (int j = 0; j < polynomialDegree + 1; j++) {
				gradient[j] = xn;
				xn *= x;
			}

			if (frequencyHarmonicTerms == null) {
				return gradient;
			}

			// Partial derivative w.r.t. amplitude, frequency, phase
			// Only the fundamental frequency terms should be optimized.
			// Harmonics must remain integer multiples of the fundamental
			int funFreqIdx = Integer.MIN_VALUE;
			int idx = polynomialDegree + 1;
			for (FrequencyHarmonicPair frequencyHarmonic : frequencyHarmonicTerms) {
				int[] harmonics = frequencyHarmonic.getHarmonics();
				double fundFreq = Double.NaN;

				for (int h = 0; h < harmonics.length; h++) {
					int harmonic = harmonics[h];
					boolean isFundamental = (harmonic == 1) ? true : false;
					final double amplitude = params[idx];
					fundFreq = (isFundamental) ? params[idx + 1] : fundFreq;
					double phi = (isFundamental) ? params[idx + 2] : params[idx + 1];

					gradient[idx++] = Math.cos(TWO_PI * harmonic * fundFreq * x + phi);

					// partial derivative w.r.t the frequency
					if (isFundamental) {
						funFreqIdx = idx;
						gradient[idx++] = -TWO_PI * x * amplitude * Math.sin(TWO_PI * fundFreq * x + phi);
					} else {
						// Harmonic frequency - add in contribution from
						// non-varying integer multiples of the fundamental
						gradient[funFreqIdx] += -TWO_PI * x * amplitude * harmonic
								* Math.sin(TWO_PI * harmonic * fundFreq * x + phi);
					}

					// partial derivative w.r.t the phase
					gradient[idx++] = -1.0 * amplitude * Math.sin(TWO_PI * harmonic * fundFreq * x + phi);
				}
			}
			return gradient;
		}

		/**
		 * Configure the polynomial degree of the parametric function
		 * 
		 * @return the polynimial degree
		 */
		public int getPolynomialDegree() {
			return polynomialDegree;
		}

		/**
		 * Get the number of parameters of the parametric function
		 * 
		 * @return the numModelParameters
		 */
		public final int getNumModelParameters() {
			return numModelParameters;
		}

		/**
		 * Get the harmonicTerms of the parametric function
		 * 
		 * @return the harmonicTerms
		 */
		public final List<FrequencyHarmonicPair> getHarmonicTerms() {
			return frequencyHarmonicTerms;
		}

		/**
		 * Validates parameters to ensure they are appropriate for the
		 * evaluation of the {@link #value(double,double[])} and
		 * {@link #gradient(double,double[])} methods.
		 *
		 * @param params
		 *            Values of the coeffricients of the model
		 * @throws NullArgumentException
		 *             if {@code param} is {@code null}.
		 * @throws DimensionMismatchException
		 *             if the size of {@code param} incorrect.
		 */
		private void validateParameters(double[] params) throws NullArgumentException, DimensionMismatchException {

			MathUtils.checkNotNull(params);

			if (params.length < 1) {
				throw new DimensionMismatchException(params.length, 1);
			}

			// A degree 0 polynomial term must be included
			if (polynomialDegree < 0) {
				throw new ModelSpecificationException(LocalizedFormats.NON_POSITIVE_POLYNOMIAL_DEGREE,
						polynomialDegree);
			}

			// The number of terms in the model, taking into account that the
			// frequencies of the
			// harmonics > 1 cannot vary. Minimum of one for the constant
			// terms.
			if (!(numModelParameters == params.length)) {
				throw new RuntimeException("NUM_MODEL_TERMS_INCORRECT");
			}
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public DerivativeStructure value(DerivativeStructure t) throws MathIllegalArgumentException {
		DerivativeStructure value = polynomialFunction.value(t);
		for (HarmonicOscillator harmonicFunction : harmonicFunctions) {
			value.add(harmonicFunction.value(t));
		}
		return value;
	}
}
