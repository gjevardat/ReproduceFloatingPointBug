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

import java.util.Collection;
import java.util.stream.DoubleStream;

import org.apache.commons.lang3.Validate;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.function.Inverse;
import org.apache.commons.math3.analysis.function.Sqrt;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.NotANumberException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.fitting.AbstractCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresFactory;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem.Evaluation;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.SimpleVectorValueChecker;
import org.apache.commons.math3.stat.regression.ModelSpecificationException;


/**
 * Create an instance of the fitter using the factory method to fit a collection
 * of observed data points. A vector of parameters will be returned Curve
 * fitting for the univariate time series model function
 * {@code TimeSeriesFunction}.
 * <p>
 * Fit a collection of observed data points to the general time series model
 * given by y(x) = a_0 + a_1 * x + a_2 * x^2 + ... + a_{n-1} x^(n-1) +
 * sin[2*PI*f*(t-t0)] + cos[2*PI*f*(t-t0)] + .... sin[(l-1) * 2*PI*l*f*(t-t0)] +
 * cos[(l-1) * 2*PI*l*f*(t-t0)]
 * <p>
 * where: x are the covariates, y are the responses a_i are the regression
 * coefficients (to be determined) and the f_i(x) are the regressors evaluated
 * in the covariates.
 * <p>
 * A vector of best fit parameters will be returned
 * 
 * <p>
 * 
 * {@code withReferenceTime} Subtract a reference time from the observation
 * times. The default is 0. If no reference time is passed, the observations
 * times will be used in the model as is.
 * <p>
 * {@code withInitialGuess} allows the user to pass a known start point. If no
 * initial start point values are provided, the ParameterGuesser will used
 * linear algebra will be used to determine a start point and the model
 * parameters.
 * <p>
 * {@code withOptimizer} A non-linear optimizer is used to fit the data to the
 * model If an initial guess {@link InitialGuess} is provided this is used as
 * the start point otherwise linear least squares is used to obtain a starting
 * point.
 * 
 * {@code withMaxIter} The maximum number of iterations {@link MaxIter} The
 * maximum number of iterations
 * 
 * {@code withMaxEval} The maximum number of evaluations {@link MaxEval} The
 * maximum number of evaluations
 * 
 * {@code withCovergenceChecker} Optional. A specific user defined convergence
 * checker can be set. If not, a default convergence checker is used. Only
 * relevant with a least squares optimizer. The initial data needed is
 * <ul>
 * <li>The polynomial degree of the model</li>
 * <li>The number of fundamental frequencies in the model</li>
 * <li>The harmonic terms associated with each fundamental frequency</li>
 * </ul>
 * 
 * 
 * @since 16.0.0
 * @author Leanne Guy &lt;leanne.guy@unige.ch&gt;
 * @version $Id$
 * 
 *          Requirements implemented by this class
 * @req
 */
public final class TimeSeriesModelFitter extends AbstractCurveFitter {


	/** The specific TimeSeriesFunction to fit */
	//private ParametricUnivariateFunction function;
	private TimeSeriesFunction.Parametric function;

	/** The input initial guess */
	private double[] initialGuess = null;

	/**
	 * The start point taken either from the input {@code initialGuess} or
	 * computed by linear decomposition
	 */
	private double[] startPoint = null;

	/** The LeastSquaresOptimizer */
	private LeastSquaresOptimizer optimizer = null;

	/** Maximum number of iterations of the optimization algorithm. */
	private int maxIter;

	/** Maximum number of evaluations of the optimization algorithm. */
	private int maxEvaluations;

	/** Return a new optimizer set up to fit a Gaussian curve to the observed points. */
	private ConvergenceChecker<Evaluation> convergenceChecker = null;

	/** Weighted fit */
	private boolean weighting;


	/** The Optimum resulting from a non-linear optimization */
	private Optimum optimum;

	/** The covariance singularity threshold  */
	private double covarianceSingularityThreshold;

	/** The weights on the data to be fit */
	private RealVector weights;

	/** The TheoreticalValuesFunction */
	private TheoreticalValuesFunction model;

	/** The default relative threshold for the convergence checker for least squares optimization */
	private static final double REL_THRESH = 1.0e-18;

	/** The default absolute threshold for the convergence checker for least squares optimization */
	private static final double ABS_THRESH = 1.0e-18;

	/** The default covariance singularity threshold */
	private static final Double DEFAULT_COVARIANCE_SINGULARITY_THRESHOLD = 1e-018;

	/**
	 * Prevent direct instantiation
	 */
	private TimeSeriesModelFitter() {
	}

	/**
	 * Construct a TimeSeriesModelFitter
	 * 
	 * @param builder
	 */
	private TimeSeriesModelFitter(Builder builder) {

		this.function = builder.function;
		this.initialGuess = (builder.initialGuess == null) ? null : builder.initialGuess;
		this.maxIter = (builder.maxIter == null) ? Integer.MAX_VALUE : builder.maxIter;
		this.maxEvaluations = (builder.maxEvaluations == null) ? Integer.MAX_VALUE : builder.maxEvaluations;
		this.optimizer = builder.optimizer;

		// Default weighting is false
		this.weighting = (builder.weighting == null) ? false : builder.weighting;

		// Create a default convergence checker if none provided
		this.convergenceChecker = (builder.convergenceChecker == null)
				? LeastSquaresFactory.evaluationChecker(new SimpleVectorValueChecker(REL_THRESH, ABS_THRESH))
						: builder.convergenceChecker;

		// Set the default covarianceSingularityThreshold is none provided
		this.covarianceSingularityThreshold = (builder.covarianceSingularityThreshold == null)
				? DEFAULT_COVARIANCE_SINGULARITY_THRESHOLD : builder.covarianceSingularityThreshold;

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public double[] fit(Collection<WeightedObservedPoint> points) {
		optimum = getOptimizer().optimize(getProblem(points)); 
		return optimum.getPoint().toArray();
	}

	/**
	 * Builder for the {@code TimeSeriesModelFitter}
	 * 
	 * @since 20.0.0
	 * @author Leanne Guy &lt;leanne.guy@unige.ch&gt;
	 * @version $Id$
	 */
	public static class Builder {

		private TimeSeriesFunction.Parametric function;
		private double[] initialGuess;
		private Integer maxIter;
		private Integer maxEvaluations;
		private Boolean weighting;
		private LeastSquaresOptimizer optimizer;
		private ConvergenceChecker<Evaluation> convergenceChecker;
		private Double covarianceSingularityThreshold; 

		/**
		 * Set the function to fit
		 * 
		 * @param function
		 * @return the TimeSeriesModelFitter builder
		 */
		public Builder withFunction(final TimeSeriesFunction.Parametric function) {
			this.function = function;
			return this;
		}

		/**
		 * Set the
		 * 
		 * @param initialGuess
		 * @return the TimeSeriesModelFitter builder
		 */
		public Builder withInitialGuess(final double[] initialGuess) {
			this.initialGuess = initialGuess;
			return this;
		}

		/**
		 * Set the maximum number of iterations
		 * 
		 * @param maxIter
		 * @return the TimeSeriesModelFitter builder
		 */
		public Builder withMaxIter(final Integer maxIter) {
			this.maxIter = maxIter;
			return this;
		}

		/**
		 * Set the maximum number of evaluations
		 * 
		 * @param maxEvaluations
		 * @return the TimeSeriesModelFitter builder
		 */
		public Builder withMaxEvaluations(final Integer maxEvaluations) {
			this.maxEvaluations = maxEvaluations;
			return this;
		}

		/**
		 * Set the maximum number of evaluations
		 * 
		 * @param weighting
		 * @return the TimeSeriesModelFitter builder
		 */
		public Builder withWeighting(final Boolean weighting) {
			this.weighting = weighting;
			return this;
		}

		/**
		 * Set the potimizer
		 * 
		 * @param optimizer
		 *            the LeastSquaresOptimizer
		 * @return the TimeSeriesModelFitter builder
		 */
		public Builder withOptimizer(final LeastSquaresOptimizer optimizer) {
			this.optimizer = optimizer;
			return this;
		}

		/**
		 * Set the convercenge checker
		 * 
		 * @param convergenceChecker
		 * @return the TimeSeriesModelFitter builder
		 */
		public Builder withConvergenceChecker(final ConvergenceChecker<Evaluation> convergenceChecker) {
			this.convergenceChecker = convergenceChecker;
			return this;
		}

		/**
		 * Set the covarianceSingularityThreshold
		 * @param covarianceSingularityThreshold
		 * @return the TimeSeriesModelFitter builder
		 */
		public Builder withCovarianceSingularityThreshold(double covarianceSingularityThreshold) {
			this.covarianceSingularityThreshold = covarianceSingularityThreshold;
			return this;

		}

		/**
		 * Build the TimeSeriesModelFitter
		 * 
		 * @return the TimeSeriesModelFitter
		 */
		public TimeSeriesModelFitter build() {
			validate();
			return new TimeSeriesModelFitter(this);
		}

		private void validate() {

			Validate.notNull(function, "Model function may not be null");

			// Start point may be null. If not it must have the same dimensions
			// as the function
			if (initialGuess != null) {
				//	int nModelParams = ((TimeSeriesFunction.Parametric) function).getNumModelParameters();
				int nModelParams = function.getNumModelParameters();
				Validate.isTrue(initialGuess.length == nModelParams,
						"initial guess does not have the correct number of values " + initialGuess.length + "/"
								+ nModelParams);
			}
			if (covarianceSingularityThreshold != null) {
				Validate.isTrue(covarianceSingularityThreshold > 0, "covarianceSingularityThreshold should be >0");
			}
			if (maxIter != null) {
				Validate.isTrue(maxIter > 0, "maxIter should be >0");
			}
			if (maxEvaluations != null) {
				Validate.isTrue(maxEvaluations > 0, "maxEvaluations should be set");
			}
		}
	}

	/**
	 * Get the builder
	 * 
	 * @return a new Builder
	 */
	public static Builder getBuilder() {
		return new TimeSeriesModelFitter.Builder();
	}

	/**
	 * Return the optimised point with any post-processing validators applied 
	 * @return the 
	 * @throws MathIllegalArgumentException 
	 */
	public RealVector getPoint() throws MathIllegalArgumentException {

		RealVector result = optimum.getPoint();
		
		// Check for NaN values - problem with convergence
        boolean isNaN = DoubleStream.of(result.toArray()).anyMatch(x-> Double.isNaN(x));
        if (isNaN) {
        	throw new NotANumberException();
        }

		// Check for negative amplitudes and correct  
		if (function.getHarmonicTerms() != null) {
			result = new NegativeAmplitudeValidator(function).validate(optimum.getPoint());	
		}
		
		return result;
	}

	/**
	 * Return the errors on the optimised point
	 * @return the errors on the optimised point
	 * @throws MathIllegalArgumentException 
	 */
	public RealVector getPointErrors() throws MathIllegalArgumentException {
		RealVector result = optimum.getSigma(covarianceSingularityThreshold);
		
		// Check for NaN values - indicates a problem with convergence
        boolean isNaN = DoubleStream.of(result.toArray()).anyMatch(x-> Double.isNaN(x));
        if (isNaN) {
        	throw new NotANumberException();
        }
		return result;
	}

	/**
	 * Get the startPoint. If a start point has been passed using
	 * {@code InitialGuess} this initialGuess value will be returned. If no
	 * initial guess is passed, the initial guess computed using linear algebra
	 * will be returned.
	 * 
	 * @return the startPoint
	 */
	public final double[] getStartPoint() {
		return startPoint;
	}

	

	/**
	 * Get the optimum
	 * @return the optimum
	 */
	public final Optimum getOptimum() {
		return optimum;
	}

	/**
	 * Get the raw residuals
	 * Note that optimum.getResiduals returns the weighted residuals not the raw residuals. 
	 * To get the raw residuals, the weighted residuals must be multiplied by the square root 
	 * of the weight
	 * @return the residuals
	 */
	public final RealVector getResiduals() {
		RealVector invWeights = weights.copy().mapToSelf(new Sqrt());
		return optimum.getResiduals().ebeDivide(invWeights) ;
	}

	/**
	 * Set the 
	 * @return the 
	 */
	public final RealVector getResidualErrors() {
		
		final RealMatrix identityM = MatrixUtils.createRealIdentityMatrix(weights.getDimension());
		RealMatrix hatMatrix = computeHatMatrix(optimum.getPoint().toArray());
		final RealMatrix identityMatrixMinusHat = identityM.subtract(hatMatrix);
		final RealMatrix sigma = new DiagonalMatrix(weights.copy().mapToSelf(new Inverse()).toArray());
		final RealMatrix sigmaIdentityMatrixMinusHat = sigma.multiply(identityMatrixMinusHat);
		final RealMatrix residualErrorsMatrix = identityMatrixMinusHat.transpose()
				.multiply(sigmaIdentityMatrixMinusHat);
	
		double[] resErrors =getDiagonal(residualErrorsMatrix);
		for (int i = 0; i < resErrors.length; i++) {
			resErrors[i] = Math.sqrt(resErrors[i]);
		}
		return MatrixUtils.createRealVector(resErrors);
	}
	
	public static double[] getDiagonal(RealMatrix matrix) {
		if (!matrix.isSquare()) {
			return null;
		}
		final double[] diag = new double[matrix.getColumnDimension()];
		for (int i=0; i<matrix.getColumnDimension(); i++) {
			diag[i] = matrix.getEntry(i, i);
		}
		return diag;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected LeastSquaresOptimizer getOptimizer() {
		return (optimizer == null) ? super.getOptimizer() : this.optimizer;
	}

	/**
	 * {@inheritDoc}
	 */
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> observations) {

		// Prepare least-squares problem.
		final int len = observations.size();
		final double[] times = new double[len];
		final double[] target = new double[len];
		final double[] wgts = new double[len]; 

		// Only because linear models uses errors - TODO
		final double[] errors = new double[len]; 
		int i = 0;
		for (WeightedObservedPoint obs : observations) {
			times[i] = obs.getX();
			target[i] = obs.getY();
			wgts[i] = obs.getWeight();

			// TODO harmonise weights/errors - not neceessarily INV_SQ and look at rounding
			errors[i] = Math.sqrt(1.0 / obs.getWeight());
			// errors[i] = (double) Math.round(Math.sqrt(1.0 / obs.getWeight())
			// * 1000d) / 1000d ;
			++i;
		}

		// The weights are stored to compute the raw residuals after the fit
		weights = MatrixUtils.createRealVector(wgts);

		// The curve fitter function
		model = new AbstractCurveFitter.TheoreticalValuesFunction(function, observations);

		// If an initial guess is provided it is used to seed the non-linear
		// optimization.
		// If an initial guess is not provided, it is used and linear modelling
		// using the exact model definition is used to determine the start point
		// Extract the data to compute the linear model
		startPoint = initialGuess;

		// Use a parameter guesser ?
		// new ParameterGuesser(observations).guess();

		return new LeastSquaresBuilder().maxEvaluations(maxEvaluations).maxIterations(maxIter).start(startPoint)
				.target(target).weight(new DiagonalMatrix(wgts)).checker(convergenceChecker)
				.model(model.getModelFunction(), model.getModelFunctionJacobian()).build();
	}

	

	/**
	 * Compute the hat matrix from the model jacobian 
	 * H = J(J^T Σ^−1J)^−1J^T Σ^−1,
	 * 
	 * TODO In Utils? 
	 * @param params
	 *            the value of the model parameters at which to
	 * @return the hat matrix
	 */
	private RealMatrix computeHatMatrix(final double[] params) {
	
		// A model needs to have been computed 
		if (model==null) {
			throw new ModelSpecificationException(LocalizedFormats.NULL_NOT_ALLOWED, model);
		}
		
		// TODO Can this be used ? 
		//optimum.getJacobian();
		
		final MultivariateMatrixFunction  jacobianFunction = model.getModelFunctionJacobian();
		final RealMatrix modelJacobian =MatrixUtils.createRealMatrix(jacobianFunction.value(params));
		final RealMatrix modelJacobianT = modelJacobian.transpose();

		//final RealMatrix sigmaInverse = model.getFit().getOmegaInverse();
		final RealMatrix sigmaInverse = new DiagonalMatrix(weights.toArray());
		final RealMatrix modelJacobianTransposeSigmaInverse = modelJacobianT.
				multiply(sigmaInverse.multiply(modelJacobian));

		final DecompositionSolver solver = new QRDecomposition(modelJacobianTransposeSigmaInverse, 1E-024).getSolver();
		if (!solver.isNonSingular()) {
			throw new SingularMatrixException();
		}
		final RealMatrix m = solver.getInverse().multiply(modelJacobianT.multiply(sigmaInverse));
		return modelJacobian.multiply(m);
	}
}