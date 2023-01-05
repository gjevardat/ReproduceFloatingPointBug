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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresFactory;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem.Evaluation;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.SimpleVectorValueChecker;

/**
 * Minimal reproducible example of the floating point issues between Java 8 and Java 17
 */
public class TimeSeriesModelFitterBugTest2  {


	@org.junit.jupiter.api.Test
	public void reproduceRegression() throws Exception {

		double[] obstimes = new double[] { 1769.19447887209822, 1793.11335157792223, 1793.18734866346449, 1824.84317446474029, 1824.9171491695015, 1862.64812786623224, 1862.82424848524033, 1862.89830169656966, 1886.82817182598365, 1886.90217904917267, 1930.56263613548481, 1930.63661741753253,
				1960.55749143200819, 1960.6315376658888, 1990.79828079462595, 2029.83639229992332, 2047.76510303999771, 2047.83911783365966, 2090.70738114505957, 2090.78138516485524, 2106.78380068309389, 2145.41454883987853, 2145.48858125141805, 2171.72794244868419, 2201.64019651627223,
				2201.71420294044719, 2242.62177696932258, 2242.69579125066502, 2266.37672127201404, 2266.45075082078665, 2308.18529315827436, 2340.10080028618131, 2340.17477593992226, 2367.59446360805805, 2367.66847315820314, 2408.12655151605168, 2424.30609532823837, 2424.38010415350982,
				2467.07194494119585, 2484.74901952488108, 2521.70939279792128, 2550.51942914891561, 2578.76164520528891, 2621.97752425096633, 2622.05149727853131, 2622.22769049269436, 2646.48207256349133, 2646.55606933698073 };
		double[] vals = new double[] { 12.0194979779388191, 12.0198270700936938, 12.0283468968662355, 11.8157307317799258, 11.8292824887879711, 11.8364549278854163, 11.8243874344468267, 11.8246456132521747, 11.9491042429398995, 11.9447542832251408, 11.8922851272789565, 11.8907090967527331,
				11.8217218832161972, 11.8292590886469497, 11.7951313891598666, 11.9691138027728634, 11.927750517204279, 11.9458558038541867, 11.7670720729569265, 11.7655584688580426, 11.8110680135031085, 12.0513566259304827, 12.0391523130277385, 11.9564206736219347, 12.0115885825541877,
				12.0233580691422297, 11.9186735465738955, 11.9257010783989568, 12.0957273964296093, 12.0741823643380677, 12.2502406453301216, 12.0663815198101343, 12.0587377784384024, 12.1459262000534327, 12.170732078122418, 12.2699171274285188, 12.2221578844271157, 12.2576334497918591,
				11.8368712431067795, 11.9058313920358803, 11.9271965273927076, 12.0537857209779133, 11.9047418265925415, 11.9554862499311607, 11.9553430520426289, 11.9471229366267888, 12.016369081699704, 12.0249000892867954 };
		double[] errs = new double[] { 0.00578021267792471621, 0.00310395835363509905, 0.00506232696862514778, 0.00266974593429538439, 0.00289622753097086441, 0.00384728068036813229, 0.00352722521407869503, 0.00217235509739079286, 0.00143385015485857407, 0.00240811057817686552,
				0.00168286614634127123, 0.00250305523405236268, 0.00183143473758006175, 0.00515516789412332789, 0.00594126983681109166, 0.00384881424177748593, 0.00123286206150668488, 0.00313741007042465155, 0.00308568961087013396, 0.00396212336382095742, 0.00199436633918759904,
				0.00342499092300945981, 0.00583258661553120051, 0.00681011827913153948, 0.00241390955578568017, 0.00244310496433679551, 0.00301603699291613835, 0.00520996150719486677, 0.00403799681136299184, 0.00342220110235065805, 0.00371660033844624978, 0.00348266209268642372,
				0.005134341946080756, 0.00251798863494797579, 0.00981085575457805135, 0.00211295762645843368, 0.00221172367572451244, 0.00484976864499936551, 0.0051466669325059165, 0.00367297646878710126, 0.00383912278176944382, 0.00501193292409431168, 0.00205133435451890743, 0.00254625889711689165,
				0.00300766411345966848, 0.00282910501115183483, 0.00391228159439286228, 0.00373353551727234481 };


		// 0.00122495113
		// Model definition Final model: Model [x0 + f0 ([1, 2, 3, 4])]
		int polynomialDegree = 1;
		double periodSearchFrequency = 0.00122495113;
		int[] harmonics = new int[] { 1};

		// Taken from the linear model results in the taker. Both tests should
		// start with the same initial guess for comparison. Weighted, ,
		// reference time = TIMESERIES_MEAN
		double[] initialGuess = new double[] { 11.980784702177685, 1.7941443371450404E-5, 0.13143925467546005, 0.00122495113, 4.837618003254914};
		double referenceTime = 2191.2714506308716;
		List<FrequencyHarmonicPair> fHarmonics = new LinkedList<FrequencyHarmonicPair>();
		fHarmonics.add(0, new FrequencyHarmonicPair(periodSearchFrequency, harmonics));
		final TimeSeriesFunction.Parametric function = new TimeSeriesFunction.Parametric(polynomialDegree, fHarmonics);
		
		
		final double relThreshold = 1.0e-20;
		final double absThreshold = 1.0e-20;
		// Set up the default convergence threshold
		ConvergenceChecker<Evaluation> convergenceChecker = LeastSquaresFactory
				.evaluationChecker(new SimpleVectorValueChecker(relThreshold, absThreshold));

		final double initialStepBoundFactor = 100.0;
		final double costRelativeToTolerance = 1.0e-30;
		final double parRelativeToTolerance = 1.0e-10;
		final double orthoTolerance = 1.0e-30;
		final double qrRankingThreshold = 1.0e-100;

		LeastSquaresOptimizer optimizer =  new LevenbergMarquardtOptimizer()
			.withCostRelativeTolerance(costRelativeToTolerance)
			.withInitialStepBoundFactor(initialStepBoundFactor).withOrthoTolerance(orthoTolerance)
			.withParameterRelativeTolerance(parRelativeToTolerance)
			.withRankingThreshold(qrRankingThreshold)			;
		
			
		final TimeSeriesModelFitter fitterNLM = new TimeSeriesModelFitter.Builder().withFunction(function)
		.withInitialGuess(initialGuess).withWeighting(false)
		.withConvergenceChecker(convergenceChecker).withOptimizer(optimizer)
		.withMaxEvaluations(1000).withMaxIter(1000).build();
		
		
		double[] resultNLM = fitterNLM.fit(getWeightedPoints(referenceTime,obstimes,vals,errs));
				
		
		System.err.println(Arrays.toString(resultNLM));
	}
	
	/**
	 * Build list of Apache WeightedObservedPoint from input data
	 * @param referenceTime
	 * @param obstimes
	 * @param vals
	 * @param errs
	 * @param weightingStrategy
	 * @return
	 */
	private List<WeightedObservedPoint>  getWeightedPoints(double referenceTime, double[] obstimes,double[] vals,double[] errs ){
		List<WeightedObservedPoint> list = new ArrayList<>() ;
		for(int i=0; i<obstimes.length; i++) {
			WeightedObservedPoint p = new WeightedObservedPoint(1.0, obstimes[i] - referenceTime, vals[i]);
			list.add(p);
		}
		return list;
	}
}

