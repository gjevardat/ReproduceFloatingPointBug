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

import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.regression.ModelSpecificationException;
import org.apache.commons.math3.util.MathUtils;


public final class NegativeAmplitudeValidator implements ParameterValidator {


	/** The parametric time series function to validate */
	private TimeSeriesFunction.Parametric function;

	/**
	 * Construct a NegativeAmplitudeValidator
	 * 
	 * @param function the parametric function to validate
	 */
	public NegativeAmplitudeValidator(TimeSeriesFunction.Parametric function) {
		this.function = function;
	}

	/**
	 * {@inheritDoc}
	 */
	public RealVector validate(RealVector params) {

		if (params.getDimension() != function.getNumModelParameters()) {
			throw new RuntimeException("Incorrect Number of terms");
		}
		
		if (function.getHarmonicTerms() == null) {
			return params.copy();
		}

		// Loop over each fundamental frequency and associated harmonic terms - check the appropriate terms
		int idx = function.getPolynomialDegree()+1;
		RealVector p = params.copy();
		for (FrequencyHarmonicPair frequencyHarmonics : function.getHarmonicTerms()) {
			int[] harmonics = frequencyHarmonics.getHarmonics();
			double fundFreq = Double.NaN;

			for (int h=0; h<harmonics.length; h++) {
				int harmonic = harmonics[h];
				boolean isFundamental = (harmonic ==1) ? true : false;
				final double amplitude = p.getEntry(idx);
				fundFreq = (isFundamental) ?  p.getEntry(idx + 1) : fundFreq;
				double phi = (isFundamental) ? p.getEntry(idx + 2) :  p.getEntry(idx + 1);
				
				// If the amplitude has gone negative, fix
				if (amplitude <0) { 
					p.setEntry(idx, -1.0*amplitude);
					phi = MathUtils.normalizeAngle((phi+ Math.PI), Math.PI);
					if (isFundamental) {
						p.setEntry(idx+2, phi);
					} else {
						p.setEntry(idx+1, phi);
					}
				}
				idx = (isFundamental) ? idx+3 : idx+2;
			}
		}
		return p;
	}
}
