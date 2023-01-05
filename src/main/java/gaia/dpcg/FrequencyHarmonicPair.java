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

import org.apache.commons.math3.util.Pair;

/**
 * This class holds a fundamental frequency term and and its corresponding harmonic terms in 
 * a model.
 * 
 * @since 20.0.0
 * @author Leanne Guy &lt;leanne.guy@unige.ch&gt;
 * @version $Id$
 */
public final class FrequencyHarmonicPair extends Pair<Double, int[]>  implements Serializable {

    /** Serializable UID. */
	private static final long serialVersionUID = 20L;
	
	/**
     * Builds a frequency and harmonic terms pair.
     *
     * @param frequency The frequency term
     * @param harmonics The harmonics to mode with the frequency
	 */
	public FrequencyHarmonicPair(final Double frequency, 
								 final int[] harmonics) {
		super(frequency, harmonics);
	}

    /**
     * Gets the frequency value.
     *
     * @return a copy of the stored point.
     */
    public double getFrequency() {
        final Double p = getKey();
        return p == null ? null : p.doubleValue();
    }

    /**
     * Gets the harmonic terms associated this frequency
     *
     * @return a copy of the stored value of the objective function.
     */
    public int[] getHarmonics() {
        return super.getValue();
    }
}
