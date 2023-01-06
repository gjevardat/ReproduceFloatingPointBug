# ReproduceFloatingPointBug

To reproduce the bug
Run the unit test named TimeSeriesModelFitterBugTest

To see the difference of behavior just switch the java.version property from 1.8 to 17 in the pom.xml.
* On JDK 1.8 it shall converge and output (syserr) the value: [11.979873441619297, 2.133235451846515E-4, 0.10993866745791972, 0.0015814425063299668, 4.931738739425212]

* On JDK 17 it shall throw the following exception:
 org.apache.commons.math3.exception.ConvergenceException
	at org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer.optimize(LevenbergMarquardtOptimizer.java:536)
	at gaia.dpcg.TimeSeriesModelFitter.fit(TimeSeriesModelFitter.java:210)
	at gaia.dpcg.TimeSeriesModelFitterBugTest.reproduceRegression(TimeSeriesModelFitterBugTest.java:109)



I can pinpoint the location in Apache math commons of why the Java 17 version throws exception: In the class LevenbergMarquardtOptimizer of Apache Math commmons 3.6.1 (line 453) 

  double actRed = -1.0;
      if (0.1 * currentCost < previousCost) {
          double r = currentCost / previousCost;
          actRed = 1.0 - r * r;
      }
At iteration 17 actRed becomes 0 because r==1 in Java 17 while it is slightly positive in Java 8 Because of the value being 0, later (line 533) the non convergenceException is thrown

 if (FastMath.abs(actRed) <= TWO_EPS &&
     preRed <= TWO_EPS && ratio <= 2.0) {
         throw new ConvergenceException(LocalizedFormats.TOO_SMALL_COST_RELATIVE_TOLERANCE,
                                               costRelativeTolerance);
            } 
