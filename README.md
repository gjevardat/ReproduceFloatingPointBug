# ReproduceFloatingPointBug

To reproduce the bug
Run the unit test named TimeSeriesModelFitterBugTest

To see the difference of behavior just switch the java.version property from 1.8 to 17 in the pom.xml.
* On JDK 1.8 it shall converge and output (syserr) the value
*   [11.979873441619297, 2.133235451846515E-4, 0.10993866745791972, 0.0015814425063299668, 4.931738739425212]

* On JDK 17 it shall throw the following exception
*   org.apache.commons.math3.exception.ConvergenceException:
	at org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer.optimize(LevenbergMarquardtOptimizer.java:536)
	at gaia.dpcg.TimeSeriesModelFitter.fit(TimeSeriesModelFitter.java:210)
	at gaia.dpcg.TimeSeriesModelFitterBugTest.reproduceRegression(TimeSeriesModelFitterBugTest.java:109)
