using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using System.IO;
using MathNetExtensions;

namespace UKF
{
    /// <summary>
    /// <para>Unscented transform estimate for the model y = Phi(x) + Nu, where</para> 
    /// <para>- x is a random variable with known mean and covariance,</para>
    /// <para>- Nu - noise with zero mean and known covariance</para>
    /// <para>Usage:</para>
    /// <para>- specify the parameters of the unscented transform in utProperty manually or by means of primitive optmization procedure: UTStaticEstimate.EstimateParameters</para>
    /// <para>- calculate the estimate: UTStaticEstimate.Estimate</para>
    /// <para>It should be noted, that the train and test sets may be different. 
    /// That is, the arrays of samples X = [x_0,...,x_N] and observations Y = [y_0,...,y_N] = [Phi(x_0) + Nu_0,...,Phi(x_N) + Nu_N] may vary 
    /// for the step of the unscented transform parameters optimization and the step of unscented transform estimate calculation.</para>
    /// </summary>
    public class UTStaticEstimate
    {
        public UTParams utParams; // parameters of the unscented transform
        UTDefinitionType utDefinitionType;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="type">The way to calculate the unscented transform parameters in the optimization procedure.</param>
        public UTStaticEstimate(UTDefinitionType type = UTDefinitionType.ImplicitAlpha)
        {
            utDefinitionType = type;
        }

        /// <summary>
        /// Provides the unscented transform estimate for the array X = [x_0,...,x_N] of random variable x samples 
        /// given the array of observations Y = [y_0,...,y_N], where y_i = Phi(x_i) + Nu_i. 
        /// Parameters of the unscented transform utParams should be initialized.
        /// </summary>
        /// <param name="Phi">Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu</param>
        /// <param name="X">Array of initial variable x samples</param>
        /// <param name="Y">Array of transformed variable y = Phi(x) + nu samples</param>
        /// <param name="mX">Mean of x</param>
        /// <param name="KX">Cov of x</param>
        /// <param name="KNu">Cov of the noize nu</param>
        /// <param name="mErr_UT">Returns: estimation error mean vector</param>
        /// <param name="KErr_UT">Returns: estimation error covariance marix</param>
        /// <param name="KErrTh_UT">Returns: estimation error theoretical covariance marix</param>
        /// <returns>Array of estimates \hat{X} = [\hat{x}_0,...,\hat{x}_N]</returns>
        public Vector<double>[] Estimate(Func<Vector<double>, Vector<double>> Phi, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu,
            out Vector<double> mErr_UT, out Matrix<double> KErr_UT, out Matrix<double> KErrTh_UT)
        {
            UnscentedTransform.Transform(x => Phi(x), mX, KX, KNu, utParams, out Vector<double> M_UT, out Matrix<double> Kxy_UT, out Matrix<double> Kyy_UT);
            Matrix<double> P_UT = Kxy_UT * ((Kyy_UT).PseudoInverse());
            KErrTh_UT = KX - P_UT * Kyy_UT * P_UT.Transpose();

            Vector<double>[] Xhat_UT = Y.Select(y => mX + P_UT * (y - M_UT)).ToArray();
            Vector<double>[] Err_UT = Xhat_UT.Subtract(X);
            mErr_UT = Err_UT.Average();
            KErr_UT = Exts.Cov(Err_UT, Err_UT);

            return Xhat_UT;
        }

        /// <summary>
        /// Calls the static unscented transform parameters "optimization" procedure and saves the result into the utParams property.
        /// The way how the random samples are transformed to the UT params is determined by the optimizationType property.
        /// </summary>
        /// <param name="N1">Number of samples on the step 1</param>
        /// <param name="N2">Number of samples on the step 2</param>
        /// <param name="Phi">Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat))  </param>
        /// <param name="X">Array of initial variable x samples</param>
        /// <param name="Y">Array of transformed variable y = Phi(x) + nu samples</param>
        /// <param name="mX">Mean of x</param>
        /// <param name="KX">Cov of x</param>
        /// <param name="KNu">Cov of the noize nu</param>
        /// <param name="outputFolder">The results are saved to this folder in file "UT_optimization_{type}.txt"</param>
        public void EstimateParameters(int N1, int N2, Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu, string outputFolder)
        {
            utParams = UTParmsOptimize(N1, N2, utDefinitionType, Phi, Crit, X, Y, mX, KX, KNu, outputFolder);
            //optimalParams = UTParmsOptimize(N1, N2, UTOptimizationType.ImplicitAlpha, Phi, Crit, X, Y, mX, KX, KNu);
            //optimalParams = UTParmsOptimize(N1, N2, UTOptimizationType.ImplicitAlphaBetaKappa, Phi, Crit, X, Y, mX, KX, KNu);
            //optimalParams = UTParmsOptimize(N1, N2, UTOptimizationType.Explicit, Phi, Crit, X, Y, mX, KX, KNu);

        }

        /// <summary>
        /// <para>Unscented transform parameters "optimization" procedure.</para>
        /// <para>- Step 1: generate N1 random samples, calculate the unscented transform estimates given the parameters determined by each random sample.</para>
        /// <para>- Step 2:choose the random sample with the best estimate quality criterion value and generate N2 random samples in closer area of this sample.</para>
        /// <para>- Step 3:again choose the random sample with the best estimate quality criterion value, save the samples ordered by the criterion value to the output file and return the best found unscented transform parameters.</para>
        /// <para>The UTOptimizationType type param determines the way how the random samples define the unscented tranform params (UTParams).</para>
        /// <para>- If type is UTOptimizationType.ImplicitAlpha, then random samples define alpha0 - scalar weight of the central points for both sample mean and cov: UTParams.SetUTParams(int, double);</para>
        /// <para>- If type is UTOptimizationType.ImplicitAlphaBetaKappa, then random samples are vectors of dim 3 and represent three parameters alpha, beta, kappa which are then transformed to the parameters of the inscented transform: UTParams.SetUTParams(int, double, double, double);</para> 
        /// <para>- If type is UTOptimizationType.Explicit, then random samples are vectors of dim 4 and explicitly define the unscented transform parameters: UTParams.SetUTParams(int, double, double, double, double). ///TODO it is not right to define the parameters of the unsctnted transform arbitraty, they have to be interdependent, so that the mean and cov would be transformed correctly.</para>
        /// </summary>
        /// <param name="N1">Number of samples on the step 1</param>
        /// <param name="N2">Number of samples on the step 2</param>
        /// <param name="type">The way how the random samples are transformed to the UT params</param>
        /// <param name="Phi">Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat))  </param>
        /// <param name="X">Array of initial variable x samples</param>
        /// <param name="Y">Array of transformed variable y = Phi(x) + nu samples</param>
        /// <param name="mX">Mean of x</param>
        /// <param name="KX">Cov of x</param>
        /// <param name="KNu">Cov of the noize nu</param>
        /// <param name="outputFolder">The results are saved to this folder in file "UT_optimization_{type}.txt"</param>
        /// <returns>Returns the parameters of the unscented transform with best estimation quality</returns>
        static UTParams UTParmsOptimize(int N1, int N2, UTDefinitionType type, Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu, string outputFolder)
        {
            return (new UTParams(1, 1));
            //int n = 0;
            //string filename = Path.Combine(outputFolder, "UT_optimization_{type}.txt");
            //switch (type)
            //{
            //    case UTDefinitionType.ImplicitAlpha: n = 1; filename = filename.Replace("{type}", "ImplicitAlpha"); break;
            //    case UTDefinitionType.ImplicitAlphaBetaKappa: n = 3; filename = filename.Replace("{type}", "ImplicitABK"); break;
            //    case UTDefinitionType.Explicit: n = 4; filename = filename.Replace("{type}", "Explicit"); break;
            //}

            //IContinuousDistribution[] distr = new IContinuousDistribution[n];
            //for (int i = 0; i < n; i++)
            //{
            //    distr[i] = new ContinuousUniform(-15, 15);
            //}

            //AsyncCalculatorPlanner acp = new AsyncCalculatorPlanner(N1, 10, () => CalculateSampleCriterion(Phi, Crit, distr, X, Y, mX, KX, KNu));

            //List<double[]> results1 = acp.DoCalculate();

            //double min1 = results1.Where(i => !double.IsNaN(i[0])).Min(e => e[0]);
            //double[] best1 = results1.First(e => e[0] == min1);

            //for (int i = 0; i < n; i++)
            //{
            //    distr[i] = new Normal(best1[i + 5], Math.Sqrt(0.05));
            //}

            //acp = new AsyncCalculatorPlanner(N2, 10, () => CalculateSampleCriterion(Phi, Crit, distr, X, Y, mX, KX, KNu));
            //List<double[]> results2 = acp.DoCalculate();

            //double min2 = results2.Where(i => !double.IsNaN(i[0])).Min(e => e[0]);
            //double[] best2 = results2.First(e => e[0] == min2);


            //using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(filename))
            //{
            //    NumberFormatInfo provider;
            //    provider = new NumberFormatInfo
            //    {
            //        NumberDecimalSeparator = "."
            //    };

            //    var results = results1.Concat(results2).Where(i => !double.IsNaN(i[0]) && i[0] < double.MaxValue);
            //    foreach (var e in results.OrderBy(e => e[0]))
            //    {
            //        outputfile.WriteLine(String.Join(",", e.Select(s => string.Format(provider, "{0}", s))));
            //    }

            //}
            //return new UTParams(mX.Count, best2.Skip(5).Take(n).ToArray());
        }

        /// <summary>
        /// Generates a random sample for the unscented transform parameters and calculates the criterion value for the unscented transform estimate.
        /// The size of the distribution parameter determines the unscented transform parameters definition method (UTParams)
        /// </summary>
        /// <param name="Phi">Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat))  </param>
        /// <param name="distribution">Array of distributions to generate random unscented transform parameters</param>
        /// <param name="X">Array of initial variable x samples</param>
        /// <param name="Y">Array of transformed variable y = Phi(x) + nu samples</param>
        /// <param name="mX">Mean of x</param>
        /// <param name="KX">Cov of x</param>
        /// <param name="KNu">Cov of the noize nu</param>
        /// <returns></returns>
        public static double[] CalculateSampleCriterion(Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, IContinuousDistribution[] distribution, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {
            double[] s = distribution.Select(d => d.Sample()).ToArray();
            int L = mX.Count;

            UTParams p = new UTParams(mX.Count, s);
            double crit = CalculateCriterionValue(Phi, Crit, p, X, Y, mX, KX, KNu);
            return (new double[] { crit, p.Lambda, p.Wc[0], p.Wm[0], p.Wm[1] }).Concat(s).ToArray();
        }

        /// <summary>
        /// Calculates the criterion value for the estimate given the particular unscented transform parameters
        /// </summary>
        /// <param name="Phi">Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat))  </param>
        /// <param name="p">Unscented transform parameters</param>
        /// <param name="X">Array of initial variable x samples</param>
        /// <param name="Y">Array of transformed variable y = Phi(x) + nu samples</param>
        /// <param name="mX">Mean of x</param>
        /// <param name="KX">Cov of x</param>
        /// <param name="KNu">Cov of the noize nu</param>
        /// <returns></returns>
        public static double CalculateCriterionValue(Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, UTParams p, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {
            double crit = 0;
            try
            {
                UnscentedTransform.Transform(x => Phi(x), mX, KX, KNu, p, out Vector<double> M_UT, out Matrix<double> Kxy_UT, out Matrix<double> Kyy_UT);
                Matrix<double> P_UT = Kxy_UT * ((Kyy_UT).PseudoInverse());
                Matrix<double> KErrTh_UT = KX - P_UT * Kyy_UT * P_UT.Transpose();

                Vector<double>[] Xhat_UT = Y.Select(y => mX + P_UT * (y - M_UT)).ToArray();
                Vector<double>[] Err_UT = Xhat_UT.Subtract(X);
                Vector<double> mErr_UT = Err_UT.Average();
                Matrix<double> KErr_UT = Exts.Cov(Err_UT, Err_UT);
                crit = Crit(KErr_UT);//.Trace();
            }
            catch { crit = double.MaxValue; }
            return crit;
        }

    }

}

    
