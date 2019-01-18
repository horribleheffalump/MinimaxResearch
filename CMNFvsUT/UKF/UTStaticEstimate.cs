using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using System.IO;
using MathNetExtensions;
using MathNet.Numerics.Optimization;

namespace UKF
{
    /// <summary>
    /// <para>Unscented transform estimate for the model y = Phi(x) + Nu, where</para> 
    /// <para>- x is a random variable with known mean and covariance,</para>
    /// <para>- Nu - noise with zero mean and known covariance</para>
    /// <para>Usage:</para>
    /// <para>- specify the parameters of the unscented transform in utProperty manually or by means of the optmization procedure UTStaticEstimate.EstimateParameters</para>
    /// <para>- calculate the estimate: UTStaticEstimate.Estimate</para>
    /// <para>It should be noted, that the train and test sets may be different. 
    /// That is, the arrays of samples X = [x_0,...,x_N] and observations Y = [y_0,...,y_N] = [Phi(x_0) + Nu_0,...,Phi(x_N) + Nu_N] may vary 
    /// for the step of the unscented transform parameters optimization and the step of unscented transform estimate calculation.</para>
    /// </summary>
    public class UTStaticEstimate
    {
        public UTParams utParams; // parameters of the unscented transform
        UTDefinitionType utDefinitionType;
        OptimizationMethod optimizationMethod;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="type">Unscented transform parameters definition type</param>
        /// <param name="method">Unscented transform parameters optimization method</param>
        public UTStaticEstimate(UTDefinitionType type = UTDefinitionType.ImplicitAlpha, OptimizationMethod method = OptimizationMethod.NelderMeed)
        {
            utDefinitionType = type;
            optimizationMethod = method;
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
        /// Calls the static unscented transform parameters optimization procedure and saves the result into the utParams property.
        /// The way how to define the UT params is determined by the optimizationType property.
        /// The optimization method is determined by the optimizationMethod property.
        /// </summary>
        /// <param name="Phi">Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat))  </param>
        /// <param name="X">Array of initial variable x samples</param>
        /// <param name="Y">Array of transformed variable y = Phi(x) + nu samples</param>
        /// <param name="mX">Mean of x</param>
        /// <param name="KX">Cov of x</param>
        /// <param name="KNu">Cov of the noize nu</param>
        /// <param name="outputFolder">The results are saved to this folder in file "UT_optimization_{type}.txt"</param>
        public void EstimateParameters(Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu, string outputFolder)
        {
            (_, utParams) = UTParmsOptimize(optimizationMethod, utDefinitionType, Phi, Crit, X, Y, mX, KX, KNu, outputFolder);

        }

        /// <summary>
        /// <para>Unscented transform parameters optimization procedure.</para>
        /// <para>The OptimizationMethod param determines the optimization method:</para>
        /// <para>- OptimizationMethod.RandomShoot - parameters are randomly sampled and the best sample is chosen as optimal;
        /// <para>- OptimizationMethod.NelderMeed - parameters are optimized with non-gradient Nelder-Meed method.</para>
        /// <para>The UTOptimizationType type param determines the relation between the optimized variable and the unscented tranform params (see UTParams and its constructors for details). </para>
        /// <para>- If type is UTOptimizationType.ImplicitAlpha, then the optimized variable is saclar [alpha0];</para>
        /// <para>- If type is UTOptimizationType.ImplicitAlphaBetaKappa, then optimized variable is a vector [alpha, beta, kappa];</para> 
        /// <para>- If type is UTOptimizationType.Explicit, then then optimized variable is a vector [lambda, wm0, wc0, wi]. ///TODO it is not correct to define the parameters of the unsctnted transform arbitraty, they have to be interdependent, so that the mean and cov would be transformed correctly.</para>
        /// </summary>
        /// <param name="method">Unscented transform parameters optimization method</param>
        /// <param name="type">Unscented transform parameters definition type</param>
        /// <param name="Phi">Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat))  </param>
        /// <param name="X">Array of initial variable x samples</param>
        /// <param name="Y">Array of transformed variable y = Phi(x) + nu samples</param>
        /// <param name="mX">Mean of x</param>
        /// <param name="KX">Cov of x</param>
        /// <param name="KNu">Cov of the noize nu</param>
        /// <param name="outputFolder">If needed, the results or the random sampling (OptimizationMethod.RandomShoot) are saved to this folder in file "UT_optimization_{type}.txt"</param>
        /// <returns>Returns touple (the best criteron value, the parameters of the unscented transform with best estimation quality)</returns>
        static (double, UTParams) UTParmsOptimize(OptimizationMethod method, UTDefinitionType type,
                                     Func<Vector<double>, Vector<double>> Phi, 
                                     Func<Matrix<double>, double> Crit, 
                                     Vector<double>[] X, 
                                     Vector<double>[] Y, 
                                     Vector<double> mX, 
                                     Matrix<double> KX, 
                                     Matrix<double> KNu, 
                                     string outputFolder = null
            )
        {

            int n;
            string filename = string.IsNullOrWhiteSpace(outputFolder) ? null : Path.Combine(outputFolder, "UT_optimization_{type}.txt");
            Vector<double> lowerBound;
            Vector<double> upperBound;
            Vector<double> initialGuess;
            switch (type)
            {
                case UTDefinitionType.ImplicitAlpha:
                    n = 1;
                    lowerBound = Exts.Vector(1 - 2 / mX.Count);
                    upperBound = Exts.Vector(1);
                    initialGuess = Exts.Vector(0.5);
                    filename = filename.Replace("{type}", "ImplicitAlpha");
                    break;
                case UTDefinitionType.ImplicitAlphaBetaKappa:
                    n = 3;
                    lowerBound = Exts.Vector(0, 0, 3.0 - mX.Count - 2.0);
                    upperBound = Exts.Vector(5, 5, 3.0 - mX.Count + 2.0);
                    initialGuess = Exts.Vector(0.5, 2.0, 3.0 - mX.Count);
                    filename = filename.Replace("{type}", "ImplicitABK"); break;
                case UTDefinitionType.Explicit:
                    n = 4;
                    lowerBound = Exts.Vector(-10, -10, -10, -10);
                    upperBound = Exts.Vector(10, 10, 10, 10);
                    initialGuess = Exts.Vector((new UTParams(mX.Count, 0.5, 2.0, 3.0 - mX.Count)).Params);
                    filename = filename.Replace("{type}", "Explicit"); break;
                default:
                    n = 0;
                    lowerBound = null;
                    upperBound = null;
                    initialGuess = null;
                    break;
            }

            double min = double.MaxValue;
            Vector<double> argmin = initialGuess;

            switch (method)
            {
                case OptimizationMethod.RandomShoot:
                    var OptimumRandom = RandomOptimizer.Minimize((p) => CalculateSampleCriterion(Phi, Crit, p, X, Y, mX, KX, KNu), lowerBound, upperBound, 1000, 1000, filename);
                    min = OptimumRandom.min;
                    argmin = OptimumRandom.argmin;
                    break;
                case OptimizationMethod.NelderMeed:
                    NelderMeadSimplex optimizer = new NelderMeadSimplex(1e-3, 100);
                    var objective = ObjectiveFunction.Value((p) => CalculateSampleCriterion(Phi, Crit, p, X, Y, mX, KX, KNu));
                    var optimumNM = optimizer.FindMinimum(objective, initialGuess);
                    min = optimumNM.FunctionInfoAtMinimum.Value;
                    argmin = optimumNM.MinimizingPoint;
                    break;

            }
            return (min, new UTParams(mX.Count, argmin.AsArray()));
        }

        /// <summary>
        /// Wrapper for CalculateCriterionValue function. 
        /// The sample vector is transformed to the form unscented transform parameters (the way it is done depends on the length of the sample vector, see UTParams and its constructors for details). 
        /// Then the criterion value given the provided unscented transform parameters is calculated. 
        /// </summary>
        /// <param name="Phi">Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat))  </param>
        /// <param name="P">Sample vector to be transformed ro unscented transfrom parameters</param>
        /// <param name="X">Array of initial variable x samples</param>
        /// <param name="Y">Array of transformed variable y = Phi(x) + nu samples</param>
        /// <param name="mX">Mean of x</param>
        /// <param name="KX">Cov of x</param>
        /// <param name="KNu">Cov of the noize nu</param>
        /// <returns>The criterion value for the unscented transfrom parameters obtained from the sample vactor</returns>
        public static double CalculateSampleCriterion(Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, Vector<double>P, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {
            int L = mX.Count;
            UTParams p = new UTParams(mX.Count, P.AsArray());
            double crit = CalculateCriterionValue(Phi, Crit, p, X, Y, mX, KX, KNu);
            return crit;
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
        /// <returns>The criterion value for the particular unscented transform parameters</returns>
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

    
