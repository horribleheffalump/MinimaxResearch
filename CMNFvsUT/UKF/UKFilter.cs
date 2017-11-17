using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using MathNetExtensions;
using System.IO;
using MathNet.Numerics.Optimization;

namespace UKF
{
    /// <summary>
    /// <para>Unscented Kalman filter for the model x_{t+1} = Phi(x_t) + W_t, y[t] = Psi(x_t) + Nu_t, where</para> 
    /// <para>- x_t is unobservable state with disturbances W_t ~ (0, R_w) in dynamics,</para>
    /// <para>- y_t - observations with noise Nu_t ~ (0, R_nu).</para>
    /// <para>Usage:</para>
    /// <para>- specify the parameters of the unscented transform on forecast and correction phases in utParamsForecast and utParamsCorrection properties
    /// manually or by means of the optmization procedure: UKFilter.EstimateParameters</para>
    /// <para>- calculate the estimate step by step with UKFilter.Step</para>
    /// <para>It should be noted, that the train and test trajectory bundles may be different. 
    /// That is, the arrays of discrete vector models may vary 
    /// for the step of the unscented transform parameters optimization and the step of unscented transform filter calculation.</para>
    /// </summary>
    public class UKFilter
    {
        public UTParams utParamsForecast; // parameters of the unscented transform on forecast phase
        public UTParams utParamsCorrection; // parameters of the unscented transform on correction phase
        UTDefinitionType utDefinitionType;
        OptimizationMethod optimizationMethod;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="type">Unscented transform parameters definition type</param>
        /// <param name="method">Unscented transform parameters optimization method</param>
        public UKFilter(UTDefinitionType type = UTDefinitionType.ImplicitAlpha, OptimizationMethod method = OptimizationMethod.NelderMeed)
        {
            utDefinitionType = type;
            optimizationMethod = method;
        }

        /// Calls the static unscented transform parameters optimization procedure and saves the result into the utParamsForecast and utParamsCorrection properties.
        /// The way how to define the UT params is determined by the optimizationType property.
        /// The optimization method is determined by the optimizationMethod property.
        /// </summary>
        /// <param name="Phi">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y[t] = Psi(x_t) + Nu_t</param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="T">The upper bound of the observation interval</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xhat0">Initial condition</param>
        /// <param name="DX0Hat">Initial condition covariance</param>
        /// <param name="outputFolder">The results are saved to this folder in file "UT_optimization_{type}.txt"</param>
        public void EstimateParameters(Func<int, Vector<double>, Vector<double>> Phi,
                                             Func<int, Vector<double>, Vector<double>> Psi,
                                             Matrix<double> Rw,
                                             Matrix<double> Rnu,
                                             Func<Matrix<double>, double> Crit,
                                             int T,
                                             DiscreteVectorModel[] models,
                                             Vector<double> xhat0,
                                             Matrix<double> DX0Hat,
                                             string outputFolder)
        {
            (_, utParamsForecast, utParamsCorrection) = UTParmsOptimize(optimizationMethod, utDefinitionType, Phi, Psi, Rw, Rnu, Crit, T, models, xhat0, DX0Hat, outputFolder);
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
        /// <param name="Phi">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y[t] = Psi(x_t) + Nu_t</param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="T">The upper bound of the observation interval</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xhat0">Initial condition</param>
        /// <param name="DX0Hat">Initial condition covariance</param>
        /// <param name="outputFolder">The results are saved to this folder in file "UT_optimization_{type}.txt"</param>
        static (double, UTParams, UTParams) UTParmsOptimize(OptimizationMethod method, UTDefinitionType type,
                                             Func<int, Vector<double>, Vector<double>> Phi,
                                             Func<int, Vector<double>, Vector<double>> Psi,
                                             Matrix<double> Rw,
                                             Matrix<double> Rnu,
                                             Func<Matrix<double>, double> Crit,
                                             int T,
                                             DiscreteVectorModel[] models,
                                             Vector<double> xhat0,
                                             Matrix<double> DX0Hat,
                                             string outputFolder)
        {

            int n = 0;
            string filename = string.IsNullOrWhiteSpace(outputFolder) ? null : Path.Combine(outputFolder, "UT_optimization_{type}.txt");
            Vector<double> lowerBound;
            Vector<double> upperBound;
            Vector<double> initialGuess;
            switch (type)
            {
                case UTDefinitionType.ImplicitAlpha:
                    n = 1;
                    lowerBound = Exts.Vector(1 - 2 / xhat0.Count);
                    upperBound = Exts.Vector(1);
                    initialGuess = Exts.Vector(0.5);
                    filename = filename.Replace("{type}", "ImplicitAlpha");
                    break;
                case UTDefinitionType.ImplicitAlphaBetaKappa:
                    n = 3;
                    lowerBound = Exts.Vector(0, 0, 0);
                    upperBound = Exts.Vector(5, 5, 5);
                    initialGuess = Exts.Vector(1, 2, 1);
                    filename = filename.Replace("{type}", "ImplicitABK"); break;
                case UTDefinitionType.Explicit:
                    n = 4;
                    lowerBound = Exts.Vector(-10, -10, -10, -10);
                    upperBound = Exts.Vector(10, 10, 10, 10);
                    initialGuess = Exts.Vector((new UTParams(xhat0.Count, 1, 2, 1)).Params);
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
                    var OptimumRandom = RandomOptimizer.Minimize((x) => CalculateSampleCriterion(Phi, Psi, Rw, Rnu, Crit, x, T, models, xhat0, DX0Hat), Exts.Stack(lowerBound, lowerBound), Exts.Stack(upperBound, upperBound), 100, 100, filename);
                    min = OptimumRandom.min;
                    argmin = OptimumRandom.argmin;
                    break;
                case OptimizationMethod.NelderMeed:
                    NelderMeadSimplex optimizer = new NelderMeadSimplex(1e-3, 100);
                    var objective = ObjectiveFunction.Value((x) => CalculateSampleCriterion(Phi, Psi, Rw, Rnu, Crit, x, T, models, xhat0, DX0Hat));
                    var optimumNM = optimizer.FindMinimum(objective, Exts.Stack(initialGuess, initialGuess));
                    min = optimumNM.FunctionInfoAtMinimum.Value;
                    argmin = optimumNM.MinimizingPoint;
                    break;

            }
            return (min, new UTParams(xhat0.Count, argmin.Take(n).ToArray()), new UTParams(xhat0.Count, argmin.Skip(n).Take(n).ToArray()));
        }

        /// <summary>
        /// Wrapper for CalculateCriterionValue function. 
        /// The sample vector is transformed to a couple of the unscented transform parameters (the way it is done depends on the length of the sample vector, see UTParams and its constructors for details). 
        /// Then the criterion value given the provided unscented transform parameters for the forecast and the correction phases is calculated. 
        /// </summary>
        /// <param name="Phi">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y[t] = Psi(x_t) + Nu_t</param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="P">Sample vector to be transformed to a couple of the unscented transfrom parameters</param>
        /// <param name="T">The upper bound of the observation interval</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xhat0">Initial condition</param>
        /// <param name="DX0Hat">Initial condition covariance</param>
        /// <returns>The criterion value for the unscented transfrom parameters obtained from the sample vactor</returns>
        public static double CalculateSampleCriterion(Func<int, Vector<double>, Vector<double>> Phi,
                                     Func<int, Vector<double>, Vector<double>> Psi,
                                     Matrix<double> Rw,
                                     Matrix<double> Rnu,
                                     Func<Matrix<double>, double> Crit,
                                     Vector<double> P,
                                     int T,
                                     DiscreteVectorModel[] models,
                                     Vector<double> xhat0,
                                     Matrix<double> DX0Hat
                                     )
        {
            int n = P.Count();
            if (n % 2 != 0)
                new ArgumentException("Distribution number must be even", "distribution");

            int L = xhat0.Count;

            double[] s1 = P.Take(n / 2).ToArray();
            UTParams p1 = new UTParams(L, s1);

            double[] s2 = P.Skip(n / 2).ToArray();
            UTParams p2 = new UTParams(L, s2);

            double crit = CalculateCriterionValue(Phi, Psi, Rw, Rnu, Crit, p1, p2, T, models, xhat0, DX0Hat);

            return crit;
        }

        /// <summary>
        /// Calculates the criterion value for the estimate given the particular unscented transform parameters
        /// </summary>
        /// <param name="Phi">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y[t] = Psi(x_t) + Nu_t</param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="p1">Unscented transfrom parameters for the forecast phase</param>
        /// <param name="p2">Unscented transfrom parameters for the correction phase</param>
        /// <param name="T">The upper bound of the observation interval</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xhat0">Initial condition</param>
        /// <param name="DX0Hat">Initial condition covariance</param>
        /// <returns>The criterion value for the particular unscented transform parameters</returns>
        public static double CalculateCriterionValue(Func<int, Vector<double>, Vector<double>> Phi,
                                                     Func<int, Vector<double>, Vector<double>> Psi,
                                                     Matrix<double> Rw,
                                                     Matrix<double> Rnu,
                                                     Func<Matrix<double>, double> Crit,
                                                     UTParams p1,
                                                     UTParams p2,
                                                     int T,
                                                     DiscreteVectorModel[] models,
                                                     Vector<double> xhat0,
                                                     Matrix<double> DX0Hat
                                                     )
        {
            double crit = 0;
            try
            {
                int N = models.Count();

                Vector<double>[] xHatU = models.Select(x => xhat0).ToArray();
                Matrix<double>[] PHatU = models.Select(x => DX0Hat).ToArray();
                //Vector<double> PHatU = Vector<double>.Build.Dense(N, DX0Hat);



                for (int t = 0; t < T; t++)
                {
                    for (int i = 0; i < N; i++)
                    {
                        Step(Phi, Psi, Rw, Rnu, p1, p2, t, models[i].Trajectory[t][1], xHatU[i], PHatU[i], out Vector<double> xHatU_i, out Matrix<double> PHatU_i);
                        xHatU[i] = xHatU_i;
                        PHatU[i] = PHatU_i;
                    }

                    Vector<double>[] states = models.Select(x => (x.Trajectory[t][0])).ToArray();
                    Matrix<double> errorUPow2 = Exts.Cov(states.Subtract(xHatU), states.Subtract(xHatU));

                    crit = Crit(errorUPow2);
                }
            }
            catch { crit = double.MaxValue; }
            return crit;
        }

        /// <summary>
        /// Performs a step of Unscented Kalman Filter given the particular unscented transform parameters
        /// for forecast and correction phases
        /// </summary>
        /// <param name="Phi">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y[t] = Psi(x_t) + Nu_t</param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="p1">Unscented transfrom parameters for the forecast phase</param>
        /// <param name="p2">Unscented transfrom parameters for the correction phase</param>
        /// <param name="t">Current step time instant</param>
        /// <param name="y">Observations on the current step</param>
        /// <param name="xHat_">Estimate on the previous step</param>
        /// <param name="P_">Approximated previous step error covariance</param>        
        /// <param name="xHat">Returns: current step estimate</param>
        /// <param name="P">Returns: approximated current step error covariance</param>        
        public static void Step(Func<int, Vector<double>, Vector<double>> Phi,
                                Func<int, Vector<double>, Vector<double>> Psi,
                                Matrix<double> Rw,
                                Matrix<double> Rnu,
                                UTParams p1,
                                UTParams p2,
                                int t,
                                Vector<double> y,
                                Vector<double> xHat_,
                                Matrix<double> P_,
                                out Vector<double> xHat,
                                out Matrix<double> PHat)
        {
            UnscentedTransform.Transform(x => Phi(t, x), xHat_, P_, Rw, p1, out Vector<double> Xtilde, out _, out Matrix<double> Ptilde);
            UnscentedTransform.Transform(x => Psi(t, x), Xtilde, Ptilde, Rnu, p2, out Vector<double> Ytilde, out Matrix<double> PXY, out Matrix<double> PYtilde);
            Matrix<double> K = PXY * PYtilde.Inverse();
            xHat = Xtilde + K * (y - Ytilde);
            PHat = Ptilde - K * PYtilde * K.Transpose();
        }

        /// <summary>
        /// Performs a step of Unscented Kalman Filter with fixed transformation parameters
        /// for forecast and correction phases (utParamsForecast and utParamsCorrection must be initialized)
        /// </summary>
        /// <param name="Phi">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y[t] = Psi(x_t) + Nu_t</param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="t">Current step time instant</param>
        /// <param name="y">Observations on the current step</param>
        /// <param name="xHat_">Estimate on the previous step</param>
        /// <param name="P_">Approximated previous step error covariance</param>        
        /// <param name="xHat">Returns: current step estimate</param>
        /// <param name="P">Returns: approximated current step error covariance</param>        
        public void Step(Func<int, Vector<double>, Vector<double>> Phi,
                                Func<int, Vector<double>, Vector<double>> Psi,
                                Matrix<double> Rw,
                                Matrix<double> Rnu,
                                int t,
                                Vector<double> y,
                                Vector<double> xHat_,
                                Matrix<double> P_,
                                out Vector<double> xHat,
                                out Matrix<double> PHat)
        {
            Step(Phi, Psi, Rw, Rnu, utParamsForecast, utParamsCorrection, t, y, xHat_, P_, out xHat, out PHat);
        }
    }
}

