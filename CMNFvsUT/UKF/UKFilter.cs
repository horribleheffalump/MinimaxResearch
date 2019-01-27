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
using System.Threading.Tasks;

namespace UKF
{
    /// <summary>
    /// <para>Unscented Kalman filter for a model x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t, y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t, W_t ~ (M_w, R_w), Nu_t ~ (M_nu, R_nu)</para>
    /// <para>Usage:</para>
    /// <para>- specify the parameters of the unscented transform on forecast and correction phases in utParamsForecast and utParamsCorrection properties
    /// manually or by means of the optmization procedure: UKFilter.EstimateParameters</para>
    /// <para>- calculate the estimate step by step with UKFilter.Step</para>
    /// <para>The train and test trajectory bundles may be different. 
    /// That is, the arrays of discrete vector models may vary 
    /// for the step of the unscented transform parameters optimization and the step of unscented transform filter calculation.</para>
    /// </summary>
    public partial class UKFilter
    {
        public UTParams utParamsForecast; // parameters of the unscented transform on forecast phase defined for the whole trajectory
        public UTParams utParamsCorrection; // parameters of the unscented transform on correction phase defined for the whole trajectory
        public UTParams[] utParamsForecastStepwise; // parameters of the unscented transform on forecast phase defined separately for each step
        public UTParams[] utParamsCorrectionStepwise; // parameters of the unscented transform on correction phase defined separately for each step
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

        /// <summary>
        /// Calls the static unscented transform parameters optimization procedure and saves the result into the utParamsForecast and utParamsCorrection properties.
        /// The way how to define the UT params is determined by the optimizationType property.
        /// The optimization method is determined by the optimizationMethod property.
        /// </summary>
        /// <param name="Phi1">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t</param>
        /// <param name="Phi2">Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi1">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Psi2">Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Mw">Mean of the noise in the dynamics equation </param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Mnu">Mean of the noise in the obseration equation </param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="T">The upper bound of the observation interval</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xhat0">Initial condition</param>
        /// <param name="DX0Hat">Initial condition covariance</param>
        /// <param name="outputFolder">The results are saved to this folder in file "UT_optimization_{type}.txt"</param>
        public void EstimateParameters(Func<int, Vector<double>, Vector<double>> Phi1,
                                        Func<int, Vector<double>, Matrix<double>> Phi2,
                                        Func<int, Vector<double>, Vector<double>> Psi1,
                                        Func<int, Vector<double>, Matrix<double>> Psi2,
                                        Vector<double> Mw,
                                        Matrix<double> Rw,
                                        Vector<double> Mnu,
                                        Matrix<double> Rnu,
                                        Func<Matrix<double>, double> Crit,
                                        int T,
                                        DiscreteVectorModel[] models,
                                        Vector<double> xhat0,
                                        Matrix<double> DX0Hat,
                                        string outputFolder)
        {
            (_, utParamsForecast, utParamsCorrection) = UTParmsOptimize(optimizationMethod, utDefinitionType, Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, Crit, T, models, xhat0, DX0Hat, outputFolder);
            //(_, var utParamsForecast2, var utParamsCorrection2) = UTParmsOptimize(optimizationMethod, utDefinitionType, Phi1, Psi1, Rw, Rnu, Crit, T, models, xhat0, DX0Hat, outputFolder);
        }

        /// <summary>
        /// Calls the static unscented transform parameters stepwise optimization procedure and saves the result into the utParamsForecast and utParamsCorrection properties.
        /// The way how to define the UT params is determined by the optimizationType property.
        /// The optimization method is determined by the optimizationMethod property.
        /// </summary>
        /// <param name="Phi1">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t</param>
        /// <param name="Phi2">Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi1">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Psi2">Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Mw">Mean of the noise in the dynamics equation </param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Mnu">Mean of the noise in the obseration equation </param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="T">The upper bound of the observation interval</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xhat0">Initial condition</param>
        /// <param name="DX0Hat">Initial condition covariance</param>
        /// <param name="outputFolder">The results are saved to this folder in file "UT_optimization_{type}.txt"</param>
        public void EstimateParametersStepwise(Func<int, Vector<double>, Vector<double>> Phi1,
                                        Func<int, Vector<double>, Matrix<double>> Phi2,
                                        Func<int, Vector<double>, Vector<double>> Psi1,
                                        Func<int, Vector<double>, Matrix<double>> Psi2,
                                        Vector<double> Mw,
                                        Matrix<double> Rw,
                                        Vector<double> Mnu,
                                        Matrix<double> Rnu,
                                        Func<Matrix<double>, double> Crit,
                                        int T,
                                        DiscreteVectorModel[] models,
                                        Vector<double> xhat0,
                                        Matrix<double> DX0Hat,
                                        string outputFolder)
        {
            (_, utParamsForecastStepwise, utParamsCorrectionStepwise) = UTParmsOptimizeStepwise(optimizationMethod, utDefinitionType, Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, Crit, T, models, xhat0, DX0Hat, outputFolder);
            //(_, var utParamsForecast2, var utParamsCorrection2) = UTParmsOptimize(optimizationMethod, utDefinitionType, Phi1, Psi1, Rw, Rnu, Crit, T, models, xhat0, DX0Hat, outputFolder);
        }


        /// <summary>
        /// Defines the number of dimensions, lower and upper bound and the initial guess and the file name to store results for the Unscented transform parameters optimization procedures.
        /// </summary>
        /// <param name="type">Unscented transform parameters definition type</param>
        /// <param name="xhat0">Initial condition</param>
        /// <param name="filename">filename name template to substitute the optimization type name"</param>
        /// <returns></returns>
        static (int, Vector<double>, Vector<double>, Vector<double>, string) DefineOptimizationParameters(
                                                         UTDefinitionType type,
                                                         Vector<double> xhat0,
                                                         string filename
            )
        {
            int n;
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
                    lowerBound = Exts.Vector(0, 0, 3.0 - xhat0.Count - 2.0);
                    upperBound = Exts.Vector(5, 5, 3.0 - xhat0.Count + 2.0);
                    initialGuess = Exts.Vector(0.5, 2.0, 3.0 - xhat0.Count);
                    filename = filename.Replace("{type}", "ImplicitABK"); break;
                case UTDefinitionType.Explicit:
                    n = 4;
                    lowerBound = Exts.Vector(-10, -10, -10, -10);
                    upperBound = Exts.Vector(10, 10, 10, 10);
                    initialGuess = Exts.Vector((new UTParams(xhat0.Count, 0.5, 2.0, 3.0 - xhat0.Count)).Params);
                    filename = filename.Replace("{type}", "Explicit"); break;
                default:
                    n = 0;
                    lowerBound = null;
                    upperBound = null;
                    initialGuess = null;
                    break;
            }
            return (n, lowerBound, upperBound, initialGuess, filename);

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
        /// <param name="Phi1">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t</param>
        /// <param name="Phi2">Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi1">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Psi2">Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Mw">Mean of the noise in the dynamics equation </param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Mnu">Mean of the noise in the obseration equation </param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="T">The upper bound of the observation interval</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xhat0">Initial condition</param>
        /// <param name="DX0Hat">Initial condition covariance</param>
        /// <param name="outputFolder">The results are saved to this folder in file "UT_optimization_{type}.txt"</param>
        static (double, UTParams, UTParams) UTParmsOptimize(OptimizationMethod method, UTDefinitionType type,
                                                            Func<int, Vector<double>, Vector<double>> Phi1,
                                                            Func<int, Vector<double>, Matrix<double>> Phi2,
                                                            Func<int, Vector<double>, Vector<double>> Psi1,
                                                            Func<int, Vector<double>, Matrix<double>> Psi2,
                                                            Vector<double> Mw,
                                                            Matrix<double> Rw,
                                                            Vector<double> Mnu,
                                                            Matrix<double> Rnu,
                                                            Func<Matrix<double>, double> Crit,
                                                            int T,
                                                            DiscreteVectorModel[] models,
                                                            Vector<double> xhat0,
                                                            Matrix<double> DX0Hat,
                                                            string outputFolder)
        {

            (int n, Vector<double> lowerBound, Vector<double> upperBound, Vector<double> initialGuess, string filename) = DefineOptimizationParameters(type, xhat0, string.IsNullOrWhiteSpace(outputFolder) ? null : Path.Combine(outputFolder, "UT_optimization_{type}.txt"));
            double min = double.MaxValue;
            Vector<double> argmin = Exts.Stack(initialGuess, initialGuess);

            switch (method)
            {
                case OptimizationMethod.RandomShoot:
                    var OptimumRandom = RandomOptimizer.Minimize((x) => CalculateSampleCriterion(Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, Crit, x, T, models, xhat0, DX0Hat), Exts.Stack(lowerBound, lowerBound), Exts.Stack(upperBound, upperBound), 100, 100, filename);
                    min = OptimumRandom.min;
                    argmin = OptimumRandom.argmin;
                    break;
                case OptimizationMethod.NelderMeed:
                    NelderMeadSimplex optimizer = new NelderMeadSimplex(1e-3, 100);
                    var objective = ObjectiveFunction.Value((x) => CalculateSampleCriterion(Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, Crit, x, T, models, xhat0, DX0Hat));
                    try
                    {
                        var optimumNM = optimizer.FindMinimum(objective, Exts.Stack(initialGuess, initialGuess));
                        min = optimumNM.FunctionInfoAtMinimum.Value;
                        argmin = optimumNM.MinimizingPoint;
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine($"Optimizer faild, using the initail guess ({e.Message})");
                        argmin = Exts.Stack(initialGuess, initialGuess);
                    }
                    break;
                default: // no optimization by default
                    break;


            }
            return (min, new UTParams(xhat0.Count, argmin.Take(n).ToArray()), new UTParams(xhat0.Count, argmin.Skip(n).Take(n).ToArray()));
        }

        /// <summary>
        /// <para>Unscented transform parameters stepwize optimization procedure.</para>
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
        /// <param name="Phi1">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t</param>
        /// <param name="Phi2">Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi1">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Psi2">Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Mw">Mean of the noise in the dynamics equation </param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Mnu">Mean of the noise in the obseration equation </param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="T">The upper bound of the observation interval</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xhat0">Initial condition</param>
        /// <param name="DX0Hat">Initial condition covariance</param>
        /// <param name="outputFolder">The results are saved to this folder in file "UT_optimization_{type}.txt"</param>
        static (double, UTParams[], UTParams[]) UTParmsOptimizeStepwise(OptimizationMethod method, UTDefinitionType type,
                                                            Func<int, Vector<double>, Vector<double>> Phi1,
                                                            Func<int, Vector<double>, Matrix<double>> Phi2,
                                                            Func<int, Vector<double>, Vector<double>> Psi1,
                                                            Func<int, Vector<double>, Matrix<double>> Psi2,
                                                            Vector<double> Mw,
                                                            Matrix<double> Rw,
                                                            Vector<double> Mnu,
                                                            Matrix<double> Rnu,
                                                            Func<Matrix<double>, double> Crit,
                                                            int T,
                                                            DiscreteVectorModel[] models,
                                                            Vector<double> xhat0,
                                                            Matrix<double> DX0Hat,
                                                            string outputFolder)
        {
            UTParams[] pForecast = new UTParams[T];
            UTParams[] pCorrect = new UTParams[T];

            (int n, Vector<double> lowerBound, Vector<double> upperBound, Vector<double> initialGuess, string filename) = DefineOptimizationParameters(type, xhat0, string.IsNullOrWhiteSpace(outputFolder) ? null : Path.Combine(outputFolder, "UT_stepwise_ptimization_{type}.txt"));

            Vector<double>[] xHatU = models.Select(x => xhat0).ToArray();
            Matrix<double>[] PHatU = models.Select(x => DX0Hat).ToArray();

            double min = double.MaxValue;
            Console.WriteLine($"UKF estimate parameters start");
            DateTime start = DateTime.Now;

            for (int t = 1; t < T; t++)
            //Parallel.For(0, T, new ParallelOptions() { MaxDegreeOfParallelism = System.Environment.ProcessorCount }, t =>
            {
                DateTime startiteration = DateTime.Now;
                min = double.MaxValue;
                Vector<double> argmin = initialGuess;

                switch (method)
                {
                    case OptimizationMethod.RandomShoot:
                        var OptimumRandom = RandomOptimizer.Minimize((x) => CalculateSampleStepwiseCriterion(Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, Crit, x, t, models, xHatU, PHatU), Exts.Stack(lowerBound, lowerBound), Exts.Stack(upperBound, upperBound), 100, 100, filename);
                        min = OptimumRandom.min;
                        argmin = OptimumRandom.argmin;
                        break;
                    case OptimizationMethod.NelderMeed:
                        NelderMeadSimplex optimizer = new NelderMeadSimplex(1e-3, 100);
                        var objective = ObjectiveFunction.Value((x) => CalculateSampleStepwiseCriterion(Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, Crit, x, t, models, xHatU, PHatU));
                        try
                        {
                            var optimumNM = optimizer.FindMinimum(objective, Exts.Stack(initialGuess, initialGuess));
                            min = optimumNM.FunctionInfoAtMinimum.Value;
                            argmin = optimumNM.MinimizingPoint;
                        }
                        catch (Exception e)
                        {
                            Console.WriteLine($"Optimizer faild, using the initail guess ({e.Message})");
                            argmin = Exts.Stack(initialGuess, initialGuess);
                        }
                        break;
                }
                pForecast[t] = new UTParams(xhat0.Count, argmin.Take(n).ToArray());
                pCorrect[t] = new UTParams(xhat0.Count, argmin.Skip(n).Take(n).ToArray());
                for (int i = 0; i < models.Count(); i++)
                {
                    (xHatU[i], PHatU[i]) = Step(Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, pForecast[t], pCorrect[t], t, models[i].Trajectory[t][1], xHatU[i], PHatU[i]);
                }
                Console.WriteLine($"UKF estimate parameters for t={t}, done in {(DateTime.Now - startiteration).ToString(@"hh\:mm\:ss\.fff")}");
            }
            //    });
            DateTime finish = DateTime.Now;
            Console.WriteLine($"UKF estimate parameters finished in {(finish - start).ToString(@"hh\:mm\:ss\.fff")}");
            return (min, pForecast, pCorrect);
        }

        /// <summary>
        /// The sample vector is transformed to a couple of the unscented transform parameters (the way it is done depends on the length of the sample vector, see UTParams and its constructors for details). 
        /// </summary>
        /// <param name="P">Sample vector to be transformed to a couple of the unscented transfrom parameters</param>
        /// <param name="dim">The state vector dimention</param>
        /// <returns></returns>
        public static (UTParams, UTParams) SampleVectorToUTParams(Vector<double> P, int dim)
        {
            int n = P.Count();
            if (n % 2 != 0)
                new ArgumentException("Sample vector size must be even", "Sample vector");

            int L = dim;

            double[] s1 = P.Take(n / 2).ToArray();
            UTParams p1 = new UTParams(L, s1);

            double[] s2 = P.Skip(n / 2).ToArray();
            UTParams p2 = new UTParams(L, s2);

            return (p1, p2);
        }

        /// <summary>
        /// Wrapper for CalculateCriterionValue function. 
        /// Calculates the criterion value given the provided unscented transform parameters for the forecast and the correction phases. 
        /// </summary>
        /// <param name="Phi1">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t</param>
        /// <param name="Phi2">Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi1">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Psi2">Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Mw">Mean of the noise in the dynamics equation </param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Mnu">Mean of the noise in the obseration equation </param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="P">Sample vector to be transformed to a couple of the unscented transfrom parameters</param>
        /// <param name="T">The upper bound of the observation interval</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xhat0">Initial condition</param>
        /// <param name="DX0Hat">Initial condition covariance</param>
        /// <returns>The criterion value for the unscented transfrom parameters obtained from the sample vactor</returns>
        public static double CalculateSampleCriterion(Func<int, Vector<double>, Vector<double>> Phi1,
                                                    Func<int, Vector<double>, Matrix<double>> Phi2,
                                                    Func<int, Vector<double>, Vector<double>> Psi1,
                                                    Func<int, Vector<double>, Matrix<double>> Psi2,
                                                    Vector<double> Mw,
                                                    Matrix<double> Rw,
                                                    Vector<double> Mnu,
                                                    Matrix<double> Rnu,
                                                    Func<Matrix<double>, double> Crit,
                                                    Vector<double> P,
                                                    int T,
                                                    DiscreteVectorModel[] models,
                                                    Vector<double> xhat0,
                                                    Matrix<double> DX0Hat
                                     )
        {
            (UTParams p1, UTParams p2) = SampleVectorToUTParams(P, xhat0.Count);
            return CalculateCriterionValue(Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, Crit, p1, p2, T, models, xhat0, DX0Hat);
        }

        /// <summary>
        /// Wrapper for CalculateStepwiseCriterionValue function. 
        /// Calculates the criterion value given the provided unscented transform parameters for the forecast and the correction phases. 
        /// </summary>
        /// <param name="Phi1">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t</param>
        /// <param name="Phi2">Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi1">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Psi2">Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Mw">Mean of the noise in the dynamics equation </param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Mnu">Mean of the noise in the obseration equation </param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="P">Sample vector to be transformed to a couple of the unscented transfrom parameters</param>
        /// <param name="t">Step</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xHat">array of the state estimates on the previous step for all the given models</param>
        /// <param name="PHat">array of the covariance matricies on the previous step for all the given models</param>
        /// <returns>The criterion value for the particular unscented transform parameters</returns>
        public static double CalculateSampleStepwiseCriterion(Func<int, Vector<double>, Vector<double>> Phi1,
                                                    Func<int, Vector<double>, Matrix<double>> Phi2,
                                                    Func<int, Vector<double>, Vector<double>> Psi1,
                                                    Func<int, Vector<double>, Matrix<double>> Psi2,
                                                    Vector<double> Mw,
                                                    Matrix<double> Rw,
                                                    Vector<double> Mnu,
                                                    Matrix<double> Rnu,
                                                    Func<Matrix<double>, double> Crit,
                                                    Vector<double> P,
                                                    int t,
                                                    DiscreteVectorModel[] models,
                                                    Vector<double>[] xHat,
                                                    Matrix<double>[] PHat
                                                     )
        {
            (UTParams p1, UTParams p2) = SampleVectorToUTParams(P, xHat[0].Count);
            var result = CalculateStepwiseCriterionValue(Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, Crit, p1, p2, t, models, xHat, PHat);
            return result;
        }


 
        /// <summary>
        /// Calculates the criterion value for the estimate given the particular unscented transform parameters
        /// </summary>
        /// <param name="Phi1">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t</param>
        /// <param name="Phi2">Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi1">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Psi2">Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Mw">Mean of the noise in the dynamics equation </param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Mnu">Mean of the noise in the obseration equation </param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="p1">Unscented transfrom parameters for the forecast phase</param>
        /// <param name="p2">Unscented transfrom parameters for the correction phase</param>
        /// <param name="T">The upper bound of the observation interval</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xhat0">Initial condition</param>
        /// <param name="DX0Hat">Initial condition covariance</param>
        /// <returns>The criterion value for the particular unscented transform parameters</returns>
        public static double CalculateCriterionValue(Func<int, Vector<double>, Vector<double>> Phi1,
                                                    Func<int, Vector<double>, Matrix<double>> Phi2,
                                                    Func<int, Vector<double>, Vector<double>> Psi1,
                                                    Func<int, Vector<double>, Matrix<double>> Psi2,
                                                    Vector<double> Mw,
                                                    Matrix<double> Rw,
                                                    Vector<double> Mnu,
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



                for (int t = 1; t < T; t++)
                {
                    for (int i = 0; i < N; i++)
                    {
                        (xHatU[i], PHatU[i]) = Step(Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, p1, p2, t, models[i].Trajectory[t][1], xHatU[i], PHatU[i]);
                    }

                    Vector<double>[] states = models.Select(x => (x.Trajectory[t][0])).ToArray();
                    Matrix<double> errorUPow2 = Exts.Cov(states.Subtract(xHatU), states.Subtract(xHatU));

                    crit = Crit(errorUPow2);
                }
            }
            catch (Exception e)
            {
                crit = double.MaxValue;
            }
            return crit;
        }

        /// <summary>
        /// Calculates the criterion value for the estimate on a particular step given the particular unscented transform parameters
        /// </summary>
        /// <param name="Phi1">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t</param>
        /// <param name="Phi2">Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi1">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Psi2">Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Mw">Mean of the noise in the dynamics equation </param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Mnu">Mean of the noise in the obseration equation </param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="Crit">Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T))  </param>
        /// <param name="p1">Unscented transfrom parameters for the forecast phase</param>
        /// <param name="p2">Unscented transfrom parameters for the correction phase</param>
        /// <param name="t">Step</param>
        /// <param name="models">Discrete vector model samples</param>
        /// <param name="xHat">array of the state estimates on the previous step for all the given models</param>
        /// <param name="PHat">array of the covariance matricies on the previous step for all the given models</param>
        /// <returns>The criterion value for the particular unscented transform parameters</returns>
        public static double CalculateStepwiseCriterionValue(Func<int, Vector<double>, Vector<double>> Phi1,
                                                    Func<int, Vector<double>, Matrix<double>> Phi2,
                                                    Func<int, Vector<double>, Vector<double>> Psi1,
                                                    Func<int, Vector<double>, Matrix<double>> Psi2,
                                                    Vector<double> Mw,
                                                    Matrix<double> Rw,
                                                    Vector<double> Mnu,
                                                    Matrix<double> Rnu,
                                                    Func<Matrix<double>, double> Crit,
                                                    UTParams p1,
                                                    UTParams p2,
                                                    int t,
                                                    DiscreteVectorModel[] models,
                                                    Vector<double>[] xHat,
                                                    Matrix<double>[] PHat
                                                     )
        {
            double crit = 0;
            try
            {
                int N = models.Count();
                Vector<double>[] xHatU = new Vector<double>[N];
                for (int i = 0; i < N; i++)
                {
                    (xHatU[i], _) = Step(Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, p1, p2, t, models[i].Trajectory[t][1], xHat[i], PHat[i]);
                }

                Vector<double>[] states = models.Select(x => (x.Trajectory[t][0])).ToArray();
                Matrix<double> errorUPow2 = Exts.Cov(states.Subtract(xHatU), states.Subtract(xHatU));

                crit = Crit(errorUPow2);
            }
            catch (Exception e)
            {
                crit = double.MaxValue;
            }
            return crit;
        }

        /// <summary>
        /// Performs a step of Unscented Kalman Filter given the particular unscented transform parameters
        /// for forecast and correction phases
        /// </summary>
        /// <param name="Phi1">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t</param>
        /// <param name="Phi2">Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi1">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Psi2">Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Mw">Mean of the noise in the dynamics equation </param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Mnu">Mean of the noise in the obseration equation </param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="p1">Unscented transfrom parameters for the forecast phase</param>
        /// <param name="p2">Unscented transfrom parameters for the correction phase</param>
        /// <param name="t">Current step time instant</param>
        /// <param name="y">Observations on the current step</param>
        /// <param name="xHat_">Estimate on the previous step</param>
        /// <param name="P_">Approximated previous step error covariance</param>        
        /// <param name="xHat">Returns: current step estimate</param>
        /// <param name="P">Returns: approximated current step error covariance</param>        
        public static (Vector<double>, Matrix<double>) Step(Func<int, Vector<double>, Vector<double>> Phi1,
                                Func<int, Vector<double>, Matrix<double>> Phi2,
                                Func<int, Vector<double>, Vector<double>> Psi1,
                                Func<int, Vector<double>, Matrix<double>> Psi2,
                                Vector<double> Mw,
                                Matrix<double> Rw,
                                Vector<double> Mnu,
                                Matrix<double> Rnu,
                                UTParams p1,
                                UTParams p2,
                                int t,
                                Vector<double> y,
                                Vector<double> xHat_,
                                Matrix<double> P_)
        {
            try
            {
                UnscentedTransform.Transform(x => Phi1(t, x) + Phi2(t, x) * Mw, xHat_, P_, Phi2(t, xHat_) * Rw * Phi2(t, xHat_).Transpose(), p1, out Vector<double> Xtilde, out _, out Matrix<double> Ptilde);
                //UnscentedTransform.Transform(x => Phi1(t, x), xHat_, P_, Rw, p1, out Vector<double> Xtilde2, out _, out Matrix<double> Ptilde2);
                UnscentedTransform.Transform(x => Psi1(t, x) + Psi2(t, x) * Mnu, Xtilde, Ptilde, Psi2(t, Xtilde) * Rnu * Psi2(t, Xtilde).Transpose(), p2, out Vector<double> Ytilde, out Matrix<double> PXY, out Matrix<double> PYtilde);
                //UnscentedTransform.Transform(x => Psi1(t, x), Xtilde, Ptilde, Rnu, p2, out Vector<double> Ytilde2, out Matrix<double> PXY2, out Matrix<double> PYtilde2);

                Matrix<double> K = PXY * PYtilde.Inverse();
                return (Xtilde + K * (y - Ytilde), Ptilde - K * PYtilde * K.Transpose());
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                return (xHat_, P_);
            }
            //Matrix<double> K2 = PXY2 * PYtilde2.Inverse();
            //Vector<double> xHat2 = Xtilde2 + K2 * (y - Ytilde2);
            //Matrix<double> PHat2 = Ptilde2 - K2 * PYtilde2 * K2.Transpose();

            //Console.WriteLine(xHat - xHat2);
            //Console.WriteLine(PHat - PHat2);

        }


        /// <summary>
        /// Performs a step of Unscented Kalman Filter with fixed transformation parameters
        /// for forecast and correction phases (utParamsForecast and utParamsCorrection must be initialized)
        /// </summary>
        /// <param name="Phi1">State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t</param>
        /// <param name="Phi2">Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t</param>
        /// <param name="Psi1">Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Psi2">Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t</param>
        /// <param name="Mw">Mean of the noise in the dynamics equation </param>
        /// <param name="Rw">Covariance matrix of the state disturbances</param>
        /// <param name="Mnu">Mean of the noise in the obseration equation </param>
        /// <param name="Rnu">Convariance matrix of the observation noise</param>
        /// <param name="t">Current step time instant</param>
        /// <param name="y">Observations on the current step</param>
        /// <param name="xHat_">Estimate on the previous step</param>
        /// <param name="P_">Approximated previous step error covariance</param>        
        /// <param name="xHat">Returns: current step estimate</param>
        /// <param name="P">Returns: approximated current step error covariance</param>        
        public (Vector<double>, Matrix<double>) Step(Func<int, Vector<double>, Vector<double>> Phi1,
                                Func<int, Vector<double>, Matrix<double>> Phi2,
                                Func<int, Vector<double>, Vector<double>> Psi1,
                                Func<int, Vector<double>, Matrix<double>> Psi2,
                                Vector<double> Mw,
                                Matrix<double> Rw,
                                Vector<double> Mnu,
                                Matrix<double> Rnu,
                                int t,
                                Vector<double> y,
                                Vector<double> xHat_,
                                Matrix<double> P_)
        {
            if (utParamsForecastStepwise != null && utParamsCorrectionStepwise != null)
                return Step(Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, utParamsForecastStepwise[t], utParamsCorrectionStepwise[t], t, y, xHat_, P_);
            else
                return Step(Phi1, Phi2, Psi1, Psi2, Mw, Rw, Mnu, Rnu, utParamsForecast, utParamsCorrection, t, y, xHat_, P_);
        }
    }
}

