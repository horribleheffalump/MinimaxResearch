using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using System.IO;

namespace UKF
{
    /// <summary>
    /// Unscented transform estimate for the model y = Phi(x) + Nu, where 
    /// - x is a random variable with known mean and covariance,
    /// - Nu - noise with zero mean and known covariance
    /// Usage:
    /// - specify the parameters of the unscented transform in utProperty manually or by means of primitive optmization procedure  <see cref="UTStaticEstimate.EstimateParameters()"/>
    /// - calculate the estimate <see cref="UTStaticEstimate.Estimate()"/>
    /// It should be noted, that the train and test sets may be different. 
    /// That is, the arrays of samples X = [x_0,...,x_N] and observations Y = [y_0,...,y_N] = [Phi(x_0) + Nu_0,...,Phi(x_N) + Nu_N]  may vary 
    /// for the step of the unscented transform parameters optimization and the step of unscented transform estimate calculation.
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
            Vector<double> M_UT;
            Matrix<double> Kxy_UT;
            Matrix<double> Kyy_UT;
            Matrix<double> KUT = K_UT(x => Phi(x), mX, KX, KNu, utParams, out M_UT, out Kxy_UT, out Kyy_UT);
            Matrix<double> P_UT = Kxy_UT * ((Kyy_UT).PseudoInverse());
            KErrTh_UT = KX - P_UT * Kyy_UT * P_UT.Transpose();

            Vector<double>[] Xhat_UT = Y.Select(y => mX + P_UT * (y - M_UT)).ToArray();
            Vector<double>[] Err_UT = Xhat_UT.Subtract(X);
            mErr_UT = Err_UT.Average();
            KErr_UT = Utils.Cov(Err_UT, Err_UT);

            return Xhat_UT;
        }

        /// <summary>
        /// The unscented transform
        /// </summary>
        /// <param name="Phi">Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu</param>
        /// <param name="mX">Mean of the transformed random variable</param>
        /// <param name="dX">Cov of the ransformed random variable</param>
        /// <param name="dY">Cov of the additive random variable</param>
        /// <param name="p">Parameters of the unscented transform</param>
        /// <param name="y">Returns: approximated mean of the transformed variable</param>
        /// <param name="Kxy">Returns: approximated cross-covariance of the initial and the transformed variable</param>
        /// <param name="Kyy">Returns: approximated covariance of the transormed variable</param>
        static Matrix<double> K_UT(Func<Vector<double>, Vector<double>> Phi, Vector<double> mX, Matrix<double> dX, Matrix<double> dY, UTParams p, out Vector<double> y, out Matrix<double> Kxy, out Matrix<double> Kyy)
        {
            int L = mX.Count;

            Matrix<double> Xi = UnscentedTransform.GenerateSigmaPoints(mX, dX, p.Lambda);

            Matrix<double> Upsilon = Phi(Xi.Column(0)).ToColumnMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {  
                Upsilon = Upsilon.Append(Phi(Xi.Column(i)).ToColumnMatrix());
            }

            y = p.Wm[0] * Upsilon.Column(0);
            for (int i = 1; i < 2 * L + 1; i++)
            {
                y = y + p.Wm[i] * Upsilon.Column(i);
            }


            Kxy = p.Wc[0] * (Xi.Column(0) - mX).ToColumnMatrix() * (Upsilon.Column(0) - y).ToRowMatrix();
            Kyy = p.Wc[0] * (Upsilon.Column(0) - y).ToColumnMatrix() * (Upsilon.Column(0) - y).ToRowMatrix();

            Matrix<double> Z = Xi.Stack(Upsilon);
            Vector<double> MZ = mX.Stack(y);
            Matrix<double> PFull = p.Wc[0] * (Z.Column(0) - MZ).ToColumnMatrix() * (Z.Column(0) - MZ).ToRowMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {
                PFull = PFull + p.Wc[i] * (Z.Column(i) - MZ).ToColumnMatrix() * (Z.Column(i) - MZ).ToRowMatrix();
                Kxy = Kxy + p.Wc[i] * (Xi.Column(i) - mX).ToColumnMatrix() * (Upsilon.Column(i) - y).ToRowMatrix();
                Kyy = Kyy + p.Wc[i] * (Upsilon.Column(i) - y).ToColumnMatrix() * (Upsilon.Column(i) - y).ToRowMatrix();
            }
            Kyy = Kyy + dY;
            //Py = Py + R;
            return PFull;
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
        /// Unscented transform parameters "optimization" procedure. 
        /// - Step 1: generate N1 random samples, calculate the unscented transform estimates given the parameters determined by each random sample.
        /// - Step 2:choose the random sample with the best estimate quality criterion value and generate N2 random samples in closer area of this sample.
        /// - Step 3:again choose the random sample with the best estimate quality criterion value, save the samples ordered by the criterion value to the output file and return the best found unscented transform parameters.
        /// The UTOptimizationType type param determines the way how the random samples define the unscented tranform params <see cref="UTParams"/>.
        /// - If type is UTOptimizationType.ImplicitAlpha, then random samples define alpha0 - scalar weight of the central points for both sample mean and cov <see cref="UTParams.SetUTParams(int, double)"/>;
        /// - If type is UTOptimizationType.ImplicitAlphaBetaKappa, then random samples are vectors of dim 3 and represent three parameters alpha, beta, kappa which are then transformed to the parameters of the inscented transform <see cref="UTParams.SetUTParams(int, double, double, double)"/>; 
        /// - If type is UTOptimizationType.Explicit, then random samples are vectors of dim 4 and explicitly define the unscented transform parameters <see cref="UTParams.SetUTParams(int, double, double, double, double)"/>. ///TODO it is not right to define the parameters of the unsctnted transform arbitraty, they have to be interdependent, so that the mean and cov would be transformed correctly.
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

            int n = 0;
            string filename = Path.Combine(outputFolder, "UT_optimization_{type}.txt");
            switch (type)
            {
                case UTDefinitionType.ImplicitAlpha: n = 1; filename = filename.Replace("{type}", "ImplicitAlpha"); break;
                case UTDefinitionType.ImplicitAlphaBetaKappa: n = 3; filename = filename.Replace("{type}", "ImplicitABK"); break;
                case UTDefinitionType.Explicit: n = 4; filename = filename.Replace("{type}", "Explicit"); break;
            }

            IContinuousDistribution[] distr = new IContinuousDistribution[n];
            for (int i = 0; i < n; i++)
            {
                distr[i] = new ContinuousUniform(-15, 15);
            }

            AsyncCalculatorPlanner acp = new AsyncCalculatorPlanner(N1, 10, () => CalculateSampleCriterion(Phi, Crit, distr, X, Y, mX, KX, KNu));

            List<double[]> results1 = acp.DoCalculate();

            double min1 = results1.Where(i => !double.IsNaN(i[0])).Min(e => e[0]);
            double[] best1 = results1.First(e => e[0] == min1);

            for (int i = 0; i < n; i++)
            {
                distr[i] = new Normal(best1[i + 5], Math.Sqrt(0.005));
            }

            acp = new AsyncCalculatorPlanner(N2, 10, () => CalculateSampleCriterion(Phi, Crit, distr, X, Y, mX, KX, KNu));
            List<double[]> results2 = acp.DoCalculate();

            double min2 = results2.Where(i => !double.IsNaN(i[0])).Min(e => e[0]);
            double[] best2 = results2.First(e => e[0] == min2);


            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(filename))
            {
                NumberFormatInfo provider;
                provider = new NumberFormatInfo();
                provider.NumberDecimalSeparator = ".";

                var results = results1.Concat(results2).Where(i => !double.IsNaN(i[0]) && i[0] < double.MaxValue);
                foreach (var e in results.OrderBy(e => e[0]))
                {
                    outputfile.WriteLine(String.Join(",", e.Select(s => string.Format(provider, "{0}", s))));
                }

            }
            return new UTParams(mX.Count, best2.Skip(1 + n).Take(n).ToArray());
        }

        /// <summary>
        /// Generates a random sample for the unscented transform parameters and calculates the criterion value for the unscented transform estimate.
        /// The size of the distribution parameter determines the unscented transform parameters definition method <see cref="UTParams"/>
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
            return (new double[] { crit, p.Lambda, p.Wm[0], p.Wc[0], p.Wm[1] }).Concat(s).ToArray();
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
            double crit = double.MaxValue;
            try
            {
                Vector<double> M_UT;
                Matrix<double> Kxy_UT;
                Matrix<double> Kyy_UT;
                Matrix<double> KUT = K_UT(x => Phi(x), mX, KX, KNu, p, out M_UT, out Kxy_UT, out Kyy_UT);
                Matrix<double> P_UT = Kxy_UT * ((Kyy_UT).PseudoInverse());
                Matrix<double> KErrTh_UT = KX - P_UT * Kyy_UT * P_UT.Transpose();

                Vector<double>[] Xhat_UT = Y.Select(y => mX + P_UT * (y - M_UT)).ToArray();
                Vector<double>[] Err_UT = Xhat_UT.Subtract(X);
                Vector<double> mErr_UT = Err_UT.Average();
                Matrix<double> KErr_UT = Utils.Cov(Err_UT, Err_UT);
                crit = Crit(KErr_UT);//.Trace();
            }
            catch { }
            return crit;
        }

    }

    /// <summary>
    /// Unscented transform parameters optimization type. 
    /// </summary>
    public enum UTDefinitionType { ImplicitAlpha, ImplicitAlphaBetaKappa, Explicit }

    /// <summary>
    /// Parameters of the unscented transform: 
    /// - Lambda - scaling parameter
    /// - Wm - weights for sample mean
    /// - Wc - weights for sample covariance
    /// Three ways to define: 
    /// - explicitly,
    /// - implicitly with one parameter alpha0, 
    /// - implicitly with three parameters alpha, beta, kappa
    /// </summary>
    public class UTParams
    {
        public double Lambda;
        public Vector<double> Wm;
        public Vector<double> Wc;
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="L">Dimension of the transformed random variable</param>
        /// <param name="p">Input params to define the unscented transform params. The appropriate method to define UT params is chosen depending on the size of the array.</param>
        public UTParams(int L, params double[] p)
        {
            int n = p.Count();
            switch (n)
            {
                case 1: SetUTParams(L, p[0]); break;
                case 3: SetUTParams(L, p[0], p[1], p[2]); break;
                case 4: SetUTParams(L, p[0], p[1], p[2], p[3]); break;
            }
        }

        /// <summary>
        /// The parameters of the unscented transform are defined by a single parameter alpha0:
        /// </summary>
        /// <param name="L">Dimension of the transformed random variable</param>
        /// <param name="alpha0">Alpha0 - weight of the central points for both sample mean and cov </param>
        public void SetUTParams(int L, double alpha0)
        {
            Lambda = 2.0 / (1 - alpha0) - L;
            Wm = Vector<double>.Build.Dense(2 * L + 1, (1.0 - alpha0) / 4.0);
            Wm[0] = alpha0;
            Wc = Vector<double>.Build.Dense(2 * L + 1, (1.0 - alpha0) / 4.0);
            Wc[0] = alpha0;
        }

        /// <summary>
        /// The parameters of the unscented transform are defined by three parameters: alpha, beta, kappa
        /// </summary>
        /// <param name="L">Dimension of the transformed random variable</param>
        /// <param name="alpha">Alpha - determines the spread of the sigma points around the mean of the transformed random variable (small positive value 0 \leq alpha \leq 10e-4)</param>
        /// <param name="beta">Beta - is used to incorporate prior knowledge of the distribution of the transformed random variable (for Gaussian b = 2 is optimal)</param>
        /// <param name="kappa">Kappa - is a secondary scaling parameter</param>
        public void SetUTParams(int L, double alpha, double beta, double kappa)
        {
            Lambda = Math.Pow(alpha, 2.0) * (L + kappa) - L;
            Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wm[0] = Lambda / (Lambda + L);

            Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wc[0] = Lambda / (Lambda + L) + 1.0 - Math.Pow(alpha, 2.0) + beta;
        }

        /// <summary>
        /// Explicit definition of the unscented transform parameters
        /// </summary>
        /// <param name="L">Dimension of the transformed random variable</param>
        /// <param name="lambda">Scaling parameter</param>
        /// <param name="wm0">Central point weight for the sample mean</param>
        /// <param name="wc0">Central point weight for the sample cov</param>
        /// <param name="wi">Non-central points weight for sample mean and cov</param>
        public void SetUTParams(int L, double lambda, double wm0, double wc0, double wi)
        {
            Lambda = lambda;
            Wm = Vector<double>.Build.Dense(2 * L + 1, wi);
            Wm[0] = wm0;

            Wc = Vector<double>.Build.Dense(2 * L + 1, wi);
            Wc[0] = wc0;
        }
    }
}

    
