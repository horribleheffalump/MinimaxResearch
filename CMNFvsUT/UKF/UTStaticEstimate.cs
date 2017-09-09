using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;

namespace UKF
{
    public class UTStaticEstimate
    {
        public UTParams optimalParams;
        UTOptimizationType OptimizationType;

        public UTStaticEstimate(UTOptimizationType type = UTOptimizationType.ImplicitAlpha)
        {
            OptimizationType = type;
        }

        public Vector<double>[] Estimate(Func<Vector<double>, Vector<double>> Phi, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu,
            out Vector<double> mErr_UT, out Matrix<double> KErr_UT, out Matrix<double> KErrTh_UT)
        {
            Vector<double> M_UT;
            Matrix<double> Kxy_UT;
            Matrix<double> Kyy_UT;
            Matrix<double> KUT = K_UT(x => Phi(x), mX, KX, KNu, optimalParams, out M_UT, out Kxy_UT, out Kyy_UT);
            Matrix<double> P_UT = Kxy_UT * ((Kyy_UT).PseudoInverse());
            KErrTh_UT = KX - P_UT * Kyy_UT * P_UT.Transpose();

            Vector<double>[] Xhat_UT = Y.Select(y => mX + P_UT * (y - M_UT)).ToArray();
            Vector<double>[] Err_UT = Xhat_UT.Subtract(X);
            mErr_UT = Err_UT.Average();
            KErr_UT = Utils.Cov(Err_UT, Err_UT);

            return Xhat_UT;
        }

        public void EstimateParameters(int N1, int N2, Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {
            optimalParams = UTParmsOptimize(N1, N2, OptimizationType, Phi, Crit, X, Y, mX, KX, KNu);
            //optimalParams = UTParmsOptimize(N1, N2, UTOptimizationType.ImplicitAlpha, Phi, Crit, X, Y, mX, KX, KNu);
            //optimalParams = UTParmsOptimize(N1, N2, UTOptimizationType.ImplicitAlphaBetaKappa, Phi, Crit, X, Y, mX, KX, KNu);
            //optimalParams = UTParmsOptimize(N1, N2, UTOptimizationType.Explicit, Phi, Crit, X, Y, mX, KX, KNu);

        }

        static UTParams UTParmsOptimize(int N1, int N2, UTOptimizationType type, Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {

            int n = 0;
            string filename = "..\\..\\..\\output\\UT_optimization_{type}.txt";
            switch (type)
            {
                case UTOptimizationType.ImplicitAlpha: n = 1; filename = filename.Replace("{type}", "ImplicitAlpha"); break;
                case UTOptimizationType.ImplicitAlphaBetaKappa: n = 3; filename = filename.Replace("{type}", "ImplicitABK"); break;
                case UTOptimizationType.Explicit: n = 4; filename = filename.Replace("{type}", "Explicit"); break;
            }

            IContinuousDistribution[] distr = new IContinuousDistribution[n];
            for (int i = 0; i< n; i++)
            {
                distr[i] = new ContinuousUniform(-15, 15);
            } 

            AsyncCalculatorPlanner acp = new AsyncCalculatorPlanner(N1, 10, () => CalculateSampleCriterion(Phi, Crit, distr, X, Y, mX, KX, KNu));

            List<double[]> results1 = acp.DoCalculate();

            double min1 = results1.Where(i => !double.IsNaN(i[0])).Min(e => e[0]);
            double[] best1 = results1.First(e => e[0] == min1);

            for (int i = 0; i < n; i++)
            {
                distr[i] = new Normal(best1[i+5], Math.Sqrt(0.005));
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

        public static double[] CalculateSampleCriterion(Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, IContinuousDistribution[] distribution, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {
            double[] s = distribution.Select(d => d.Sample()).ToArray();
            int L = mX.Count;

            UTParams p = new UTParams(mX.Count, s);
            double crit = CalculateCriterionValue(Phi, Crit, p, X, Y, mX, KX, KNu);
            return (new double[] { crit, p.Lambda, p.Wm[0], p.Wc[0], p.Wm[1] }).Concat(s).ToArray();
        }


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

    }

    public enum UTOptimizationType { ImplicitAlpha, ImplicitAlphaBetaKappa, Explicit }

    public class UTParams
    {
        public double Lambda;
        public Vector<double> Wm;
        public Vector<double> Wc;

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

        public void SetUTParams(int L, double alpha0)
        {
            Lambda = 2.0 / (1 - alpha0) - L;
            Wm = Vector<double>.Build.Dense(2 * L + 1, (1.0 - alpha0) / 4.0);
            Wm[0] = alpha0;
            Wc = Vector<double>.Build.Dense(2 * L + 1, (1.0 - alpha0) / 4.0);
            Wc[0] = alpha0;
        }

        public void SetUTParams(int L, double alpha, double beta, double kappa)
        {
            Lambda = Math.Pow(alpha, 2.0) * (L + kappa) - L;
            Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wm[0] = Lambda / (Lambda + L);

            Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wc[0] = Lambda / (Lambda + L) + 1.0 - Math.Pow(alpha, 2.0) + beta;
        }

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


