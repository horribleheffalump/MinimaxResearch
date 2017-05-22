using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;

namespace UKF
{
    public class UTStaticEstimate
    {
        public double Alpha0;
        public double _alpha = 1e-3;
        public double _beta = 2.0;
        public double _kappa = 3.0;

        public void Estimate(Func<Vector<double>, Vector<double>> Phi,Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu, 
            out Vector<double> mErr_UT, out Matrix<double> KErr_UT, out Matrix<double> KErrTh_UT)
        {
            int L = mX.Count;
            //double _alpha = 1e-3;
            //double _beta = 2.0;
            //double _kappa = 3.0;

            //double _Lambda = Math.Pow(_alpha, 2.0) * (L + _kappa) - L;
            Vector<double> _Wm;
            Vector<double> _Wc;
            _Wm = Vector<double>.Build.Dense(2 * L + 1, (1.0 - Alpha0) / 4.0);
            _Wm[0] = Alpha0;
            _Wc = Vector<double>.Build.Dense(2 * L + 1, (1.0 - Alpha0) / 4.0);
            _Wc[0] = Alpha0;
            double _Lambda = 2.0 / (1 - Alpha0) - L;

            //_Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (_Lambda + L));
            //_Wm[0] = _Lambda / (_Lambda + L);

            //_Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (_Lambda + L));
            //_Wc[0] = _Lambda / (_Lambda + L) + 1.0 - Math.Pow(_alpha, 2.0) + _beta;



            //Matrix<double> KUT = K_UT(x => Phi(x), mX, KX, _Lambda, _Wm, _Wc, out M_UT);
            //Matrix<double> P_UT = KUT.SubMatrix(0, L, L, L) * ((KUT.SubMatrix(L, L, L, L) + KNu).PseudoInverse());
            //Matrix<double> KErrTh_UT = KX - P_UT * (KUT.SubMatrix(L, L, L, L) + KNu) * P_UT.Transpose();

            Vector<double> M_UT;
            Matrix<double> Kxy_UT;
            Matrix<double> Kyy_UT;
            Matrix<double> KUT = K_UT(x => Phi(x), mX, KX, KNu, _Lambda, _Wm, _Wc, out M_UT, out Kxy_UT, out Kyy_UT);
            Matrix<double> P_UT = Kxy_UT * ((Kyy_UT).PseudoInverse());
            KErrTh_UT = KX - P_UT * Kyy_UT * P_UT.Transpose();

            Vector<double>[] Xhat_UT = Y.Select(y => mX + P_UT * (y - M_UT)).ToArray();
            Vector<double>[] Err_UT = Xhat_UT.Subtract(X);
            mErr_UT = Err_UT.Average();
            KErr_UT = Extensions.Cov(Err_UT, Err_UT);

        }

        public void Estimate3(Func<Vector<double>, Vector<double>> Phi, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu,
    out Vector<double> mErr_UT, out Matrix<double> KErr_UT, out Matrix<double> KErrTh_UT)
        {
            int L = mX.Count;

            double _Lambda = Math.Pow(_alpha, 2.0) * (L + _kappa) - L;
            Vector<double> _Wm;
            Vector<double> _Wc;
            _Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (_Lambda + L));
            _Wm[0] = _Lambda / (_Lambda + L);

            _Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (_Lambda + L));
            _Wc[0] = _Lambda / (_Lambda + L) + 1.0 - Math.Pow(_alpha, 2.0) + _beta;



            //Matrix<double> KUT = K_UT(x => Phi(x), mX, KX, _Lambda, _Wm, _Wc, out M_UT);
            //Matrix<double> P_UT = KUT.SubMatrix(0, L, L, L) * ((KUT.SubMatrix(L, L, L, L) + KNu).PseudoInverse());
            //Matrix<double> KErrTh_UT = KX - P_UT * (KUT.SubMatrix(L, L, L, L) + KNu) * P_UT.Transpose();

            Vector<double> M_UT;
            Matrix<double> Kxy_UT;
            Matrix<double> Kyy_UT;
            Matrix<double> KUT = K_UT(x => Phi(x), mX, KX, KNu, _Lambda, _Wm, _Wc, out M_UT, out Kxy_UT, out Kyy_UT);
            Matrix<double> P_UT = Kxy_UT * ((Kyy_UT).PseudoInverse());
            KErrTh_UT = KX - P_UT * Kyy_UT * P_UT.Transpose();

            Vector<double>[] Xhat_UT = Y.Select(y => mX + P_UT * (y - M_UT)).ToArray();
            Vector<double>[] Err_UT = Xhat_UT.Subtract(X);
            mErr_UT = Err_UT.Average();
            KErr_UT = Extensions.Cov(Err_UT, Err_UT);

        }
        public void EstimateParameters(Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {
            Alpha0 = UTParms(Phi, Crit, X, Y, mX, KX, KNu);
        }

        public void EstimateParameters3(Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {
            double[] p = UTParms3(Phi, Crit, X, Y, mX, KX, KNu);
            _alpha = p[0];
            _beta = p[1];
            _kappa = p[2];
        }

        static Matrix<double> K_UT(Func<Vector<double>, Vector<double>> Phi, Vector<double> mX, Matrix<double> dX, Matrix<double> dY, double lambda, Vector<double> wm, Vector<double> wc, out Vector<double> y, out Matrix<double> Kxy, out Matrix<double> Kyy)
        {
            int L = mX.Count;

            Matrix<double> Xi = UnscentedTransform.GenerateSigmaPoints(mX, dX, lambda);

            Matrix<double> Upsilon = Phi(Xi.Column(0)).ToColumnMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {
                Upsilon = Upsilon.Append(Phi(Xi.Column(i)).ToColumnMatrix());
            }

            y = wm[0] * Upsilon.Column(0);
            for (int i = 1; i < 2 * L + 1; i++)
            {
                y = y + wm[i] * Upsilon.Column(i);
            }


            Kxy = wc[0] * (Xi.Column(0) - mX).ToColumnMatrix() * (Upsilon.Column(0) - y).ToRowMatrix();
            Kyy = wc[0] * (Upsilon.Column(0) - y).ToColumnMatrix() * (Upsilon.Column(0) - y).ToRowMatrix();

            Matrix<double> Z = Xi.Stack(Upsilon);
            Vector<double> MZ = mX.Stack(y);
            Matrix<double> PFull = wc[0] * (Z.Column(0) - MZ).ToColumnMatrix() * (Z.Column(0) - MZ).ToRowMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {
                PFull = PFull + wc[i] * (Z.Column(i) - MZ).ToColumnMatrix() * (Z.Column(i) - MZ).ToRowMatrix();
                Kxy = Kxy + wc[i] * (Xi.Column(i) - mX).ToColumnMatrix() * (Upsilon.Column(i) - y).ToRowMatrix();
                Kyy = Kyy + wc[i] * (Upsilon.Column(i) - y).ToColumnMatrix() * (Upsilon.Column(i) - y).ToRowMatrix();
            }
            Kyy = Kyy + dY;
            //Py = Py + R;
            return PFull;
        }

        static double UTParms(Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {
            ContinuousUniform UnifAlpha = new ContinuousUniform(-10, 10);

            AsyncCalculatorPlanner acp = new AsyncCalculatorPlanner(100, 10, () => CalculateCriterionValue(Phi, Crit, UnifAlpha, X, Y, mX, KX, KNu));

            List<double[]> results1 = acp.DoCalculate();

            double min1 = results1.Where(i => !double.IsNaN(i[0])).Min(e => e[0]);
            double[] best1 = results1.First(e => e[0] == min1);

            Normal NormalAlpha = new Normal(best1[1], Math.Sqrt(0.005));

            acp = new AsyncCalculatorPlanner(100, 10, () => CalculateCriterionValue(Phi, Crit, NormalAlpha, X, Y, mX, KX, KNu));
            List<double[]> results2 = acp.DoCalculate();

            double min2 = results2.Where(i => !double.IsNaN(i[0])).Min(e => e[0]);
            double[] best2 = results2.First(e => e[0] == min2);

            double Alpha_0 = best2[1];
            return Alpha_0;
        }

        static double[] UTParms3(Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {
            ContinuousUniform UnifAlpha = new ContinuousUniform(-10, 10);
            ContinuousUniform UnifBeta = new ContinuousUniform(-10, 10);
            ContinuousUniform UnifKappa = new ContinuousUniform(-10, 10);

            AsyncCalculatorPlanner acp = new AsyncCalculatorPlanner(100, 10, () => CalculateCriterionValue3(Phi, Crit, new IContinuousDistribution[] { UnifAlpha, UnifBeta, UnifKappa }, X, Y, mX, KX, KNu));

            List<double[]> results1 = acp.DoCalculate();

            double min1 = results1.Where(i => !double.IsNaN(i[0])).Min(e => e[0]);
            double[] best1 = results1.First(e => e[0] == min1);

            Normal NormalAlpha = new Normal(best1[1], Math.Sqrt(0.005));
            Normal NormalBeta = new Normal(best1[2], Math.Sqrt(0.005));
            Normal NormalKappa = new Normal(best1[3], Math.Sqrt(0.005));

            acp = new AsyncCalculatorPlanner(100, 10, () => CalculateCriterionValue3(Phi, Crit, new IContinuousDistribution[] { NormalAlpha, NormalBeta, NormalKappa }, X, Y, mX, KX, KNu));
            List<double[]> results2 = acp.DoCalculate();

            double min2 = results2.Where(i => !double.IsNaN(i[0])).Min(e => e[0]);
            double[] best2 = results2.First(e => e[0] == min2);
            return new double[] { best2[1] , best2[2] , best2[3] };
        }
        public static double[] CalculateCriterionValue3(Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, IContinuousDistribution[] dirtribution, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {
            double _alpha = dirtribution[0].Sample();
            double _beta = dirtribution[1].Sample(); ;
            double _kappa = dirtribution[2].Sample(); ;

            int L = mX.Count;
            double _Lambda = Math.Pow(_alpha, 2.0) * (L + _kappa) - L;
            Vector<double> _Wm;
            Vector<double> _Wc;
            _Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (_Lambda + L));
            _Wm[0] = _Lambda / (_Lambda + L);

            _Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (_Lambda + L));
            _Wc[0] = _Lambda / (_Lambda + L) + 1.0 - Math.Pow(_alpha, 2.0) + _beta;


            double crit = double.MaxValue;
            try
            {
                //Matrix<double> KUT = K_UT(x => Phi(x), mX, KX, _Lambda, _Wm, _Wc, out M_UT);
                //Matrix<double> P_UT = KUT.SubMatrix(0, L, L, L) * ((KUT.SubMatrix(L, L, L, L) + KNu).PseudoInverse());
                //Matrix<double> KErrTh_UT = KX - P_UT * (KUT.SubMatrix(L, L, L, L) + KNu) * P_UT.Transpose();

                Vector<double> M_UT;
                Matrix<double> Kxy_UT;
                Matrix<double> Kyy_UT;
                Matrix<double> KUT = K_UT(x => Phi(x), mX, KX, KNu, _Lambda, _Wm, _Wc, out M_UT, out Kxy_UT, out Kyy_UT);
                Matrix<double> P_UT = Kxy_UT * ((Kyy_UT).PseudoInverse());
                Matrix<double> KErrTh_UT = KX - P_UT * Kyy_UT * P_UT.Transpose();

                Vector<double>[] Xhat_UT = Y.Select(y => mX + P_UT * (y - M_UT)).ToArray();
                Vector<double>[] Err_UT = Xhat_UT.Subtract(X);
                Vector<double> mErr_UT = Err_UT.Average();
                Matrix<double> KErr_UT = Extensions.Cov(Err_UT, Err_UT);
                crit = Crit(KErr_UT);//.Trace();
            }
            catch { }
            return new double[] { crit, _alpha, _beta, _kappa };
        }

        public static double[] CalculateCriterionValue(Func<Vector<double>, Vector<double>> Phi, Func<Matrix<double>, double> Crit, IContinuousDistribution dirtribution, Vector<double>[] X, Vector<double>[] Y, Vector<double> mX, Matrix<double> KX, Matrix<double> KNu)
        {
            double Alpha_0 = dirtribution.Sample();
            int L = mX.Count;

            Vector<double> _Wm;
            Vector<double> _Wc;
            _Wm = Vector<double>.Build.Dense(2 * L + 1, (1.0 - Alpha_0) / 4.0);
            _Wm[0] = Alpha_0;
            _Wc = Vector<double>.Build.Dense(2 * L + 1, (1.0 - Alpha_0) / 4.0);
            _Wc[0] = Alpha_0;
            double _Lambda = 2.0 / (1 - Alpha_0) - L;

            double crit = double.MaxValue;
            try
            {
                //Matrix<double> KUT = K_UT(x => Phi(x), mX, KX, _Lambda, _Wm, _Wc, out M_UT);
                //Matrix<double> P_UT = KUT.SubMatrix(0, L, L, L) * ((KUT.SubMatrix(L, L, L, L) + KNu).PseudoInverse());
                //Matrix<double> KErrTh_UT = KX - P_UT * (KUT.SubMatrix(L, L, L, L) + KNu) * P_UT.Transpose();

                Vector<double> M_UT;
                Matrix<double> Kxy_UT;
                Matrix<double> Kyy_UT;
                Matrix<double> KUT = K_UT(x => Phi(x), mX, KX, KNu, _Lambda, _Wm, _Wc, out M_UT, out Kxy_UT, out Kyy_UT);
                Matrix<double> P_UT = Kxy_UT * ((Kyy_UT).PseudoInverse());
                Matrix<double> KErrTh_UT = KX - P_UT * Kyy_UT * P_UT.Transpose();

                Vector<double>[] Xhat_UT = Y.Select(y => mX + P_UT * (y - M_UT)).ToArray();
                Vector<double>[] Err_UT = Xhat_UT.Subtract(X);
                Vector<double> mErr_UT = Err_UT.Average();
                Matrix<double> KErr_UT = Extensions.Cov(Err_UT, Err_UT);
                crit = Crit(KErr_UT);//.Trace();
            }
            catch { }
            return new double[] { crit, Alpha_0 };
        }

    }
}
