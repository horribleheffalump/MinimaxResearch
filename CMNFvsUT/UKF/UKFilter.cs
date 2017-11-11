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
    public class UKFilter
    {
        private int L; //state dim
        public double Alpha; //spread param
        public double Beta;  //distribution param
        public double Kappa; //secondary scale param 3 - L;
        private double Lambda; //scale param

        private Vector<double> Wm;
        private Vector<double> Wc;

        public Func<int, Vector<double>, Vector<double>> Phi;
        public Func<int, Vector<double>, Vector<double>> Psi;
        public Matrix<double> Rw;
        public Matrix<double> Rnu;

        ContinuousUniform UnifAlpha;
        ContinuousUniform UnifBeta;
        ContinuousUniform UnifKappa;

        Normal NormalAlpha;
        Normal NormalBeta;
        Normal NormalKappa;

        private int t = 0;

        public UKFilter(int l, Func<int, Vector<double>, Vector<double>> phi, Func<int, Vector<double>, Vector<double>> psi, Matrix<double> Rw, Matrix<double> Rnu, double alpha = 1e-3, double beta = 2.0, double kappa = double.NaN)
        {
            L = l;
            Phi = phi;
            Psi = psi;
            this.Rw = Rw;
            this.Rnu = Rnu;

            Alpha = alpha;
            Beta = beta;
            if (double.IsNaN(kappa))
            {
                kappa = 3.0 - L;
                //kappa = 0;
            }
            Kappa = kappa;
            Lambda = Math.Pow(Alpha, 2.0) * (L + Kappa) - L;

            Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wm[0] = Lambda / (Lambda + L);

            Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wc[0] = Lambda / (Lambda + L) + 1.0 - Math.Pow(Alpha, 2.0) + Beta;
        }

        public void RecalcParams()
        {
            Lambda = Math.Pow(Alpha, 2.0) * (L + Kappa) - L;

            Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wm[0] = Lambda / (Lambda + L);

            Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wc[0] = Lambda / (Lambda + L) + 1.0 - Math.Pow(Alpha, 2.0) + Beta;
        }

        public void EstimateParametersRandom(DiscreteVectorModel[] models, int T, Vector<double> xhat0, Matrix<double> DX0Hat, string fileName)
        {
            UnifAlpha = new ContinuousUniform(-10, 10);
            UnifBeta = new ContinuousUniform(-10, 10);
            UnifKappa = new ContinuousUniform(-10, 10);

            AsyncCalculatorPlanner acp = new AsyncCalculatorPlanner(100, 10, () => CalculateCriterionValueAtRandomUniform(models, T, xhat0, DX0Hat));
            List<double[]> results1 = acp.DoCalculate();

            double min1 = results1.Min(e => e[0]);
            double[] best1 = results1.First(e => e[0] == min1);

            NormalAlpha = new Normal(best1[1], Math.Sqrt(0.5));
            NormalBeta = new Normal(best1[2], Math.Sqrt(0.5));
            NormalKappa = new Normal(best1[3], Math.Sqrt(0.5));

            acp = new AsyncCalculatorPlanner(100, 10, () => CalculateCriterionValueAtRandomNormal(models, T, xhat0, DX0Hat));
            List<double[]> results2 = acp.DoCalculate();

            double min2 = results2.Min(e => e[0]);
            double[] best2 = results2.First(e => e[0] == min2);


            Alpha = best2[1];
            Beta = best2[2];
            Kappa = best2[3];

            RecalcParams();

            NumberFormatInfo provider;
            provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
            {
                foreach (var e in results2.OrderBy(e => e[0]))
                {
                    outputfile.WriteLine(string.Format(provider, "{0} {1} {2} {3}", e[0], e[1], e[2], e[3]));
                }
                //outputfile.WriteLine("=======================");
                foreach (var e in results1.OrderBy(e => e[0]))
                {
                    outputfile.WriteLine(string.Format(provider, "{0} {1} {2} {3}", e[0], e[1], e[2], e[3]));
                }



                outputfile.Close();
            }
        }

        public void EstimateParameters(DiscreteVectorModel[] models, int T, Vector<double> xhat0, Matrix<double> DX0Hat, string fileName)
        {
            //double[] alpha = new double[] { 1e-4, 0.5e-3, 1e-3, 0.5e-2, 1e-2, 0.5-1, 1e-1, 0.5, 1 };
            double[] alpha = new double[] { 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5 };
            //double[] beta = new double[] { 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4 };
            //double[] kappa = new double[] { 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4 };
            //double[] alpha = new double[] { 0.5, 0.6, 0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
            double[] beta = new double[] { 2};
            double[] kappa = new double[] { 2};

            Dictionary<double, double[]> results = new Dictionary<double, double[]>();

            for (int i = 0; i < alpha.Length; i++)
                for (int j = 0; j < beta.Length; j++)
                    for (int k = 0; k < kappa.Length; k++)
                    {
                        try
                        {
                            results.Add(CalculateCriterionValue(models, T, xhat0, DX0Hat, alpha[i], beta[j], kappa[k]), new double[] { alpha[i], beta[j], kappa[k] });
                        }
                        catch { }//for the same keys
                    }

            double min = results.Min(e => e.Key);
            double[] best = results[min];
            Alpha = best[0];
            Beta = best[1];
            Kappa = best[2];

            RecalcParams();

            NumberFormatInfo provider;
            provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
            {
                foreach (var e in results.OrderBy(x => x.Key))
                {
                    outputfile.WriteLine(string.Format(provider, "{0} {1} {2} {3}", e.Key, e.Value[0], e.Value[1], e.Value[2]));
                }
                outputfile.Close();
            }



        }

        public double[] CalculateCriterionValueAtRandomUniform(DiscreteVectorModel[] models, int T, Vector<double> xhat0, Matrix<double> DX0Hat)
        {
            double _alpha = UnifAlpha.Sample();
            double _beta = UnifBeta.Sample();
            double _kappa = UnifKappa.Sample();
            return new double[] { CalculateCriterionValue(models, T, xhat0, DX0Hat, _alpha, _beta, _kappa), _alpha, _beta, _kappa };
        }

        public double[] CalculateCriterionValueAtRandomNormal(DiscreteVectorModel[] models, int T, Vector<double> xhat0, Matrix<double> DX0Hat)
        {
            double _alpha = NormalAlpha.Sample();
            double _beta = NormalBeta.Sample();
            double _kappa = NormalKappa.Sample();
            return new double[] { CalculateCriterionValue(models, T, xhat0, DX0Hat, _alpha, _beta, _kappa), _alpha, _beta, _kappa };
        }


        public double CalculateCriterionValue(DiscreteVectorModel[] models, int T, Vector<double> xhat0, Matrix<double> DX0Hat, double _alpha, double _beta, double _kappa)
        {
            double result = 0;
            try
            {
                double _Lambda = Math.Pow(_alpha, 2.0) * (L + _kappa) - L;
                Vector<double> _Wm;
                Vector<double> _Wc;
                _Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (_Lambda + L));
                _Wm[0] = _Lambda / (_Lambda + L);

                _Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (_Lambda + L));
                _Wc[0] = _Lambda / (_Lambda + L) + 1.0 - Math.Pow(_alpha, 2.0) + _beta;

                int N = models.Count();

                Vector<double>[] xHatU = models.Select(x => xhat0).ToArray();
                Matrix<double>[] PHatU = models.Select(x => DX0Hat).ToArray();
                //Vector<double> PHatU = Vector<double>.Build.Dense(N, DX0Hat);



                for (int t = 0; t < T; t++)
                {
                    for (int i = 0; i < N; i++)
                    {
                        Vector<double> xHatU_i;
                        Matrix<double> PHatU_i;
                        Step(models[i].Trajectory[t][1],
                            xHatU[i],
                            PHatU[i],
                            _Lambda, _Wm, _Wc,
                            out xHatU_i, out PHatU_i);
                        xHatU[i] = xHatU_i;
                        PHatU[i] = PHatU_i;
                    }

                    Vector<double>[] states = models.Select(x => (x.Trajectory[t][0])).ToArray();
                    Matrix<double> errorUPow2 = Utils.Cov(states.Subtract(xHatU), states.Subtract(xHatU));

                    //result = errorUPow2.Trace();
                    result = errorUPow2[1,1];
                    //result += errorUPow2.Average();
                }
            }
            catch { result = double.MaxValue; }
            return result;
        }

        public double CalculateCriterionValue(DiscreteScalarModel[] models, int T, double xhat0, double DX0Hat, double _alpha, double _beta, double _kappa)
        {
            double result = 0;
            try
            {
                double _Lambda = Math.Pow(_alpha, 2.0) * (L + _kappa) - L;
                Vector<double> _Wm;
                Vector<double> _Wc;
                _Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (_Lambda + L));
                _Wm[0] = _Lambda / (_Lambda + L);

                _Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (_Lambda + L));
                _Wc[0] = _Lambda / (_Lambda + L) + 1.0 - Math.Pow(_alpha, 2.0) + _beta;


                int N = models.Count();

                Vector<double> xHatU = Vector<double>.Build.Dense(N, xhat0);
                Vector<double> PHatU = Vector<double>.Build.Dense(N, DX0Hat);



                for (int t = 0; t < T; t++)
                {
                    //Vector<double> y = Vector<double>.Build.Dense(N, i => models[i].Step());
                    //Vector<double> x = Vector<double>.Build.Dense(N, i => models[i].State);
                    Vector<double> y = Vector<double>.Build.Dense(N, i => models[i].Trajectory[t][1]);
                    Vector<double> x = Vector<double>.Build.Dense(N, i => models[i].Trajectory[t][0]);

                    for (int i = 0; i < N; i++)
                    {
                        Vector<double> xHatU_i;
                        Matrix<double> PHatU_i;
                        Step(Vector<double>.Build.Dense(1, y[i]),
                            Vector<double>.Build.Dense(1, xHatU[i]),
                            Matrix<double>.Build.Dense(1, 1, PHatU[i]),
                            _Lambda, _Wm, _Wc,
                            out xHatU_i, out PHatU_i);
                        xHatU[i] = xHatU_i[0];
                        PHatU[i] = PHatU_i[0, 0];
                    }

                    Vector<double> errorUPow2 = (x - xHatU).PointwiseMultiply(x - xHatU);
                    result = errorUPow2.Average();
                    //result += errorUPow2.Average();
                }
            }
            catch { result = double.MaxValue; }
            return result; // / T;
        }

        public Matrix<double> GenerateSigmaPoints(Vector<double> x, Matrix<double> P)
        {
            return UnscentedTransform.GenerateSigmaPoints(x, P, Lambda);
        }



        public void UT(Func<Vector<double>, Vector<double>> f, Matrix<double> Xi, Matrix<double> R, Vector<double> wm, Vector<double> wc, out Matrix<double> Upsilon, out Vector<double> y, out Matrix<double> Py)
        {
            Upsilon = f(Xi.Column(0)).ToColumnMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {
                Upsilon = Upsilon.Append(f(Xi.Column(i)).ToColumnMatrix());
            }

            y = wm[0] * Upsilon.Column(0);
            for (int i = 1; i < 2 * L + 1; i++)
            {
                y = y + wm[i] * Upsilon.Column(i);
            }

            Py = wc[0] * (Upsilon.Column(0) - y).ToColumnMatrix() * (Upsilon.Column(0) - y).ToRowMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {
                Py = Py + wc[i] * (Upsilon.Column(i) - y).ToColumnMatrix() * (Upsilon.Column(i) - y).ToRowMatrix();
            }
            Py = Py + R;

        }

        public void UT(Func<Vector<double>, Vector<double>> f, Matrix<double> Xi, Matrix<double> R, out Matrix<double> Upsilon, out Vector<double> y, out Matrix<double> Py)
        {
            UT(f, Xi, R, Wm, Wc, out Upsilon, out y, out Py);
        }

        public void Step(Vector<double> y, Vector<double> xHat_, Matrix<double> P_, double lambda, Vector<double> wm, Vector<double> wc, out Vector<double> xHat, out Matrix<double> PHat)
        {
            Matrix<double> Xi_ = UnscentedTransform.GenerateSigmaPoints(xHat_, P_, lambda);

            Matrix<double> XiStar;
            Vector<double> Xtilde;
            Matrix<double> Ptilde;

            UT(x => Phi(t,x), Xi_, Rw, wm, wc, out XiStar, out Xtilde, out Ptilde);

            //Matrix<double> Xiplus = GenerateSigmaPoints(XiStar.Column(0), Rw);
            Matrix<double> Xi = UnscentedTransform.GenerateSigmaPoints(Xtilde, Ptilde, lambda);

            Matrix<double> Upsilon;
            Vector<double> Ytilde;
            Matrix<double> PYtilde;

            UT(x => Psi(t,x), Xi, Rnu, wm, wc, out Upsilon, out Ytilde, out PYtilde);

            Matrix<double> PXY = wc[0] * (Xi.Column(0) - Xtilde).ToColumnMatrix() * (Upsilon.Column(0) - Ytilde).ToRowMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {
                PXY = PXY + wc[i] * (Xi.Column(i) - Xtilde).ToColumnMatrix() * (Upsilon.Column(i) - Ytilde).ToRowMatrix();
            }

            Matrix<double> K = PXY * PYtilde.Inverse();

            xHat = Xtilde + K * (y - Ytilde);
            PHat = Ptilde - K * PYtilde * K.Transpose();

            t++;
        }

        public void Step(Vector<double> y, Vector<double> xHat_, Matrix<double> P_, out Vector<double> xHat, out Matrix<double> PHat)
        {
            Step(y, xHat_, P_, Lambda, Wm, Wc, out xHat, out PHat);
            //Matrix<double> Xi_ = GenerateSigmaPoints(xHat_, P_);

            //Matrix<double> XiStar;
            //Vector<double> Xtilde;
            //Matrix<double> Ptilde;

            //UT(Phi, Xi_, Rw, out XiStar, out Xtilde, out Ptilde);

            ////Matrix<double> Xiplus = GenerateSigmaPoints(XiStar.Column(0), Rw);
            //Matrix<double> Xi = GenerateSigmaPoints(Xtilde, Ptilde);

            //Matrix<double> Upsilon;
            //Vector<double> Ytilde;
            //Matrix<double> PYtilde;

            //UT(Psi, Xi, Rnu, out Upsilon, out Ytilde, out PYtilde);

            //Matrix<double> PXY = Wc[0] * (Xi.Column(0) - Xtilde).ToColumnMatrix() * (Upsilon.Column(0) - Ytilde).ToRowMatrix();
            //for (int i = 1; i < 2 * L + 1; i++)
            //{
            //    PXY = PXY + Wc[i] * (Xi.Column(i) - Xtilde).ToColumnMatrix() * (Upsilon.Column(i) - Ytilde).ToRowMatrix();
            //}

            //Matrix<double> K = PXY * PYtilde.Inverse();

            //xHat = Xtilde + K * (y - Ytilde);
            //PHat = Ptilde - K * PYtilde * K.Transpose();
        }


    }

    public static class UnscentedTransform
    {
        public static Matrix<double> GenerateSigmaPoints(Vector<double> x, Matrix<double> P, double lambda)
        {
            int L = x.Count;

            Matrix<double> Sqrt = (Math.Sqrt(lambda + L)) * P.Cholesky().Factor.Transpose();

            Matrix<double> Xi = x.ToColumnMatrix();
            for (int i = 0; i < L; i++)
            {
                Xi = Xi.Append((x + Sqrt.Column(i)).ToColumnMatrix());
            }
            for (int i = 0; i < L; i++)
            {
                Xi = Xi.Append((x - Sqrt.Column(i)).ToColumnMatrix());
            }

            return Xi;
        }
    }
}
