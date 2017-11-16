using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using MathNetExtensions;

namespace CMNF
{
    public class CMNFilter
    {
        Func<int, Vector<double>, Vector<double>> Xi;
        Func<int, Vector<double>, Vector<double>, Vector<double>> Zeta;

        public Dictionary<int, Matrix<double>> FHat;
        public Dictionary<int, Vector<double>> fHat;
        public Dictionary<int, Matrix<double>> HHat;
        public Dictionary<int, Vector<double>> hHat;
        public Dictionary<int, Matrix<double>> KHat;

        public CMNFilter(Func<int, Vector<double>, Vector<double>> xi, Func<int, Vector<double>, Vector<double>, Vector<double>> zeta)
        {
            Xi = xi;
            Zeta = zeta;

            FHat = new Dictionary<int, Matrix<double>>();
            fHat = new Dictionary<int, Vector<double>>();
            HHat = new Dictionary<int, Matrix<double>>();
            hHat = new Dictionary<int, Vector<double>>();
            KHat = new Dictionary<int, Matrix<double>>();
        }

        public void EstimateParameters(DiscreteVectorModel[] models, Vector<double> xhat0, int T)
        {
            int n = models.Count();

            Vector<double>[] xHat = Enumerable.Repeat(xhat0, n).ToArray();

            for (int t = 0; t < T; t++)
            {
                Vector<double>[] x = new Vector<double>[n];
                Vector<double>[] y = new Vector<double>[n];
                Vector<double>[] xiHat = new Vector<double>[n];
                for (int i = 0; i < n; i++)
                {
                    //models[i].Step();
                    //x[i] = models[i].State;
                    //y[i] = models[i].Obs;
                    x[i] = models[i].Trajectory[t][0];
                    y[i] = models[i].Trajectory[t][1];
                    xiHat[i] = Xi(t, xHat[i]);
                }

                Matrix<double> F = Exts.Cov(x, xiHat) * (Exts.Cov(xiHat, xiHat).PseudoInverse());
                Vector<double> f = x.Average() - F * xiHat.Average();

                Vector<double>[] xTilde = new Vector<double>[n];
                Vector<double>[] zetaTilde = new Vector<double>[n];
                for (int i = 0; i < n; i++)
                {
                    xTilde[i] = F * xiHat[i] + f;
                    zetaTilde[i] = Zeta(t, xTilde[i], y[i]);
                }

                Matrix<double> H = Exts.Cov(x.Subtract(xTilde), zetaTilde) * (Exts.Cov(zetaTilde, zetaTilde).PseudoInverse());
                Vector<double> h = -H * zetaTilde.Average();

                for (int i = 0; i < n; i++)
                {
                    xHat[i] = F * xiHat[i] + f + H * zetaTilde[i] + h;
                }
                FHat.Add(t, F);
                fHat.Add(t, f);
                HHat.Add(t, H);
                hHat.Add(t, h);

                KHat.Add(t, Exts.Cov(x, x) - Exts.Cov(x, xiHat) * F - Exts.Cov(x.Subtract(xTilde), zetaTilde) * H);
               // KHat.Add(t, cov(x, x) - cov(x, xiHat) * F - cov(x - xTilde, zetaTilde) * H);

            }

        }

        public Vector<double> Step(int t, Vector<double> y, Vector<double> xHat_)
        {
            Vector<double> xTilde = FHat[t] * Xi(t, xHat_) + fHat[t];
            Vector<double> xHat = xTilde + HHat[t] * Zeta(t, xTilde, y) + hHat[t];
            return xHat;
        }

        //private double cov(Vector<double> x, Vector<double> y)
        //{

        //    double r1 = ((x - x.Average()).PointwiseMultiply(y - y.Average())).Average();
        //    double r2 = (1.0 / x.Count) * x.DotProduct(y) - x.Average() * y.Average();

        //    double n = x.Count;
        //    double r3 = (x.DotProduct(y) * n - x.Sum() * y.Sum()) / n / n;

        //    return r2;
        //}

    }


}
