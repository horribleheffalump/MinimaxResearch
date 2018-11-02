using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using MathNetExtensions;
using System.Threading.Tasks;

namespace CMNF
{
    public class CMNFilter
    {
        Func<int, Vector<double>, Vector<double>> Xi;
        Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> Zeta;

        public Dictionary<int, Matrix<double>> FHat;
        public Dictionary<int, Vector<double>> fHat;
        public Dictionary<int, Matrix<double>> HHat;
        public Dictionary<int, Vector<double>> hHat;
        public Dictionary<int, Matrix<double>> KTilde;
        public Dictionary<int, Matrix<double>> KHat;

        public CMNFilter(Func<int, Vector<double>, Vector<double>> xi, Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> zeta)
        {
            Xi = xi;
            Zeta = zeta;

            FHat = new Dictionary<int, Matrix<double>>();
            fHat = new Dictionary<int, Vector<double>>();
            HHat = new Dictionary<int, Matrix<double>>();
            hHat = new Dictionary<int, Vector<double>>();
            KTilde = new Dictionary<int, Matrix<double>>();
            KHat = new Dictionary<int, Matrix<double>>();
        }

        public void EstimateParameters(DiscreteVectorModel[] models, Vector<double> xhat0, int T)
        {
            int n = models.Count();

            Vector<double>[] xHat = Enumerable.Repeat(xhat0, n).ToArray();
            Console.WriteLine($"CMNF estimate parameters start");
            DateTime start = DateTime.Now;
            //for (int t = 0; t < T; t++)
            Parallel.For(0, T, t =>
            {
                DateTime startiteration = DateTime.Now;
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
                Matrix<double> kTilde = Exts.Cov(x, x) - Exts.Cov(x, xiHat) * F.Transpose();

                Vector<double>[] xTilde = new Vector<double>[n];
                Vector<double>[] zetaTilde = new Vector<double>[n];
                for (int i = 0; i < n; i++)
                {
                    xTilde[i] = F * xiHat[i] + f;
                    zetaTilde[i] = Zeta(t, xTilde[i], y[i], kTilde);
                }

                Matrix<double> H = Exts.Cov(x.Subtract(xTilde), zetaTilde) * (Exts.Cov(zetaTilde, zetaTilde).PseudoInverse());
                Vector<double> h = -H * zetaTilde.Average();

                Matrix<double> kHat = kTilde - Exts.Cov(x.Subtract(xTilde), zetaTilde) * H.Transpose();
                for (int i = 0; i < n; i++)
                {
                    //xHat[i] = F* xiHat[i] +f + H * zetaTilde[i] + h;
                    xHat[i] = xTilde[i] + H * zetaTilde[i] + h;
                }
                FHat.Add(t, F);
                fHat.Add(t, f);
                HHat.Add(t, H);
                hHat.Add(t, h);


                KTilde.Add(t, kTilde);
                KHat.Add(t, kHat);
                //KHat.Add(t, Exts.Cov(x, x) - Exts.Cov(x, xiHat) * F - Exts.Cov(x.Subtract(xTilde), zetaTilde) * H);
                // KHat.Add(t, cov(x, x) - cov(x, xiHat) * F - cov(x - xTilde, zetaTilde) * H);
                Console.WriteLine($"CMNF estimate parameters for t={t}, done in {(DateTime.Now - startiteration).ToString(@"hh\:mm\:ss\.fff")}");
                x = null;
                y = null;
                xiHat = null;
            });
            DateTime finish = DateTime.Now;
            Console.WriteLine($"CMNF estimate parameters finished in {(finish - start).ToString(@"hh\:mm\:ss\.fff")}");


        }

        public (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat_)
        {
            Vector<double> xTilde = FHat[t] * Xi(t, xHat_) + fHat[t];
            Vector<double> xHat = xTilde + HHat[t] * Zeta(t, xTilde, y, KTilde[t]) + hHat[t];
            return (xHat, KHat[t]);
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
