using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;


namespace CMNF
{
    public class CMNScalarFilter
    {
        Func<double, double> Xi;
        Func<double, double, double> Zeta;

        public Dictionary<int, double> FHat;
        public Dictionary<int, double> fHat;
        public Dictionary<int, double> HHat;
        public Dictionary<int, double> hHat;
        public Dictionary<int, double> KHat;

        public CMNScalarFilter(Func<double, double> xi, Func<double, double, double> zeta)
        {
            Xi = xi;
            Zeta = zeta;

            FHat = new Dictionary<int, double>();
            fHat = new Dictionary<int, double>();
            HHat = new Dictionary<int, double>();
            hHat = new Dictionary<int, double>();
            KHat = new Dictionary<int, double>();
        }

        public void EstimateParameters(DiscreteScalarModel[] models, double xhat0, int T)
        {
            int n = models.Count();
            Vector<double> xHat = Vector<double>.Build.Dense(n, xhat0);
            for (int t = 0; t < T; t++)
            {
                for (int i = 0; i < n; i++)
                {
                    models[i].Step();
                }
                Vector<double> x = Vector<double>.Build.Dense(n, (i) => models[i].State);
                Vector<double> y = Vector<double>.Build.Dense(n, (i) => models[i].Obs);

                Vector<double> xiHat = Vector<double>.Build.Dense(n, (i) => Xi(xHat[i]));

                double F = cov(x, xiHat) / cov(xiHat, xiHat);
                double f = x.Average() - F * xiHat.Average();

                Vector<double> xTilde = F * xiHat + f;

                Vector<double> zetaTilde = Vector<double>.Build.Dense(n, (i) => Zeta(xTilde[i], y[i]));

                double H = cov(x - xTilde, zetaTilde) / cov(zetaTilde, zetaTilde);
                double h = -H * zetaTilde.Average();

                xHat = Vector<double>.Build.Dense(n, (i) => F*xiHat[i] + f + H*zetaTilde[i] + h);

                FHat.Add(t, F);
                fHat.Add(t, f);
                HHat.Add(t, H);
                hHat.Add(t, h);

                KHat.Add(t, cov(x, x) - cov(x, xiHat) * F - cov(x - xTilde, zetaTilde) * H);
            }

        }

        public void EstimateParameters(DiscreteVectorModel[] models, Vector<double> xhat0, int T)
        {
            int n = models.Count();
            Vector<double> xHat = Vector<double>.Build.Dense(n, xhat0[0]);
            //Console.WriteLine($"xHat={xHat}");
            for (int t = 0; t < T; t++)
            {
                //for (int i = 0; i < n; i++)
                //{
                //    models[i].Step();
                //}
                //Vector<double> x = Vector<double>.Build.Dense(n, (i) => models[i].State[0]);
                //Vector<double> y = Vector<double>.Build.Dense(n, (i) => models[i].Obs[0]);

                Vector<double> x = Vector<double>.Build.Dense(n, (i) => models[i].Trajectory[t][0][0]);
                Vector<double> y = Vector<double>.Build.Dense(n, (i) => models[i].Trajectory[t][1][0]);

                Vector<double> xiHat = Vector<double>.Build.Dense(n, (i) => Xi(xHat[i]));



                double F = cov(x, xiHat) / cov(xiHat, xiHat);
                double f = x.Average() - F * xiHat.Average();

                Vector<double> xTilde = F * xiHat + f;

                Vector<double> zetaTilde = Vector<double>.Build.Dense(n, (i) => Zeta(xTilde[i], y[i]));

                double H = cov(x - xTilde, zetaTilde) / cov(zetaTilde, zetaTilde);
                double h = -H * zetaTilde.Average();

                xHat = Vector<double>.Build.Dense(n, (i) => F * xiHat[i] + f + H * zetaTilde[i] + h);

                FHat.Add(t, F);
                fHat.Add(t, f);
                HHat.Add(t, H);
                hHat.Add(t, h);

                KHat.Add(t, cov(x, x) - cov(x, xiHat) * F - cov(x - xTilde, zetaTilde) * H);

                //Console.WriteLine($"x={x}");
                //Console.WriteLine($"y={y}");
                //Console.WriteLine($"xiHat={xiHat}");
                //Console.WriteLine($"F={F}");
                //Console.WriteLine($"f={f}");
                //Console.WriteLine($"xTilde={xTilde}");
                //Console.WriteLine($"zetaTilde={zetaTilde}");
                //Console.WriteLine($"H={H}");
                //Console.WriteLine($"h={h}");
                //Console.WriteLine($"xHat={xHat}");
                //Console.WriteLine($"KHat={KHat}");



            }

        }


        public double Step(int t, double y, double xHat_)
        {
            double xTilde = FHat[t] * Xi(xHat_) + fHat[t];
            double xHat = xTilde + HHat[t] * Zeta(xTilde, y) + hHat[t];
            return xHat;
        }

        private double cov(Vector<double> x, Vector<double> y)
        {

            double r1 = ((x - x.Average()).PointwiseMultiply(y - y.Average())).Sum() / (x.Count - 1.0);
            double r2 = (1.0 / (x.Count - 1.0)) * x.DotProduct(y) - x.Average() * y.Average();

            double n = x.Count;
            double r3 = (x.DotProduct(y) * n - x.Sum() * y.Sum()) / n / n;

            return r1;
        }
    }
}
