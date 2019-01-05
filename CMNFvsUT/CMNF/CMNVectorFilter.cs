using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using MathNetExtensions;
using System.Xml.Serialization;
using System.Xml;
using System.Xml.Schema;

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
            Func<Vector<double>, bool> Sifter = (x) => Math.Sqrt(x[0] * x[0] + x[1] * x[1]) > 1000.0;

            int n_total = models.Count();

            Vector<double>[] xHat = Enumerable.Repeat(xhat0, n_total).ToArray();
            //bool[] inUse = Enumerable.Repeat(true, n_total).ToArray();
            Console.WriteLine($"CMNF estimate parameters start");
            DateTime start = DateTime.Now;
            for (int t = 1; t < T; t++) // start from 1 because for 0 we do not have observations
            {
                //int n = inUse.Where(e => e).Count();
                DateTime startiteration = DateTime.Now;
                Vector<double>[] x = new Vector<double>[n_total];
                Vector<double>[] y = new Vector<double>[n_total];
                Vector<double>[] xiHat = new Vector<double>[n_total];
                //int k = 0;
                for (int i = 0; i < n_total; i++)
                {
                    //if (inUse[i])
                    //{
                        x[i] = models[i].Trajectory[t][0];
                        y[i] = models[i].Trajectory[t][1];
                        xiHat[i] = Xi(t, xHat[i]);
                       // k++;
                    //}
                }

                Matrix<double> CovXiHat = Exts.Cov(xiHat, xiHat);
                Matrix<double> InvCovXiHat = Matrix<double>.Build.Dense(CovXiHat.RowCount, CovXiHat.ColumnCount, 0.0);
                if (CovXiHat.FrobeniusNorm() > 0)
                    try
                    {
                        //InvCovXiHat = CovXiHat.PseudoInverse();
                        InvCovXiHat = CovXiHat.Inverse(1e-32, 1e-32);
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine("Can't inverse XiHat");
                        Console.WriteLine(CovXiHat.ToString());
                        Console.WriteLine(e.Message);
                    }
                Matrix<double> F = Exts.Cov(x, xiHat) * InvCovXiHat;
                Vector<double> f = x.Average() - F * xiHat.Average();
                Matrix<double> kTilde = Exts.Cov(x, x) - Exts.Cov(x, xiHat) * F.Transpose();

                Vector<double>[] xTilde = new Vector<double>[n_total];
                Vector<double>[] zetaTilde = new Vector<double>[n_total];
                for (int i = 0; i < n_total; i++)
                {
                    xTilde[i] = F * xiHat[i] + f;
                    zetaTilde[i] = Zeta(t, xTilde[i], y[i], kTilde);
                }


                Matrix<double> CovZetaTilde = Exts.Cov(zetaTilde, zetaTilde);
                Matrix<double> InvCovZetaTilde = Matrix<double>.Build.Dense(CovZetaTilde.RowCount, CovZetaTilde.ColumnCount, 0.0);
                try
                {
                    //InvCovZetaTilde = CovZetaTilde.PseudoInverse();
                    InvCovZetaTilde = CovZetaTilde.Inverse(1e-32, 1e-32);
                }
                catch (Exception e)
                {
                    Console.WriteLine("Can't inverse ZetaTilde");
                    Console.WriteLine(CovZetaTilde.ToString());
                    Console.WriteLine(e.Message);
                }
                Matrix<double> H = Exts.Cov(x.Subtract(xTilde), zetaTilde) * InvCovZetaTilde;
                Vector<double> h = -H * zetaTilde.Average();

                Matrix<double> kHat = kTilde - Exts.Cov(x.Subtract(xTilde), zetaTilde) * H.Transpose();

                //k = 0;
                for (int i = 0; i < n_total; i++)
                {

                    //if (inUse[i])
                    //{
                        xHat[i] = xTilde[i] + H * zetaTilde[i] + h;
                    //    inUse[i] = !Sifter(x[k] - xHat[k]);
                    //    k++;
                    //}
                }
                FHat.Add(t, F);
                fHat.Add(t, f);
                HHat.Add(t, H);
                hHat.Add(t, h);


                KTilde.Add(t, kTilde);
                KHat.Add(t, kHat);
                //KHat.Add(t, Exts.Cov(x, x) - Exts.Cov(x, xiHat) * F - Exts.Cov(x.Subtract(xTilde), zetaTilde) * H);
                // KHat.Add(t, cov(x, x) - cov(x, xiHat) * F - cov(x - xTilde, zetaTilde) * H);

                Matrix<double> I = Matrix<double>.Build.DenseIdentity(CovXiHat.RowCount);

                //Console.WriteLine();
                //Console.WriteLine($"cov_xi_{t} = {CovXiHat.ToMatlab()}");
                //Console.WriteLine();
                //Console.WriteLine($"inv_by_cov_xi_{t} = {(InvCovXiHat * CovXiHat).ToMatlab()}");
                //Console.WriteLine();
                //Console.WriteLine($"cov_zeta_{t} = {CovZetaTilde.ToMatlab()}");
                //Console.WriteLine();
                //Console.WriteLine($"inv_by_cov_zeta_{t} = {(InvCovZetaTilde * CovZetaTilde).ToMatlab()}");
                //Console.WriteLine();

                //Console.WriteLine();
                //Console.WriteLine($"$cov(\\hat{{\\xi}}_{{{t}}}, \\hat{{\\xi}}_{{{t}}}) = {CovXiHat.ToLatex("E3")}$");
                //Console.WriteLine();
                //Console.WriteLine($"$cov(\\hat{{\\xi}}_{{{t}}}, \\hat{{\\xi}}_{{{t}}}) \\cdot (cov(\\hat{{\\xi}}_{{{t}}}, \\hat{{\\xi}}_{{{t}}}))^{{-1}}  = {(InvCovXiHat * CovXiHat).ToLatex("F3")}$");
                //Console.WriteLine();
                //double diff_xi = (I - InvCovXiHat * CovXiHat).L2Norm();
                //Console.WriteLine($"$|| I - cov(\\hat{{\\xi}}_{{{t}}}, \\hat{{\\xi}}_{{{t}}}) \\cdot (cov(\\hat{{\\xi}}_{{{t}}}, \\hat{{\\xi}}_{{{t}}}))^{{-1}}||  = {diff_xi}$");
                //Console.WriteLine();
                //if (diff_xi > 0.001)
                //    Console.WriteLine("!!!");

                //Matrix<double> II = Matrix<double>.Build.DenseIdentity(CovZetaTilde.RowCount);

                //Console.WriteLine($"$cov(\\tilde{{\\zeta}}_{{{t}}}, \\tilde{{\\zeta}}_{{{t}}}) = {CovZetaTilde.ToLatex("E3")}$");
                //Console.WriteLine();
                //Console.WriteLine($"$cov(\\tilde{{\\zeta}}_{{{t}}}, \\tilde{{\\zeta}}_{{{t}}}) \\cdot (cov(\\tilde{{\\zeta}}_{{{t}}}, \\tilde{{\\zeta}}_{{{t}}}))^{{-1}}  = {(InvCovZetaTilde * CovZetaTilde).ToLatex("F3")}$");
                //Console.WriteLine();
                //double diff_zeta = (II - InvCovZetaTilde * CovZetaTilde).L2Norm();
                //Console.WriteLine($"$|| I - cov(\\tilde{{\\zeta}}_{{{t}}}, \\tilde{{\\zeta}}_{{{t}}}) \\cdot (cov(\\tilde{{\\zeta}}_{{{t}}}, \\tilde{{\\zeta}}_{{{t}}}))^{{-1}}||  = {diff_zeta}$");
                //Console.WriteLine();
                //if (diff_zeta > 0.001)
                //    Console.WriteLine("!!!");

                //Console.WriteLine($"$F_{{{t}}} = {F.ToLatex()}$");
                //Console.WriteLine($"$f_{{{t}}} = {f.ToLatex()}$");
                //Console.WriteLine($"$H_{{{t}}} = {H.ToLatex()}$");
                //Console.WriteLine($"$h_{{{t}}} = {h.ToLatex()}$");
                //Console.WriteLine($"$\\tilde{{K}}_{{{t}}} = {kTilde.ToLatex()}$");
                //Console.WriteLine($"$\\hat{{K}}_{{{t}}} = {kHat.ToLatex()}$");
                Console.WriteLine($"CMNF estimate parameters for t={t}, done in {(DateTime.Now - startiteration).ToString(@"hh\:mm\:ss\.fff")}");
                x = null;
                y = null;
                xiHat = null;
            }
            DateTime finish = DateTime.Now;
            Console.WriteLine($"CMNF estimate parameters finished in {(finish - start).ToString(@"hh\:mm\:ss\.fff")}");


        }

        public (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat_)
        {
            Vector<double> xTilde = FHat[t] * Xi(t, xHat_) + fHat[t];
            Vector<double> xHat = xTilde + HHat[t] * Zeta(t, xTilde, y, KTilde[t]) + hHat[t];
            return (xHat, KHat[t]);
        }
    }


}
