using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using MathNetExtensions;
using MathNet.Numerics.Distributions;
using System.Threading.Tasks;

namespace CMNF
{
    public class AdaptiveCMNFilter
    {
        private Func<int, Vector<double>, Vector<double>> Phi1; // Phi1(t, X)
        private Func<int, Vector<double>, Matrix<double>> Phi2;
        private Func<int, Vector<double>, Vector<double>> Psi1;
        private Func<int, Vector<double>, Matrix<double>> Psi2;
        private Func<int, Vector<double>> W;
        private Func<int, Vector<double>> Nu;

        private DiscreteVectorModel[] Models;



        Func<int, Vector<double>, Vector<double>> Xi;
        Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> Zeta;

        public Dictionary<int, Matrix<double>> FHat;
        public Dictionary<int, Vector<double>> fHat;
        public Dictionary<int, Matrix<double>> HHat;
        public Dictionary<int, Vector<double>> hHat;
        public Dictionary<int, Matrix<double>> KTilde;
        public Dictionary<int, Matrix<double>> KHat;
        public Vector<double>[][] xiHat;

        public AdaptiveCMNFilter(Func<int, Vector<double>, Vector<double>> xi,
                                    Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> zeta,
                                    Func<int, Vector<double>, Vector<double>> phi1,
                                    Func<int, Vector<double>, Matrix<double>> phi2,
                                    Func<int, Vector<double>, Vector<double>> psi1,
                                    Func<int, Vector<double>, Matrix<double>> psi2,
                                    Func<int, Vector<double>> w,
                                    Func<int, Vector<double>> nu)
        {
            Xi = xi;
            Zeta = zeta;
            Phi1 = phi1;
            Phi2 = phi2;
            Psi1 = psi1;
            Psi2 = psi2;
            W = w;
            Nu = nu;

            FHat = new Dictionary<int, Matrix<double>>();
            fHat = new Dictionary<int, Vector<double>>();
            HHat = new Dictionary<int, Matrix<double>>();
            hHat = new Dictionary<int, Vector<double>>();
            KTilde = new Dictionary<int, Matrix<double>>();
            KHat = new Dictionary<int, Matrix<double>>();
        }

        public void EstimateParameters(DiscreteVectorModel[] models, Vector<double> xhat0, int T)
        {
            Models = models;
            int n = models.Count();
            xiHat = new Vector<double>[T][];

            Vector<double>[] xHat = Enumerable.Repeat(xhat0, n).ToArray();
            Console.WriteLine($"CMNF estimate parameters start");
            DateTime start = DateTime.Now;
            for (int t = 1; t < T; t++) // start from 1 because for 0 we do not have observations
            {
                DateTime startiteration = DateTime.Now;
                Vector<double>[] x = new Vector<double>[n];
                Vector<double>[] y = new Vector<double>[n];
                xiHat[t] = new Vector<double>[n];
                for (int i = 0; i < n; i++)
                {
                    x[i] = models[i].Trajectory[t][0];
                    y[i] = models[i].Trajectory[t][1];
                    xiHat[t][i] = Xi(t, xHat[i]);
                }

                Matrix<double> CovXiHat = Exts.Cov(xiHat[t], xiHat[t]);
                Matrix<double> InvCovXiHat = Matrix<double>.Build.Dense(CovXiHat.RowCount, CovXiHat.ColumnCount, 0.0);
                if (CovXiHat.FrobeniusNorm() > 0)
                    try
                    {
                        InvCovXiHat = CovXiHat.PseudoInverse();
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine("Can't inverse XiHat");
                        Console.WriteLine(CovXiHat.ToString());
                        Console.WriteLine(e.Message);
                    }
                Matrix<double> F = Exts.Cov(x, xiHat[t]) * InvCovXiHat;
                Vector<double> f = x.Average() - F * xiHat[t].Average();
                Matrix<double> kTilde = Exts.Cov(x, x) - Exts.Cov(x, xiHat[t]) * F.Transpose();

                Vector<double>[] xTilde = new Vector<double>[n];
                Vector<double>[] zetaTilde = new Vector<double>[n];
                for (int i = 0; i < n; i++)
                {
                    xTilde[i] = F * xiHat[t][i] + f;
                    zetaTilde[i] = Zeta(t, xTilde[i], y[i], kTilde);
                }


                Matrix<double> CovZetaTilde = Exts.Cov(zetaTilde, zetaTilde);
                Matrix<double> InvCovZetaTilde = Matrix<double>.Build.Dense(CovZetaTilde.RowCount, CovZetaTilde.ColumnCount, 0.0);
                try
                {
                    InvCovZetaTilde = CovZetaTilde.PseudoInverse();
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
                //xiHat = null;
            }
            DateTime finish = DateTime.Now;
            Console.WriteLine($"CMNF estimate parameters finished in {(finish - start).ToString(@"hh\:mm\:ss\.fff")}");


        }



        public (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat_, Matrix<double> kHat_)
        {
            Vector<double>[] x_mod = new Vector<double>[Models.Count()];
            Vector<double>[] y_mod = new Vector<double>[Models.Count()];
            Vector<double>[] xiHat_mod = new Vector<double>[Models.Count()];

            int selected = 0;
            Vector<double> Bound = 5.0 * kHat_.Diagonal().PointwiseSqrt();
            for (int i = 0; i < Models.Count(); i++)
            {
                var x = Models[i].Trajectory[t][0];
                if (x.InBounds(xHat_ - Bound, xHat_ + Bound))
                {
                    x_mod[selected] = x;
                    y_mod[selected] = Models[i].Trajectory[t][1];
                    xiHat_mod[selected] = xiHat[t][i];
                    selected++;
                }
            }

            Console.WriteLine($"selected: {selected}");
            if (selected == Models.Count()) // all test trajectories selected => use CMNF
            {
                Vector<double> xTilde = FHat[t] * Xi(t, xHat_) + fHat[t];
                Vector<double> xHat = xTilde + HHat[t] * Zeta(t, xTilde, y, KTilde[t]) + hHat[t];
                return (xHat, KHat[t]);
            }
            else
            {
                //Parallel.For(0, n, new ParallelOptions() {MaxDegreeOfParallelism = System.Environment.ProcessorCount }, i =>
                RandomVector<Normal> xHatDistr = new RandomVector<Normal>(xHat_, kHat_); 
                for (int i = selected; i < Models.Count(); i++)
                {
                    x_mod[i] = xHatDistr.Sample();

                    if (t > 0)
                    {
                        x_mod[i] = Phi1(t, x_mod[i]) + Phi2(t, x_mod[i]) * W(t);
                    }
                    y_mod[i] = Psi1(t, x_mod[i]) + Psi2(t, x_mod[i]) * Nu(t);

                    xiHat_mod[i] = Xi(t, xHatDistr.Sample());
                }
                //});

                if (selected == 0)  // no test trajectories selected => use MCMNF
                {
                    Vector<double> f = x_mod.Average();
                    Matrix<double> kTilde = Exts.Cov(x_mod, x_mod);

                    Vector<double>[] zetaTilde = new Vector<double>[Models.Count()];
                    for (int i = 0; i < Models.Count(); i++)
                    {
                        zetaTilde[i] = Zeta(t, f, y_mod[i], kTilde);
                    }


                    Matrix<double> CovZetaTilde = Exts.Cov(zetaTilde, zetaTilde);
                    Matrix<double> InvCovZetaTilde = Matrix<double>.Build.Dense(CovZetaTilde.RowCount, CovZetaTilde.ColumnCount, 0.0);
                    try
                    {
                        InvCovZetaTilde = CovZetaTilde.PseudoInverse();
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine("Can't inverse ZetaTilde");
                        Console.WriteLine(CovZetaTilde.ToString());
                        Console.WriteLine(e.Message);
                    }
                    Matrix<double> H = Exts.Cov(x_mod.Subtract(f), zetaTilde) * InvCovZetaTilde;
                    Vector<double> h = -H * zetaTilde.Average();

                    Matrix<double> kHat = kTilde - Exts.Cov(x_mod.Subtract(f), zetaTilde) * H.Transpose();

                    Vector<double> xHat__ = f + H * Zeta(t, f, y, kTilde) + h;

                    return (xHat__, kHat);
                }
                else
                {
                    Matrix<double> CovXiHat = Exts.Cov(xiHat_mod, xiHat_mod);
                    Matrix<double> InvCovXiHat = Matrix<double>.Build.Dense(CovXiHat.RowCount, CovXiHat.ColumnCount, 0.0);
                    if (CovXiHat.FrobeniusNorm() > 0)
                        try
                        {
                            InvCovXiHat = CovXiHat.PseudoInverse();
                        }
                        catch (Exception e)
                        {
                            Console.WriteLine("Can't inverse XiHat");
                            Console.WriteLine(CovXiHat.ToString());
                            Console.WriteLine(e.Message);
                        }
                    Matrix<double> F = Exts.Cov(x_mod, xiHat_mod) * InvCovXiHat;
                    Vector<double> f = x_mod.Average() - F * xiHat_mod.Average();
                    Matrix<double> kTilde = Exts.Cov(x_mod, x_mod) - Exts.Cov(x_mod, xiHat_mod) * F.Transpose();

                    Vector<double>[] xTilde = new Vector<double>[Models.Count()];
                    Vector<double>[] zetaTilde = new Vector<double>[Models.Count()];
                    for (int i = 0; i < Models.Count(); i++)
                    {
                        xTilde[i] = F * xiHat_mod[i] + f;
                        zetaTilde[i] = Zeta(t, xTilde[i], y_mod[i], kTilde);
                    }

                    Matrix<double> CovZetaTilde = Exts.Cov(zetaTilde, zetaTilde);
                    Matrix<double> InvCovZetaTilde = Matrix<double>.Build.Dense(CovZetaTilde.RowCount, CovZetaTilde.ColumnCount, 0.0);
                    try
                    {
                        InvCovZetaTilde = CovZetaTilde.PseudoInverse();
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine("Can't inverse ZetaTilde");
                        Console.WriteLine(CovZetaTilde.ToString());
                        Console.WriteLine(e.Message);
                    }
                    Matrix<double> H = Exts.Cov(x_mod.Subtract(f), zetaTilde) * InvCovZetaTilde;
                    Vector<double> h = -H * zetaTilde.Average();

                    Matrix<double> kHat = kTilde - Exts.Cov(x_mod.Subtract(f), zetaTilde) * H.Transpose();

                    Vector<double> xTilde__ = F * Xi(t, xHat_) + f;
                    Vector<double> xHat__ = xTilde__ + H * Zeta(t, xTilde__, y, kTilde) + h;

                    return (xHat__, kHat);
                }
            }
        }
    }
}
