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
        private Matrix<double> DNu;
        private Vector<double> SigmaNu;

        private DiscreteVectorModel[] Models;
        private int N;


        Func<int, Vector<double>, Vector<double>> Xi;
        Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> Zeta;

        Vector<double>[] xHat;

        public AdaptiveCMNFilter(Func<int, Vector<double>, Vector<double>> xi,
                                    Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> zeta,
                                    Func<int, Vector<double>, Vector<double>> phi1,
                                    Func<int, Vector<double>, Matrix<double>> phi2,
                                    Func<int, Vector<double>, Vector<double>> psi1,
                                    Func<int, Vector<double>, Matrix<double>> psi2,
                                    Func<int, Vector<double>> w,
                                    Func<int, Vector<double>> nu,
                                    Matrix<double> dNu)
        {
            Xi = xi;
            Zeta = zeta;
            Phi1 = phi1;
            Phi2 = phi2;
            Psi1 = psi1;
            Psi2 = psi2;
            W = w;
            Nu = nu;
            DNu = dNu;
            SigmaNu = 10.0 * dNu.Diagonal().PointwiseSqrt();
        }

        public void Initialize(int n, Func<Vector<double>> X0, Vector<double> XHat0)
        {
            N = n;
            Models = new DiscreteVectorModel[N];
            for (int i = 0; i < N; i++)
            {
                Models[i] = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
                Models[i].Step();
            }
            xHat = Enumerable.Repeat(XHat0, n).ToArray();
        }

        public (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat_, Matrix<double> kHat_)
        {
            Vector<double>[] x_mod = new Vector<double>[Models.Count()];
            Vector<double>[] y_mod = new Vector<double>[Models.Count()];
            Vector<double>[] xiHat_mod = new Vector<double>[Models.Count()];

            RandomVector<Normal> xHatDistr = new RandomVector<Normal>(xHat_, kHat_);

            int selected = 0;

            for (int i = 0; i < Models.Count(); i++)
            {
                //Models[i].Step();                
                //if (Models[i].Trajectory[t][1].InBounds(y - SigmaNu, y + SigmaNu))
                //{
                //    x_mod[i] = Models[i].Trajectory[t][0];
                //    y_mod[i] = Models[i].Trajectory[t][1];
                //    xiHat_mod[i] = Xi(t, xHat[i]);
                //    selected++;
                //}
                //else
                //{
                    x_mod[i] = xHatDistr.Sample();
                    if (t > 0)
                    {
                        x_mod[i] = Phi1(t, x_mod[i]) + Phi2(t, x_mod[i]) * W(t);
                    }
                    y_mod[i] = Psi1(t, x_mod[i]) + Psi2(t, x_mod[i]) * Nu(t);
                    xiHat_mod[i] = Xi(t, xHatDistr.Sample());
                //    Models[i].Trajectory[t][0] = x_mod[i];
                //    Models[i].Trajectory[t][1] = y_mod[i];
                //}
            }

            int good = 0;

            for (int i = 0; i < Models.Count(); i++)
            {
                if (y_mod[i].InBounds(y - SigmaNu, y + SigmaNu))
                {
                    good++;
                }
            }
            //Console.WriteLine($"selected: {selected}, total good: {good}");

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

            for (int i = 0; i < Models.Count(); i++)
            {
                xHat[i] = xTilde[i] + H * zetaTilde[i] + h;
            }

            return (xHat__, kHat);
        }
    }
}
