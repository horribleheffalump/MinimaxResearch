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
    public class ModifiedCMNFilter
    {
        private Func<int, Vector<double>, Vector<double>> Phi1; // Phi1(t, X)
        private Func<int, Vector<double>, Matrix<double>> Phi2;
        private Func<int, Vector<double>, Vector<double>> Psi1;
        private Func<int, Vector<double>, Matrix<double>> Psi2;
        private Func<int, Vector<double>> W;
        private Func<int, Vector<double>> Nu;


        Func<int, Vector<double>, Vector<double>> Xi;
        Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> Zeta;

        public ModifiedCMNFilter(Func<int, Vector<double>, Vector<double>> xi,
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
        }

        public (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat_, Matrix<double> kHat_, int n)
        {
            Vector<double>[] x_mod = new Vector<double>[n];
            Vector<double>[] y_mod = new Vector<double>[n];
            Parallel.For(0, n, new ParallelOptions() {MaxDegreeOfParallelism = System.Environment.ProcessorCount }, i =>
            {
                double[] x_arr = new double[xHat_.Count];
                for (int j = 0; j < xHat_.Count; j++)
                {
                    x_arr[j] = Normal.Sample(xHat_[j], Math.Sqrt(kHat_[j, j]));
                }
                x_mod[i] = Exts.Vector(x_arr);

                if (t > 0)
                {
                    x_mod[i] = Phi1(t, x_mod[i]) + Phi2(t, x_mod[i]) * W(t);
                }
                y_mod[i] = Psi1(t, x_mod[i]) + Psi2(t, x_mod[i]) * Nu(t);
            });

            Vector<double> f = x_mod.Average(); 
            Matrix<double> kTilde = Exts.Cov(x_mod, x_mod);

            Vector<double>[] zetaTilde = new Vector<double>[n];
            for (int i = 0; i < n; i++)
            {
                zetaTilde[i] = Zeta(t, f, y_mod[i], kTilde);
            }

            Matrix<double> H = Exts.Cov(x_mod.Subtract(f), zetaTilde) * (Exts.Cov(zetaTilde, zetaTilde).PseudoInverse());
            Vector<double> h = -H * zetaTilde.Average();

            Matrix<double> kHat = kTilde - Exts.Cov(x_mod.Subtract(f), zetaTilde) * H.Transpose();

            Vector<double> xHat__ = f + H * Zeta(t, f, y, kTilde) + h;
            
            return (xHat__, kHat);
        }
    }
}
