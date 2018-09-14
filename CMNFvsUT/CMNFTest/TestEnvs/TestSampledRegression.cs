using System;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using NonlinearSystem;

namespace CMNFTest
{
    class TestSampledRegression : TestEnvironmentVector
    {
        public TestSampledRegression()
        {
            {
                TestName = "Модель семплированной регрессии";
                TestFileName = "SampledRegression";

                Vector<double> a = Exts.Vector(0.3, 0.4, 0.7);
                Vector<double> b = Exts.Vector(1.4, 3.0, 3.0);
                Vector<double> c = Exts.Vector(0.9, 1.5, 2.5);

                Func<double, int> I = x =>
                {
                    if (x < 3) return 0;
                    else if (x < 7) return 1;
                    else return 2;
                };

                Vector<double> mW = Exts.Vector(0); Matrix<double> dW = Exts.Diag(1.0);
                Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(1.0);
                Vector<double> mEta = Exts.Vector(0); Matrix<double> dEta = Exts.Diag(1.0);
                Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(a[I(x[0])] * x[0] + b[I(x[0])]);
                Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Matrix(c[I(x[0])]);
                Func<int, Vector<double>, Vector<double>> psi = (s, x) => Exts.Vector(x[0]);

                Phi1_latex = new string[] { @"a^T e(x_t) x_t + b^T e(x_t)" };
                Phi2_latex = new string[][] { new string[] { @"c^T e(x_t)" } };
                Psi1_latex = new string[] { @"x_t" };

                P_W = @"\mathcal{N}\left(" + mW.ToLatex() + ", " + dW.ToLatex() + @"\right)";
                P_Nu = @"\mathcal{N}\left(" + mNu.ToLatex() + ", " + dNu.ToLatex() + @"\right)";
                P_Eta = @"\mathcal{N}\left(" + mEta.ToLatex() + ", " + dEta.ToLatex() + @"\right)";

                Normal[] NormalW = new Normal[1] { new Normal(mW[0], Math.Sqrt(dW[0, 0])) };
                Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])) };
                Normal[] NormalEta = new Normal[1] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])) };

                //Expression<Func<int, Vector<double>, Vector<double>>> expr = (s, x) => Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1])); ;

                Phi1 = phi1;
                Phi2 = phi2;
                Psi1 = psi;
                Xi = (s, x) => phi1(s, x) + phi2(s, x) * mW;
                Zeta = (s, x, y, k) => y - psi(s, x) - mNu;
                W = (s) => Exts.Vector(NormalW[0].Sample());
                Nu = (s) => Exts.Vector(NormalNu[0].Sample());
                DW = dW;
                DNu = dNu;
                X0 = () => Exts.Vector(NormalEta[0].Sample());
                X0Hat = mEta;
                DX0Hat = dEta;
            }
        }
    }
}

