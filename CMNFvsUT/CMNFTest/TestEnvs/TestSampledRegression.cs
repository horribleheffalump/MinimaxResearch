using System;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using NonlinearSystem;

namespace CMNFTest
{
    class TestSampledRegression : TestEnvironmentVector
    {
        public TestSampledRegression(double _dnu)
        {
            {
                TestName = "Модель семплированной регрессии";
                TestFileName = "SampledRegression";

                Vector<double> a = Exts.Vector(0.3, 0.4, 0.7);
                Vector<double> b = Exts.Vector(1.4, 3.0, 3.0);
                Vector<double> c = Exts.Vector(0.9, 1.5, 2.5);
                Vector<double> d = Exts.Vector(0.33, 0.37, 0.3);
                Vector<double> m = Exts.Vector(b[0] / (1 - a[0]), b[1] / (1 - a[1]), b[2] / (1 - a[2]));
                Vector<double> S = Exts.Vector(c[0] / Math.Sqrt(1 - a[0] * a[0]), c[1] / Math.Sqrt(1 - a[1] * a[1]), c[2] / Math.Sqrt(1 - a[2] * a[2]));

                Func<double, int> I = x =>
                {
                    if (x < 3) return 0;
                    else if (x < 7) return 1;
                    else return 2;
                };

                Vector<double> mW = Exts.Vector(0); Matrix<double> dW = Exts.Diag(1.0);
                Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(_dnu);
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
                //Xi = (s, x) =>
                //{
                //    double num = 0;
                //    double den = 0;
                //    for (int i = 0; i < a.Count; i++)
                //    {
                //        num += d[i] * Normal.PDF(m[i], S[i], x[0]) * (a[i] * x[0] + b[i]);
                //        den += d[i] * Normal.PDF(m[i], S[i], x[0]);
                //    }
                //    return Exts.Vector(num / den);
                //};
                Zeta = (s, x, y, k) => y - psi(s, x) - mNu;
                W = (s) => Exts.Vector(NormalW[0].Sample());
                Nu = (s) => Exts.Vector(NormalNu[0].Sample());
                DW = dW;
                DNu = dNu;
                X0 = () => Exts.Vector(NormalEta[0].Sample());
                X0Hat = mEta;
                DX0Hat = dEta;

                useSimpleModel = false;
            }
        }
    }

    class TestSampledRegressionObs : TestEnvironmentVector
    {
        public TestSampledRegressionObs(double _dnu)
        {
            {
                TestName = "Модель семплированной регрессии";
                TestFileName = "SampledRegression";

                double a = 0.8;
                double b = 0.2;
                double c = 6.0;
                Vector<double> d = Exts.Vector(10.0, 1.0, 1.0);
                Vector<double> sig = Exts.Vector(1.0, 1.0, 10.0);
                double m = b / (1 - a);
                double S = Math.Sqrt(c * c / (1 - a * a));


                double l1 = m - S/2;
                double l2 = m + S/2;

                Vector<double> f = Exts.Vector(Normal.CDF(m, S, l1), Normal.CDF(m, S, l2) - Normal.CDF(m, S, l1), 1.0 - Normal.CDF(m, S, l2));

                Func<double, int> I = x =>
                {
                    if (x < l1) return 0;
                    else if (x < l2) return 1;
                    else return 2;
                };

                Vector<double> mW = Exts.Vector(0); Matrix<double> dW = Exts.Diag(1.0);
                Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(_dnu);
                Vector<double> mEta = Exts.Vector(0); Matrix<double> dEta = Exts.Diag(1.0);
                Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(a * x[0] + b);
                Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Matrix(c);
                Func<int, Vector<double>, Vector<double>> psi1 = (s, x) => Exts.Vector(d[I(x[0])] * x[0]);
                Func<int, Vector<double>, Matrix<double>> psi2 = (s, x) => Exts.Matrix(sig[I(x[0])]);

                Phi1_latex = new string[] { @"a x_t + b" };
                Phi2_latex = new string[][] { new string[] { @"c" } };
                Psi1_latex = new string[] { @"d^T * e(x_t) * x_t" };
                Psi2_latex = new string[][] { new string[] { @"\sigma^T * e(x_t)" } };

                P_W = @"\mathcal{N}\left(" + mW.ToLatex() + ", " + dW.ToLatex() + @"\right)";
                P_Nu = @"\mathcal{N}\left(" + mNu.ToLatex() + ", " + dNu.ToLatex() + @"\right)";
                P_Eta = @"\mathcal{N}\left(" + mEta.ToLatex() + ", " + dEta.ToLatex() + @"\right)";

                Normal[] NormalW = new Normal[1] { new Normal(mW[0], Math.Sqrt(dW[0, 0])) };
                Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])) };
                Normal[] NormalEta = new Normal[1] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])) };

                //Expression<Func<int, Vector<double>, Vector<double>>> expr = (s, x) => Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1])); ;

                Phi1 = phi1;
                Phi2 = phi2;
                Psi1 = psi1;
                Psi2 = psi2;
                Xi = (s, x) => phi1(s, x) + phi2(s, x) * mW;
                //Zeta = (s, x, y, k) => y - psi1(s, x) - psi2(s, x) * mNu;

                Zeta = (s, x, y, k) =>
                {
                    double num = 0;
                    double den = 0;
                    for (int i = 0; i < d.Count; i++)
                    {
                        num += f[i] * Normal.PDF(y[0], d[i] * x[0], d[i] * d[i] * k[0,0] + sig[i] * sig[i]) * x[0];
                        den += d[i] * Normal.PDF(y[0], d[i] * x[0], d[i] * d[i] * k[0, 0] + sig[i] * sig[i]);
                    }
                    return Exts.Vector(num / den);
                };


                W = (s) => Exts.Vector(NormalW[0].Sample());
                Nu = (s) => Exts.Vector(NormalNu[0].Sample());
                DW = dW;
                DNu = dNu;
                X0 = () => Exts.Vector(NormalEta[0].Sample());
                X0Hat = mEta;
                DX0Hat = dEta;

                useSimpleModel = false;
            }
        }
    }
}

