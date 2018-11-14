using System;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using NonlinearSystem;

namespace TestEnvironments
{
    public class TestSwitchingObservations : TestEnvironmentVector
    {
        public TestSwitchingObservations(double _dnu)
        {
            {
                TestName = "Модель с переключающимися каналами наблюдений";
                TestFileName = "SwitchingObservations";

                double a = 0.8;
                double b = 0.2;
                double c = 6.0;
                //Vector<double> d = Exts.Vector(10.0, 1.0, 1.0);
                //Vector<double> sig = Exts.Vector(1.0, 1.0, 10.0);
                Vector<double> d = Exts.Vector(4.0, 1.0, 0.5);
                Vector<double> sig = Exts.Vector(1.0, 1.0, 3.0);
                double m = b / (1 - a);
                double S = Math.Sqrt(c * c / (1 - a * a));


                double l1 = -0.6745;
                double l2 = 0.6745;

                Func<double, int> I = x =>
                {
                    if (x < l1) return 0;
                    else if (x < l2) return 1;
                    else return 2;
                };

                //Vector<double> f = Exts.Vector(Normal.CDF(m, S, l1), Normal.CDF(m, S, l2) - Normal.CDF(m, S, l1), 1.0 - Normal.CDF(m, S, l2));
                Vector<double> f = Exts.Vector(0.25, 0.5, 0.25);



                Vector<double> mW = Exts.Vector(0,0.0); Matrix<double> dW = Exts.Diag(1.0, 1.0);
                Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(_dnu);
                Vector<double> mEta = Exts.Vector(0, 0.0); Matrix<double> dEta = Exts.Diag(1.0, 1.0); // FOR AIT
                //Vector<double> mEta = Exts.Vector(1.0, 0.0); Matrix<double> dEta = Exts.Diag(100.0, 1.0); // FOR IEOPR (no transit)
                Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(a * x[0] + b, 0);
                Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Diag(c, 1.0);
                Func<int, Vector<double>, Vector<double>> psi1 = (s, x) => Exts.Vector(d[I(x[1])] * x[0]);
                Func<int, Vector<double>, Matrix<double>> psi2 = (s, x) => Exts.Matrix(sig[I(x[1])]);

                Phi1_latex = new string[] { @"a x_t + b", "0" };
                Phi2_latex = new string[][] { new string[] { @"c", "0" }, new string[] { @"0", "1" }};
                Psi1_latex = new string[] { @"d^T * e(x^1_t) * x_t" };
                Psi2_latex = new string[][] { new string[] { @"\sigma^T * e(x^1_t)" } };

                P_W = @"\mathcal{N}\left(" + mW.ToLatex() + ", " + dW.ToLatex() + @"\right)";
                P_Nu = @"\mathcal{N}\left(" + mNu.ToLatex() + ", " + dNu.ToLatex() + @"\right)";
                P_Eta = @"\mathcal{N}\left(" + mEta.ToLatex() + ", " + dEta.ToLatex() + @"\right)";

                Normal[] NormalW = new Normal[2] { new Normal(mW[0], Math.Sqrt(dW[0, 0])), new Normal(mW[1], Math.Sqrt(dW[1, 1])) };
                Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])) };
                Normal[] NormalEta = new Normal[2] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])), new Normal(mEta[1], Math.Sqrt(dEta[1, 1])) };

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
                        double xi = k[0, 0] * d[i] / (d[i] * d[i] * k[0, 0] + sig[i] * sig[i]) * (y[0] - d[i] * x[0]);
                        num += f[i] * Normal.PDF(d[i] * x[0], Math.Sqrt(d[i] * d[i] * k[0, 0] + sig[i] * sig[i]), y[0]) * xi;
                        den += f[i] * Normal.PDF(d[i] * x[0], Math.Sqrt(d[i] * d[i] * k[0, 0] + sig[i] * sig[i]), y[0]);
                    }
                    return Exts.Vector(num / den);
                };

                W = (s) => Exts.Vector(NormalW[0].Sample(), NormalW[1].Sample());
                Nu = (s) => Exts.Vector(NormalNu[0].Sample());
                DW = dW;
                DNu = dNu;
                X0 = () => Exts.Vector(NormalEta[0].Sample(), NormalEta[1].Sample());
                X0Hat = mEta;
                DX0Hat = dEta;
            }
        }
    }

    public class TestSwitchingObservationsIdentification : TestEnvironmentVector
    {
        public TestSwitchingObservationsIdentification(double _dnu)
        {
            {
                TestName = "Модель с переключающимися каналами наблюдений и идентификацией";
                TestFileName = "SwitchingObservationsIdentification";

                double a = 0.8;
                double b = 0.2;
                double c = 6.0;
                //Vector<double> d = Exts.Vector(10.0, 1.0, 1.0);
                //Vector<double> sig = Exts.Vector(1.0, 1.0, 10.0);
                Vector<double> d = Exts.Vector(4.0, 1.0, 0.5);
                Vector<double> sig = Exts.Vector(1.0, 1.0, 3.0);
                double m = b / (1 - a);
                double S = Math.Sqrt(c * c / (1 - a * a));


                double l1 = -0.6745;
                double l2 = 0.6745;

                Func<double, int> I = x =>
                {
                    if (x < l1) return 0;
                    else if (x < l2) return 1;
                    else return 2;
                };

                //Vector<double> f = Exts.Vector(Normal.CDF(m, S, l1), Normal.CDF(m, S, l2) - Normal.CDF(m, S, l1), 1.0 - Normal.CDF(m, S, l2));
                Vector<double> f = Exts.Vector(0.25, 0.5, 0.25);



                Vector<double> mW = Exts.Vector(0, 0, a); Matrix<double> dW = Exts.Diag(1.0, 1.0, 0.04/3.0);
                Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(_dnu);
                Vector<double> mEta = Exts.Vector(0, 0.0, a); Matrix<double> dEta = Exts.Diag(1.0, 1.0, 0.04/3.0); 
                Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(x[2] * x[0] + b, 0, 0);
                Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Diag(c, 1.0, 1.0);
                Func<int, Vector<double>, Vector<double>> psi1 = (s, x) => Exts.Vector(d[I(x[1])] * x[0]);
                Func<int, Vector<double>, Matrix<double>> psi2 = (s, x) => Exts.Matrix(sig[I(x[1])]);

                Normal[] NormalW = new Normal[2] { new Normal(mW[0], Math.Sqrt(dW[0, 0])), new Normal(mW[1], Math.Sqrt(dW[1, 1])) };
                Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])) };
                Normal[] NormalEta = new Normal[2] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])), new Normal(mEta[1], Math.Sqrt(dEta[1, 1])) };
                ContinuousUniform UniformW = new ContinuousUniform(a - 0.2, a + 0.2);
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
                        double xi = k[0, 0] * d[i] / (d[i] * d[i] * k[0, 0] + sig[i] * sig[i]) * (y[0] - d[i] * x[0]);
                        num += f[i] * Normal.PDF(d[i] * x[0], Math.Sqrt(d[i] * d[i] * k[0, 0] + sig[i] * sig[i]), y[0]) * xi;
                        den += f[i] * Normal.PDF(d[i] * x[0], Math.Sqrt(d[i] * d[i] * k[0, 0] + sig[i] * sig[i]), y[0]);
                    }
                    return Exts.Vector(num / den);
                };

                W = (s) => Exts.Vector(NormalW[0].Sample(), NormalW[1].Sample(), UniformW.Sample());
                Nu = (s) => Exts.Vector(NormalNu[0].Sample());
                DW = dW;
                DNu = dNu;
                X0 = () => Exts.Vector(NormalEta[0].Sample(), NormalEta[1].Sample(), UniformW.Sample());
                X0Hat = mEta;
                DX0Hat = dEta;
            }
        }
    }

    public class AnotherTestSwitchingObservationsIdentification : TestEnvironmentVector
    {
        public AnotherTestSwitchingObservationsIdentification(double _dnu)
        {
            {
                TestName = "Еще одна модель с переключающимися каналами наблюдений и идентификацией";
                TestFileName = "AnotherSwitchingObservationsIdentification";

                Vector<double> sig = Exts.Vector(1.0, 4.0, 10.0);


                double l1 = -0.6745;
                double l2 = 0.6745;

                Func<double, int> I = x =>
                {
                    if (x < l1) return 0;
                    else if (x < l2) return 1;
                    else return 2;
                };

                //Vector<double> f = Exts.Vector(Normal.CDF(m, S, l1), Normal.CDF(m, S, l2) - Normal.CDF(m, S, l1), 1.0 - Normal.CDF(m, S, l2));
                Vector<double> f = Exts.Vector(0.25, 0.5, 0.25);



                Vector<double> mW = Exts.Vector(0, 0); Matrix<double> dW = Exts.Diag(0.0, 0.0);
                Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(_dnu);
                Vector<double> mEta = Exts.Vector(0, 0.0); Matrix<double> dEta = Exts.Diag(0.27, 0.0);
                Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(x[0], x[0]*x[1]);
                Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Diag(0.0, 0.0);
                Func<int, Vector<double>, Vector<double>> psi1 = (s, x) => Exts.Vector(x[1]);
                Func<int, Vector<double>, Matrix<double>> psi2 = (s, x) => Exts.Matrix(sig[I(x[1])]);

                Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])) };
                ContinuousUniform UniformW = new ContinuousUniform(-0.9, 0.9);
                Normal NormalEta = new Normal(mEta[0], dEta[0,0]);

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
                    for (int i = 0; i < sig.Count; i++)
                    {
                        double xi = k[0, 0] / (k[0, 0] + sig[i] * sig[i]) * (y[0] - x[1]);
                        num += f[i] * Normal.PDF(x[1], Math.Sqrt(k[0, 0] + sig[i] * sig[i]), y[0]) * xi;
                        den += f[i] * Normal.PDF(x[1], Math.Sqrt(k[0, 0] + sig[i] * sig[i]), y[0]);
                    }
                    return Exts.Vector(num / den);
                };

                W = (s) => Exts.Vector(0, 0);
                Nu = (s) => Exts.Vector(NormalNu[0].Sample());
                DW = dW;
                DNu = dNu;
                X0 = () => Exts.Vector(UniformW.Sample(), 1.0);
                X0Hat = mEta;
                DX0Hat = dEta;
            }
        }
    }

    public class YetAnotherTestSwitchingObservationsIdentification : TestEnvironmentVector
    {
        public YetAnotherTestSwitchingObservationsIdentification(double _dnu)
        {
            {
                TestName = "И еще одна модель с переключающимися каналами наблюдений и идентификацией";
                TestFileName = "YetAnotherSwitchingObservationsIdentification";

                Vector<double> sig = Exts.Vector(1.0, 4.0, 10.0);


                double l1 = -0.6745;
                double l2 = 0.6745;

                Func<double, int> I = x =>
                {
                    if (x < l1) return 0;
                    else if (x < l2) return 1;
                    else return 2;
                };

                Func<double, double> F = x =>
                {
                    if (x < -1) return -1;
                    else if (x > 1) return 1;
                    else return x;
                };

                //Vector<double> f = Exts.Vector(Normal.CDF(m, S, l1), Normal.CDF(m, S, l2) - Normal.CDF(m, S, l1), 1.0 - Normal.CDF(m, S, l2));
                Vector<double> f = Exts.Vector(0.25, 0.5, 0.25);



                Vector<double> mW = Exts.Vector(0, 0); Matrix<double> dW = Exts.Diag(0.0, 0.0);
                Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(_dnu);
                Vector<double> mEta = Exts.Vector(0, 0.0); Matrix<double> dEta = Exts.Diag(0.27, 0.0);
                Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(x[0], x[0] * F(x[1]));
                Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Diag(0.0, 0.0);
                Func<int, Vector<double>, Vector<double>> psi1 = (s, x) => Exts.Vector(x[1]);
                Func<int, Vector<double>, Matrix<double>> psi2 = (s, x) => Exts.Matrix(sig[I(x[1])]);
                
                Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])) };
                ContinuousUniform UniformW = new ContinuousUniform(-0.9, 0.9);
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
                    for (int i = 0; i < sig.Count; i++)
                    {
                        double xi = k[0, 0] / (k[0, 0] + sig[i] * sig[i]) * (y[0] - x[1]);
                        num += f[i] * Normal.PDF(x[1], Math.Sqrt(k[0, 0] + sig[i] * sig[i]), y[0]) * xi;
                        den += f[i] * Normal.PDF(x[1], Math.Sqrt(k[0, 0] + sig[i] * sig[i]), y[0]);
                    }
                    return Exts.Vector(num / den);
                };

                W = (s) => Exts.Vector(0, 0);
                Nu = (s) => Exts.Vector(NormalNu[0].Sample());
                DW = dW;
                DNu = dNu;
                X0 = () => Exts.Vector(UniformW.Sample(), 1.0);
                X0Hat = mEta;
                DX0Hat = dEta;
            }
        }
    }

    public class HopefullyTheLastTestSwitchingObservationsIdentification : TestEnvironmentVector
    {
        public HopefullyTheLastTestSwitchingObservationsIdentification(double _dnu)
        {
            {
                TestName = "Надеюсь, последняя модель с переключающимися каналами наблюдений и идентификацией";
                TestFileName = "HopefullyTheLastSwitchingObservationsIdentification";

                Vector<double> sig = Exts.Vector(1.0, 4.0, 10.0);


                double l1 = -0.6745;
                double l2 = 0.6745;

                Func<double, int> I = x =>
                {
                    if (x < l1) return 0;
                    else if (x < l2) return 1;
                    else return 2;
                };

                Func<double, double> F = x =>
                {
                    return 0.4 - x * Math.Exp(-0.3 * x * x);
                };

                //Vector<double> f = Exts.Vector(Normal.CDF(m, S, l1), Normal.CDF(m, S, l2) - Normal.CDF(m, S, l1), 1.0 - Normal.CDF(m, S, l2));
                Vector<double> f = Exts.Vector(0.25, 0.5, 0.25);



                Vector<double> mW = Exts.Vector(0, 0); Matrix<double> dW = Exts.Diag(0.0, 0.0);
                Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(_dnu);
                Vector<double> mEta = Exts.Vector(0, 0.0); Matrix<double> dEta = Exts.Diag(0.27, 0.0);
                Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(x[0], x[0] * F(x[1]));
                Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Diag(0.0, 0.0);
                Func<int, Vector<double>, Vector<double>> psi1 = (s, x) => Exts.Vector(x[1]);
                Func<int, Vector<double>, Matrix<double>> psi2 = (s, x) => Exts.Matrix(sig[I(x[1])]);
                
                Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])) };
                ContinuousUniform UniformW = new ContinuousUniform(-0.9, 0.9);
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
                    for (int i = 0; i < sig.Count; i++)
                    {
                        double xi = k[0, 0] / (k[0, 0] + sig[i] * sig[i]) * (y[0] - x[1]);
                        num += f[i] * Normal.PDF(x[1], Math.Sqrt(k[0, 0] + sig[i] * sig[i]), y[0]) * xi;
                        den += f[i] * Normal.PDF(x[1], Math.Sqrt(k[0, 0] + sig[i] * sig[i]), y[0]);
                    }
                    return Exts.Vector(num / den);
                };

                W = (s) => Exts.Vector(0, 0);
                Nu = (s) => Exts.Vector(NormalNu[0].Sample());
                DW = dW;
                DNu = dNu;
                X0 = () => Exts.Vector(UniformW.Sample(), 1.0);
                X0Hat = mEta;
                DX0Hat = dEta;
            }
        }
    }
}

