using System;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using NonlinearSystem;

namespace CMNFTest
{
    class TestLogisticModelScalar : TestEnvironmentVector
    {
        public TestLogisticModelScalar(double bound, double _dw, double _dnu)
        {
            TestName = "Логистическая модель";
            TestFileName = "LogisticModel";

            Vector<double> mW = Exts.Vector(0); Matrix<double> dW = Exts.Diag(_dw);
            Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(_dnu);
            Vector<double> mEta = Exts.Vector(0.1); Matrix<double> dEta = Exts.Diag(0.01); // small values are for regularization
            Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(Math.Max(-bound, Math.Min(bound, x[0] * (1 - x[0]))));
            Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Diag(1.0);
            Func<int, Vector<double>, Vector<double>> psi = (s, x) => Exts.Vector(x[0]);

            Phi1_latex = new string[] { @"max(-" + bound.ToString() + @", min(" + bound.ToString() + @",x(1-x)))" };
            Phi2_latex = new string[][] { new string[] { "1" } };
            Psi_latex = new string[] { @"x_t" };

            P_W = @"\mathcal{N}\left(" + mW.ToLatex() + ", " + dW.ToLatex() + @"\right)";
            P_Nu = @"\mathcal{N}\left(" + mNu.ToLatex() + ", " + dNu.ToLatex() + @"\right)";
            P_Eta = @"\mathcal{N}\left(" + mEta.ToLatex() + ", " + dEta.ToLatex() + @"\right)";

            Normal[] NormalW = new Normal[1] { new Normal(mW[0], Math.Sqrt(dW[0, 0])) };
            Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])) };
            Normal[] NormalEta = new Normal[1] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])) };

            //Expression<Func<int, Vector<double>, Vector<double>>> expr = (s, x) => Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1])); ;

            Phi1 = phi1;
            Phi2 = phi2;
            Psi = psi;
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

    class TestLogisticModelZeroScalar : TestEnvironmentVector
    {
        public TestLogisticModelZeroScalar(double bound, double _dw, double _dnu)
        {
            TestName = "Логистическая модель с возвратом";
            TestFileName = "LogisticModelZero";

            Vector<double> mW = Exts.Vector(0); Matrix<double> dW = Exts.Diag(_dw);
            Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(_dnu);
            Vector<double> mEta = Exts.Vector(0.1); Matrix<double> dEta = Exts.Diag(0.01); // small values are for regularization
            Func<int, Vector<double>, Vector<double>> phi1 = (s, x) =>
            {
                if (Math.Abs(x[0] * (1 - x[0])) < bound)
                    return Exts.Vector(x[0] * (1 - x[0]));
                else
                    return Exts.Vector(0);
            };
            Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Diag(1.0);
            Func<int, Vector<double>, Vector<double>> psi = (s, x) => Exts.Vector(x[0]);

            Phi1_latex = new string[] { @"x(1-x) if abs(x(1-x)) < " + bound.ToString() +"; 0 else" };
            Phi2_latex = new string[][] { new string[] { "1" } };
            Psi_latex = new string[] { @"x_t" };

            P_W = @"\mathcal{N}\left(" + mW.ToLatex() + ", " + dW.ToLatex() + @"\right)";
            P_Nu = @"\mathcal{N}\left(" + mNu.ToLatex() + ", " + dNu.ToLatex() + @"\right)";
            P_Eta = @"\mathcal{N}\left(" + mEta.ToLatex() + ", " + dEta.ToLatex() + @"\right)";

            Normal[] NormalW = new Normal[1] { new Normal(mW[0], Math.Sqrt(dW[0, 0])) };
            Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])) };
            Normal[] NormalEta = new Normal[1] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])) };

            //Expression<Func<int, Vector<double>, Vector<double>>> expr = (s, x) => Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1])); ;

            Phi1 = phi1;
            Phi2 = phi2;
            Psi = psi;
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


    class TestLogisticModelUniformNoiseScalar : TestEnvironmentVector
    {
        public TestLogisticModelUniformNoiseScalar()
        {
            TestName = "Логистическая модель с равномерным шумом";
            TestFileName = "LogisticModelUniform";

            Vector<double> mW = Exts.Vector(1e-5); Matrix<double> dW = Exts.Diag(1.0 / 3.0);
            Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(1);
            Vector<double> mEta = Exts.Vector(0.5); Matrix<double> dEta = Exts.Diag(0.01); // small values are for regularization
            Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(3.0 * x[0] * (1 - x[0]));
            Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Diag(x[0] * (1 - x[0]));
            Func<int, Vector<double>, Vector<double>> psi = (s, x) => Exts.Vector(x[0]);

            Phi1_latex = new string[] { @"???" };
            Phi2_latex = new string[][] { new string[] { "1" } };
            Psi_latex = new string[] { @"x_t" };

            P_W = @"\mathcal{R}\left(0,1\right)";
            P_Nu = @"\mathcal{N}\left(" + mNu.ToLatex() + ", " + dNu.ToLatex() + @"\right)";
            P_Eta = @"\mathcal{N}\left(" + mEta.ToLatex() + ", " + dEta.ToLatex() + @"\right)";

            ContinuousUniform[] UniformW = new ContinuousUniform[1] { new ContinuousUniform(-1 + 1e-5, 1 + 1e-5) };
            Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])) };
            Normal[] NormalEta = new Normal[1] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])) };

            //Expression<Func<int, Vector<double>, Vector<double>>> expr = (s, x) => Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1])); ;

            Phi1 = phi1;
            Phi2 = phi2;
            Psi = psi;
            Xi = (s, x) => phi1(s, x) + phi2(s, x) * mW;
            Zeta = (s, x, y, k) => y - psi(s, x) - mNu;
            W = (s) => Exts.Vector(UniformW[0].Sample());
            Nu = (s) => Exts.Vector(NormalNu[0].Sample());
            DW = dW;
            DNu = dNu;
            X0 = () => Exts.Vector(NormalEta[0].Sample());
            X0Hat = mEta;
            DX0Hat = dEta;
        }
    }

}