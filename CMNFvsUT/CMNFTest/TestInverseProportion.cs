using System;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using NonlinearSystem;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CMNFTest
{
    class TestInverseProportionGoodScalar : TestEnvironmentVector
    {
        public TestInverseProportionGoodScalar()
        {
            TestName = "Обратнопропорциональная зависимость (уст)";
            TestFileName = "InverseProportionGood";

            Vector<double> mW = Exts.Vector(0); Matrix<double> dW = Exts.Diag(10);
            Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(1);
            Vector<double> mEta = Exts.Vector(100.0); Matrix<double> dEta = Exts.Diag(100.0);
            Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(1.0 / (Math.Sign(x[0]) * Math.Pow(Math.Abs(x[0]), 1.0 / 3.0)));
            Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Diag(1.0);
            Func<int, Vector<double>, Vector<double>> psi = (s, x) => Exts.Vector(x[0]);

            Phi1_latex = new string[] { @"\frac{1}{\sqrt[3]{x_t}}" };
            Phi2_latex = new string[][] { new string[] { "1" } };
            Psi_latex = new string[] { @"x_t" };

            P_W = @"\mathcal{N}\left(" + mW.ToLatex() + ", " + mW.ToLatex() + @"\right)";
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
            Zeta = (s, x, y) => y - psi(s, x) - mNu;
            W = (s) => Exts.Vector(NormalW[0].Sample());
            Nu = (s) => Exts.Vector(NormalNu[0].Sample());
            DW = dW;
            DNu = dNu;
            X0 = () => Exts.Vector(NormalEta[0].Sample());
            X0Hat = mEta;
            DX0Hat = dEta;
        }
    }

    class TestInverseProportionBadScalar : TestEnvironmentVector
    {
        public TestInverseProportionBadScalar()
        {
            TestName = "Обратнопропорциональная зависимость (неуст)";
            TestFileName = "InverseProportionBad";

            Vector<double> mW = Exts.Vector(0); Matrix<double> dW = Exts.Diag(10);
            Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(1);
            Vector<double> mEta = Exts.Vector(100.0); Matrix<double> dEta = Exts.Diag(100.0);
            Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(Math.Min(1e5, 1.0 / Math.Pow(Math.Abs(x[0]), 2.0)));
            Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Diag(1.0);
            Func<int, Vector<double>, Vector<double>> psi = (s, x) => Exts.Vector(x[0]);

            Phi1_latex = new string[] { @"min(10^5,\frac{1}{x_t^2})" };
            Phi2_latex = new string[][] { new string[] { "1" } };
            Psi_latex = new string[] { @"x_t" };

            P_W = @"\mathcal{N}\left(" + mW.ToLatex() + ", " + mW.ToLatex() + @"\right)";
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
            Zeta = (s, x, y) => y - psi(s, x) - mNu;
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