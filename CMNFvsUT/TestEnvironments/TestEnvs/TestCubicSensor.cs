using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using NonlinearSystem;
using PythonInteract;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TestEnvironments
{
    public class TestCubicSensor : TestEnvironmentVector
    {
        public TestCubicSensor(double _dw, double _dnu)
        {
            TestName = "Кубический сенсор";
            TestFileName = "CubicSensor";

            Vector<double> mW = Exts.Vector(0, 0); Matrix<double> dW = Exts.Diag(_dw, _dw);
            Vector<double> mNu = Exts.Vector(0, 0); Matrix<double> dNu = Exts.Diag(_dnu, _dnu);
            //Vector<double> mEta = Exts.Vector(100, 100); Matrix<double> dEta = Exts.Diag(100, 100);
            Vector<double> mEta = Exts.Vector(0, 0); Matrix<double> dEta = Exts.Diag(1, 1);
            Func<int, Vector<double>, Vector<double>> phi = (s, x) => Exts.Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1]));
            Func<int, Vector<double>, Vector<double>> psi = (s, x) => Exts.Vector(Math.Pow(x[0], 3) + Math.Pow(x[0], 1), Math.Pow(x[1], 3) + Math.Pow(x[1], 1));
            Func<int, Vector<double>, Matrix<double>> psi_test = (s, x) => Matrix<double>.Build.Dense(1, 1, 1.0);
            Psi2 = psi_test;

            //Phi1_latex = new string[] { @"\frac{x_0}{1 + x_0^2}", @"\frac{x_1}{1 + x_1^2}" };
            //Psi1_latex = new string[] { @"x_0^3+x_0", @"x_1^3+x_1" };

            //P_W = @"\mathcal{N}\left(\mathbf{0}, \mathbf{E}\right)";
            //P_Nu = @"\mathcal{N}\left(\mathbf{0}, \mathbf{E}\right)";
            //P_Eta = @"\mathcal{N}\left(" + mEta.ToLatex() + ", " + dEta.ToLatex() + @"\right)";

            Normal[] NormalW = new Normal[2] { new Normal(mW[0], Math.Sqrt(dW[0, 0])), new Normal(mW[1], Math.Sqrt(dW[1, 1])) };
            Normal[] NormalNu = new Normal[2] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])), new Normal(mNu[1], Math.Sqrt(dNu[1, 1])) }; ;
            Normal[] NormalEta = new Normal[2] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])), new Normal(mEta[1], Math.Sqrt(dEta[1, 1])) }; ;

            //Expression<Func<int, Vector<double>, Vector<double>>> expr = (s, x) => Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1])); ;

            Phi1 = phi;
            Psi1 = psi;
            Xi = (s, x) => phi(s, x) + mW;
            Zeta = (s, x, y, k) => y - psi(s, x) - mNu;
            W = (s) => Exts.Vector(NormalW[0].Sample(), NormalW[1].Sample());
            Nu = (s) => Exts.Vector(NormalNu[0].Sample(), NormalNu[1].Sample());
            DW = dW;
            DNu = dNu;
            X0 = () => Exts.Vector(NormalEta[0].Sample(), NormalEta[1].Sample());
            X0Hat = mEta;
            DX0Hat = dEta;
        }
    }
    public class TestCubicSensorScalar : TestEnvironmentVector
    {
        public TestCubicSensorScalar(double _dw, double _dnu)
        {
            TestName = "Кубический сенсор";
            TestFileName = "CubicSensor";

            Vector<double> mW = Exts.Vector(0); Matrix<double> dW = Exts.Diag(_dw);
            Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(_dnu);
            Vector<double> mEta = Exts.Vector(0.1); Matrix<double> dEta = Exts.Diag(1); // FOR AIT
            //Vector<double> mEta = Exts.Vector(0.1); Matrix<double> dEta = Exts.Diag(1.16); // FOR IEOPR
            Func<int, Vector<double>, Vector<double>> phi = (s, x) => Exts.Vector(x[0] / (1 + x[0] * x[0]));
            Func<int, Vector<double>, Vector<double>> psi = (s, x) => Exts.Vector(Math.Pow(x[0], 3) + Math.Pow(x[0], 1));

            Func<int, Vector<double>, Matrix<double>> dpsi = (s, x) => Exts.Matrix(3.0 * Math.Pow(x[0], 2) + 1.0);

            //Phi1_latex = new string[] { @"\frac{x_t}{1 + x_t^2}"};
            //Psi1_latex = new string[] { @"x_t^3+x_t"};

            //P_W = @"\mathcal{N}\left(" + mW.ToLatex() + ", " + mW.ToLatex() + @"\right)";
            //P_Nu = @"\mathcal{N}\left(" + mNu.ToLatex() + ", " + dNu.ToLatex() + @"\right)";
            //P_Eta = @"\mathcal{N}\left(" + mEta.ToLatex() + ", " + dEta.ToLatex() + @"\right)";

            Normal[] NormalW = new Normal[1] { new Normal(mW[0], Math.Sqrt(dW[0, 0])) };
            Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0]))}; 
            Normal[] NormalEta = new Normal[1] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0]))};

            //Expression<Func<int, Vector<double>, Vector<double>>> expr = (s, x) => Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1])); ;

            Phi1 = phi;
            Psi1 = psi;

            dPhi = (s, x) => Exts.Diag(0.9);
            dPsi = (s, x) => Exts.Diag(1.0);

            Xi = (s, x) => phi(s, x) + mW;
            Zeta = (s, x, y, k) => y - psi(s, x) - mNu;
            //Zeta = (s, x, y, k) => k * dpsi(s,x).Transpose() * (dpsi(s,x) * k * dpsi(s,x).Transpose() + dNu ).PseudoInverse() * (y - psi(s, x) - mNu);
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
