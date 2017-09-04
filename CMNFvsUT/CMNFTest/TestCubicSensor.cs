using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using PythonInteract;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CMNFTest
{
    public class TestCubicSensor : TestEnvironmentVector
    {
        public TestCubicSensor()
        {
            TestName = "Кубический сенсор";
            Vector <double> mW = Extensions.Vector(0, 0); Matrix<double> dW = Extensions.Diag(1, 1);
            Vector<double> mNu = Extensions.Vector(0, 0); Matrix<double> dNu = Extensions.Diag(1, 1);
            Vector<double> mEta = Extensions.Vector(100, 100); Matrix<double> dEta = Extensions.Diag(100, 100);
            Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Extensions.Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1]));
            Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Extensions.Diag(1.0, 1.0);
            Func<int, Vector<double>, Vector<double>> psi = (s, x) => Extensions.Vector(Math.Pow(x[0], 3) + Math.Pow(x[0], 1), Math.Pow(x[1], 3) + Math.Pow(x[1], 1));

            Phi1_latex = @"\frac{x_0}{1 + x_0^2}, \frac{x_1}{1 + x_1^2}";
            Phi2_latex = @"Diag(1.0, 1.0)";
            Psi_latex = @"x_0^3+x_0, x_1^3+x_1";

            Normal[] NormalW = new Normal[2] { new Normal(mW[0], Math.Sqrt(dW[0, 0])), new Normal(mW[1], Math.Sqrt(dW[1, 1])) };
            Normal[] NormalNu = new Normal[2] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])), new Normal(mNu[1], Math.Sqrt(dNu[1, 1])) }; ;
            Normal[] NormalEta = new Normal[2] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])), new Normal(mEta[1], Math.Sqrt(dEta[1, 1])) }; ;

            //Expression<Func<int, Vector<double>, Vector<double>>> expr = (s, x) => Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1])); ;

            Phi1 = phi1;
            Phi2 = phi2;
            Psi = psi;
            Xi = (s, x) => phi1(s, x) + phi2(s, x) * mW;
            Zeta = (s, x, y) => y - psi(s, x) - mNu;
            W = (s) => Extensions.Vector(NormalW[0].Sample(), NormalW[1].Sample());
            Nu = (s) => Extensions.Vector(NormalNu[0].Sample(), NormalNu[1].Sample());
            DW = dW;
            DNu = dNu;
            X0 = () => Extensions.Vector(NormalEta[0].Sample(), NormalEta[1].Sample());
            X0Hat = mEta;
            DX0Hat = dEta;
        }
    }
}
