using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestEnvironments
{
    public class TestSimpleIdentification : TestEnvironmentVector
    {
        public TestSimpleIdentification()
        {
            TestName = "Простая модель с идентификацией";
            TestFileName = "TestSimpleIdentification";


            Vector<double> mW = Exts.Vector(0, 0); Matrix<double> dW = Exts.Diag(0, 4.0);
            Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(0.0004);
            Vector<double> mEta = Exts.Vector(0, 0.0); Matrix<double> dEta = Exts.Diag(0.27, 0);
            Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Exts.Vector(x[0], x[0] * x[1]);
            Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Exts.Diag(0.0, 1.0);
            Func<int, Vector<double>, Vector<double>> psi1 = (s, x) => Exts.Vector(x[1]);
            Func<int, Vector<double>, Matrix<double>> psi2 = (s, x) => Exts.Matrix(1.0);

            RandomVector<Normal> NormalW = new RandomVector<Normal>(mW, dW);
            RandomVector<Normal> NormalNu = new RandomVector<Normal>(mNu, dNu);
            ContinuousUniform UniformW = new ContinuousUniform(-0.9, 0.9);

            Phi1 = phi1;
            Phi2 = phi2;
            Psi1 = psi1;
            Psi2 = psi2;
            Xi = (s, x) => phi1(s, x) + phi2(s, x) * mW;
            Zeta = (s, x, y, k) => y - psi1(s, x) - psi2(s, x) * mNu;


            W = (s) => NormalW.Sample();
            Nu = (s) => NormalNu.Sample();
            DW = dW;
            DNu = dNu;
            X0 = () => Exts.Vector(UniformW.Sample(), 0.0);
            X0Hat = mEta;
            DX0Hat = dEta;
        }
    }
}
