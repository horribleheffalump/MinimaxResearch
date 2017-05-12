using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using CMNF;
using CMNFTest.Properties;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using UKF;

namespace CMNFTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //int T = 25;
            //int N = 1000;
            //double mW = 0; double dW = 1;
            //double mNu = 0; double dNu = 1;
            //double mEta = 100; double dEta = 100;
            //Func<double, double> phi1 = new Func<double, double>(x => x / (1 + x * x));
            //Func<double, double> phi2 = new Func<double, double>(x => 1.0);
            //Func<double, double> psi = new Func<double, double>(x => Math.Pow(x, 3) + Math.Pow(x, 1));

            //Normal NormalW = new Normal(mW, Math.Sqrt(dW));
            //Normal NormalNu = new Normal(mNu, Math.Sqrt(dNu));
            //Normal NormalEta = new Normal(mEta, Math.Sqrt(dEta));

            //TestEnvironmentVector test1 = new TestEnvironmentVector()
            //{
            //    Phi1 = (s, x) => Vector(phi1(x[0])),
            //    Phi2 = (s, x) => Matrix(phi2(x[0])),
            //    Psi = (s, x) => Vector(psi(x[0])),
            //    Xi = (s, x) => Vector(phi1(x[0]) + phi2(x[0]) * mW),
            //    Zeta = (s, x, y) => Vector(y[0] - psi(x[0]) - mNu),
            //    W = (s) => Vector(NormalW.Sample()),
            //    Nu = (s) => Vector(NormalNu.Sample()),
            //    DW = Matrix(dW),
            //    DNu = Matrix(dNu),
            //    X0 = () => Vector(NormalEta.Sample()),
            //    X0Hat = Vector(mEta),
            //    DX0Hat = Matrix(dEta)
            //};

            //test1.Initialize(true, T, N);
            //test1.GenerateBundle(1000, Path.Combine(Settings.Default.OutputFolder, "test1_estimateAvg_{0}.txt"), true);






            //TestEnvironmentScalar test1 = new TestEnvironmentScalar(true, T, N, phi1, phi2, psi, mW, mNu, mEta, DW, DNu, DEta, x => phi1(x) + phi2(x) * mW, (x, y) => y - psi(x) - mNu);

            //test1.GenerateOne(Path.Combine(Settings.Default.OutputFolder, "test1_estimate.txt"));
            //test1.GenerateBundle(N, Path.Combine(Settings.Default.OutputFolder, "test1_estimateAvg.txt"), true);

            //int T = 25;
            //int N = 1000;
            //double mW = 0; double DW = 1;
            //double mNu = 0; double DNu = 1;
            //double mEta = 100; double DEta = 100;
            //Func<double, double> phi1 = new Func<double, double> (x => 1.0 / ( Math.Sign(x) * Math.Pow(Math.Abs(x), 1.0/3.0)));//(x => Math.Abs(x) / (x * x + 0.01) + 1.0);
            //Func<double, double> phi2 = new Func<double, double>(x => 1.0);
            //Func<double, double> psi = new Func<double, double>(x => Math.Pow(x, 3));

            //TestEnvironment test2 = new TestEnvironment(true, T, N, phi1, phi2, psi, mW, mNu, mEta, DW, DNu, DEta, x => phi1(x) + phi2(x) * mW, (x, y) => y - psi(x) - mNu);

            ////test2.GenerateOne(Path.Combine(Settings.Default.OutputFolder, "test2_estimate.txt"));
            //test2.GenerateBundle(1000, Path.Combine(Settings.Default.OutputFolder, "test2_estimateAvg.txt"), false);


            int T = 25;
            int N = 1000;
            Vector<double> mW = Vector(0, 0); Matrix<double> dW = Diag(1, 1);
            Vector<double> mNu = Vector(0, 0); Matrix<double> dNu = Diag(1, 1);
            Vector<double> mEta = Vector(100, 100); Matrix<double> dEta = Diag(100, 100);
            Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1]));
            Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Diag(1.0, 1.0);
            Func<int, Vector<double>, Vector<double>> psi = (s, x) => Vector(Math.Pow(x[0], 3) + Math.Pow(x[0], 1), Math.Pow(x[1], 3) + Math.Pow(x[1], 1));

            Normal[] NormalW = new Normal[2] { new Normal(mW[0], Math.Sqrt(dW[0, 0])), new Normal(mW[1], Math.Sqrt(dW[1, 1])) };
            Normal[] NormalNu = new Normal[2] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])), new Normal(mNu[1], Math.Sqrt(dNu[1, 1])) }; ;
            Normal[] NormalEta = new Normal[2] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])), new Normal(mEta[1], Math.Sqrt(dEta[1, 1])) }; ;

            TestEnvironmentVector test1 = new TestEnvironmentVector()
            {
                Phi1 = phi1,
                Phi2 = phi2,
                Psi = psi,
                Xi = (s, x) => phi1(s, x) + phi2(s, x) * mW,
                Zeta = (s, x, y) => y - psi(s, x) - mNu,
                W = (s) => Vector(NormalW[0].Sample(), NormalW[1].Sample()),
                Nu = (s) => Vector(NormalNu[0].Sample(), NormalNu[1].Sample()),
                DW = dW,
                DNu = dNu,
                X0 = () => Vector(NormalEta[0].Sample(), NormalEta[1].Sample()),
                X0Hat = mEta,
                DX0Hat = dEta
            };

            test1.Initialize(true, T, N);
            test1.GenerateBundle(1000, Path.Combine(Settings.Default.OutputFolder, "test1_estimateAvg_{0}.txt"), true);

        }

        static Vector<double> Vector(params double[] val)
        {
            return Vector<double>.Build.Dense(val);
        }
        static Matrix<double> Matrix(double val)
        {
            return Matrix<double>.Build.Dense(1, 1, val);
        }
        static Matrix<double> Diag(params double[] val)
        {
            return Matrix<double>.Build.DenseDiagonal(val.Length, val.Length, (i) => val[i]);
        }

    }


}
