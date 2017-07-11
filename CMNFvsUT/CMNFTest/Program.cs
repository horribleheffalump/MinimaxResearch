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
            NumberFormatInfo provider;
            provider = new NumberFormatInfo()
            {
                NumberDecimalSeparator = "."
            };

            //Vector<double>[] coords = new Vector<double>[27];
            //coords[0] = Vector(0, 0, 0);
            //coords[1] = Vector(0, 0, 1);
            //coords[2] = Vector(0, 0, -1);
            //coords[3] = Vector(0, 1, 0);
            //coords[4] = Vector(0, 1, 1);
            //coords[5] = Vector(0, 1, -1);
            //coords[6] = Vector(0, -1, 0);
            //coords[7] = Vector(0, -1, 1);
            //coords[8] = Vector(0, -1, -1);

            //coords[9] = Vector(1, 0, 0);
            //coords[10] = Vector(1, 0, 1);
            //coords[11] = Vector(1, 0, -1);
            //coords[12] = Vector(1, 1, 0);
            //coords[13] = Vector(1, 1, 1);
            //coords[14] = Vector(1, 1, -1);
            //coords[15] = Vector(1, -1, 0);
            //coords[16] = Vector(1, -1, 1);
            //coords[17] = Vector(1, -1, -1);

            //coords[18] = Vector(-1, 0, 0);
            //coords[19] = Vector(-1, 0, 1);
            //coords[20] = Vector(-1, 0, -1);
            //coords[21] = Vector(-1, 1, 0);
            //coords[22] = Vector(-1, 1, 1);
            //coords[23] = Vector(-1, 1, -1);
            //coords[24] = Vector(-1, -1, 0);
            //coords[25] = Vector(-1, -1, 1);
            //coords[26] = Vector(-1, -1, -1);

            //for (int i = 1; i < 27; i++)
            //{
            //    Vector<double> sp = Extensions.cart2sphere(coords[i]);
            //    Vector<double> inv = Extensions.sphere2cart(sp);
            //    Console.WriteLine($"({coords[i][0]}, {coords[i][1]}, {coords[i][2]}) => ({sp[0]}, {sp[1]}, {sp[2]}) => ({inv[0]}, {inv[1]}, {inv[2]})");

            //}
            //Console.WriteLine();

            //dimension 1
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
            //dimension 1




            //dimension 1 scalar
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
            //dimension 1 scalar

            //dimension 2
            //int T = 25;
            //int N = 1000;
            //Vector<double> mW = Vector(0, 0); Matrix<double> dW = Diag(1, 1);
            //Vector<double> mNu = Vector(0, 0); Matrix<double> dNu = Diag(1, 1);
            //Vector<double> mEta = Vector(100, 100); Matrix<double> dEta = Diag(100, 100);
            //Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Vector(x[0] / (1 + x[0] * x[0]), x[1] / (1 + x[1] * x[1]));
            //Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Diag(1.0, 1.0);
            //Func<int, Vector<double>, Vector<double>> psi = (s, x) => Vector(Math.Pow(x[0], 3) + Math.Pow(x[0], 1), Math.Pow(x[1], 3) + Math.Pow(x[1], 1));

            //Normal[] NormalW = new Normal[2] { new Normal(mW[0], Math.Sqrt(dW[0, 0])), new Normal(mW[1], Math.Sqrt(dW[1, 1])) };
            //Normal[] NormalNu = new Normal[2] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])), new Normal(mNu[1], Math.Sqrt(dNu[1, 1])) }; ;
            //Normal[] NormalEta = new Normal[2] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])), new Normal(mEta[1], Math.Sqrt(dEta[1, 1])) }; ;

            //TestEnvironmentVector test1 = new TestEnvironmentVector()
            //{
            //    Phi1 = phi1,
            //    Phi2 = phi2,
            //    Psi = psi,
            //    Xi = (s, x) => phi1(s, x) + phi2(s, x) * mW,
            //    Zeta = (s, x, y) => y - psi(s, x) - mNu,
            //    W = (s) => Vector(NormalW[0].Sample(), NormalW[1].Sample()),
            //    Nu = (s) => Vector(NormalNu[0].Sample(), NormalNu[1].Sample()),
            //    DW = dW,
            //    DNu = dNu,
            //    X0 = () => Vector(NormalEta[0].Sample(), NormalEta[1].Sample()),
            //    X0Hat = mEta,
            //    DX0Hat = dEta
            //};

            //test1.Initialize(true, T, N);
            //test1.GenerateBundle(1000, Path.Combine(Settings.Default.OutputFolder, "test1_estimateAvg_{0}.txt"), true);
            //dimension 2


            //polar???
            //int T = 1;
            //int N = 1000;

            //Vector<double> mW = Vector(30, 40); Matrix<double> dW = Diag(30*30, 30*30);
            //Vector<double> mNu = Vector(0, 0); Matrix<double> dNu = Diag(Math.Pow(5*Math.PI/ 180.0, 2.0), 30*30);
            //Vector<double> mEta = Vector(30, 40); Matrix<double> dEta = Diag(30 * 30, 30 * 30);
            //Func<int, Vector<double>, Vector<double>> phi1 = (s, x) => Vector(0.0, 0.0);
            //Func<int, Vector<double>, Matrix<double>> phi2 = (s, x) => Diag(1.0, 1.0);
            //Func<int, Vector<double>, Vector<double>> psi = (s, x) => Extentions.cart2pol(x);

            //Normal[] NormalW = new Normal[2] { new Normal(mW[0], Math.Sqrt(dW[0, 0])), new Normal(mW[1], Math.Sqrt(dW[1, 1])) };
            //Normal[] NormalNu = new Normal[2] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])), new Normal(mNu[1], Math.Sqrt(dNu[1, 1])) }; ;
            //Normal[] NormalEta = new Normal[2] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0])), new Normal(mEta[1], Math.Sqrt(dEta[1, 1])) }; ;

            //TestEnvironmentVector test1 = new TestEnvironmentVector()
            //{
            //    Phi1 = phi1,
            //    Phi2 = phi2,
            //    Psi = psi,
            //    Xi = (s, x) => phi1(s, x) + phi2(s, x) * mW,
            //    Zeta = (s, x, y) => y - psi(s, x) - mNu,
            //    W = (s) => Vector(NormalW[0].Sample(), NormalW[1].Sample()),
            //    Nu = (s) => Vector(NormalNu[0].Sample(), NormalNu[1].Sample()),
            //    DW = dW,
            //    DNu = dNu,
            //    X0 = () => Vector(NormalEta[0].Sample(), NormalEta[1].Sample()),
            //    X0Hat = mEta,
            //    DX0Hat = dEta
            //};
            //test1.Initialize(true, T, N);
            //test1.GenerateBundle(1000, Path.Combine(Settings.Default.OutputFolder, "test1_estimateAvg_{0}.txt"), true);
            //polar???

            //polar!

            //int N = 10000;
            //Vector<double> mX = Vector(30, 40); Matrix<double> KX = Diag(30 * 30, 30 * 30);
            //Vector<double> mNu = Vector(0, 0); Matrix<double> KNu = Diag(Math.Pow(5 * Math.PI / 180.0, 2.0), 30 * 30);
            //Normal[] NormalX = new Normal[2] { new Normal(mX[0], Math.Sqrt(KX[0, 0])), new Normal(mX[1], Math.Sqrt(KX[1, 1])) };
            //Normal[] NormalNu = new Normal[2] { new Normal(mNu[0], Math.Sqrt(KNu[0, 0])), new Normal(mNu[1], Math.Sqrt(KNu[1, 1])) }; ;

            ////Console.WriteLine(mX.ToLine());

            //TestEnvironmentStatic testPolar = new TestEnvironmentStatic
            //{
            //    Phi = x => Extensions.cart2pol(x),
            //    InvPhi = y => Extensions.pol2cart(y),
            //    W = () => Vector(NormalX[0].Sample(), NormalX[1].Sample()),
            //    Nu = () => Vector(NormalNu[0].Sample(), NormalNu[1].Sample()),
            //    MX = mX,
            //    KX = KX,
            //    KNu = KNu,
            //    utOptimizationType = UTOptimizationType.ImplicitAlphaBetaKappa
            //};


            //int N = 100000;
            //Vector<double> mX = Vector(30); Matrix<double> KX = Diag(30 * 30);
            //Vector<double> mNu = Vector(0); Matrix<double> KNu = Diag(30 * 30);
            //Normal NormalX = new Normal(mX[0], Math.Sqrt(KX[0, 0]));
            //Normal NormalNu = new Normal(mNu[0], Math.Sqrt(KNu[0, 0]));

            //TestEnvironmentStatic testPolar = new TestEnvironmentStatic
            //{
            //    Phi = x => Vector(Math.Pow(x[0], 3.0)),
            //    InvPhi = y => Vector(Math.Pow(y[0], -3.0)),
            //    W = () => Vector(NormalX.Sample()),
            //    Nu = () => Vector(NormalNu.Sample()),
            //    MX = mX,
            //    KX = KX,
            //    KNu = KNu
            //};

            //int N = 100000;
            //Vector<double> mX = Vector(10); Matrix<double> KX = Diag(1 * 1);
            //Vector<double> mNu = Vector(0); Matrix<double> KNu = Diag(1 * 1);
            //Normal NormalX = new Normal(mX[0], Math.Sqrt(KX[0, 0]));
            //Normal NormalNu = new Normal(mNu[0], Math.Sqrt(KNu[0, 0]));

            //TestEnvironmentStatic testPolar = new TestEnvironmentStatic
            //{
            //    Phi = x => Vector(Math.Pow(x[0], 3.0)),
            //    InvPhi = y => Vector(Math.Pow(Math.Abs(y[0]), 1.0/3.0)* Math.Sign(y[0])),
            //    W = () => Vector(NormalX.Sample()),
            //    Nu = () => Vector(NormalNu.Sample()),
            //    MX = mX,
            //    KX = KX,
            //    KNu = KNu
            //};

            //sphere
            int N = 1000;
            Vector<double> mX = Vector(30, 40, 100); Matrix<double> KX = Diag(30 * 30, 30 * 30, 30 * 30);
            Vector<double> mNu = Vector(0, 0, 0); Matrix<double> KNu = Diag(30 * 30, Math.Pow(5 * Math.PI / 180.0 / 60.0, 2.0), Math.Pow(5 * Math.PI / 180.0 / 60.0, 2.0));
            Normal[] NormalX = new Normal[3] { new Normal(mX[0], Math.Sqrt(KX[0, 0])), new Normal(mX[1], Math.Sqrt(KX[1, 1])), new Normal(mX[2], Math.Sqrt(KX[2, 2])) };
            Normal[] NormalNu = new Normal[3] { new Normal(mNu[0], Math.Sqrt(KNu[0, 0])), new Normal(mNu[1], Math.Sqrt(KNu[1, 1])), new Normal(mNu[2], Math.Sqrt(KNu[2, 2])) }; ;

            TestEnvironmentStatic testPolar = new TestEnvironmentStatic
            {
                Phi = x => Extensions.cart2sphere(x),
                InvPhi = y => Extensions.sphere2cart(y),
                W = () => Vector(NormalX[0].Sample(), NormalX[1].Sample(), NormalX[2].Sample()),
                Nu = () => Vector(NormalNu[0].Sample(), NormalNu[1].Sample(), NormalNu[2].Sample()),
                MX = mX,
                KX = KX,
                KNu = KNu,
                utOptimizationType = UTOptimizationType.ImplicitAlphaBetaKappa
            };

            testPolar.Initialize(N, 100, 100);
            Vector<double> mErr;
            Matrix<double> KErr;
            Matrix<double> KErrTh;
            Vector<double> mErr_inv;
            Matrix<double> KErr_inv;
            Matrix<double> KErrTh_inv;
            Vector<double> mErr_lin;
            Matrix<double> KErr_lin;
            Matrix<double> KErrTh_lin;
            Vector<double> mErr_UT;
            Matrix<double> KErr_UT;
            Matrix<double> KErrTh_UT;

            string fileName_alldata = Path.Combine(Settings.Default.OutputFolder, "test_polar_alldata.txt");
            testPolar.GenerateBundle(N, out mErr, out KErr, out KErrTh, out mErr_inv, out KErr_inv, out KErrTh_inv, out mErr_lin, out KErr_lin, out KErrTh_lin, out mErr_UT, out KErr_UT, out KErrTh_UT, fileName_alldata);


            string fileName = Path.Combine(Settings.Default.OutputFolder, "test_polar.txt");
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
            {
                //outputfile.WriteLine($"P = {P}");
                outputfile.WriteLine($"mErr = {mErr}");
                outputfile.WriteLine($"KErr = {KErr}");
                outputfile.WriteLine($"KErrTh = {KErrTh}");

                //outputfile.WriteLine($"P_inv = {P_inv}");
                outputfile.WriteLine($"mErr_inv = {mErr_inv}");
                outputfile.WriteLine($"KErr_inv = {KErr_inv}");
                outputfile.WriteLine($"KErrTh_inv = {KErrTh_inv}");

                //outputfile.WriteLine($"P_lin = {P_lin}");
                outputfile.WriteLine($"mErr_lin = {mErr_lin}");
                outputfile.WriteLine($"KErr_lin = {KErr_lin}");
                outputfile.WriteLine($"KErrTh_lin = {KErrTh_lin}");

                //outputfile.WriteLine($"P_UT = {P_UT}");
                outputfile.WriteLine($"mErr_UT = {mErr_UT}");
                outputfile.WriteLine($"KErr_UT = {KErr_UT}");
                outputfile.WriteLine($"KErrTh_UT = {KErrTh_UT}");

                outputfile.Close();
            }


            //string fileName = Path.Combine(Settings.Default.OutputFolder, "test_polar.txt");
            //using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
            //{
            //    NumberFormatInfo provider;
            //    provider = new NumberFormatInfo();
            //    provider.NumberDecimalSeparator = ".";
            //    for (int i = 0; i < X.Count(); i++)
            //    {
            //        outputfile.WriteLine(string.Format(provider, "{0} {1} {2} {3}", X[i][0], X[i][1], Xinv[i][0], Xinv[i][1]));
            //    }
            //    outputfile.Close();
            //}
            //polar!


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
