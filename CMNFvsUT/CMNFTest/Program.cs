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
using PythonInteract;
using System.Linq.Expressions;

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


            //polar!

            //sphere
            //int N = 100000;
            //Vector<double> mX = Vector(30, 40, 100); Matrix<double> KX = Diag(30 * 30, 30 * 30, 30 * 30);
            //Vector<double> mNu = Vector(0, 0, 0); Matrix<double> KNu = Diag(30 * 30, Math.Pow(5 * Math.PI / 180.0 / 60.0, 2.0), Math.Pow(5 * Math.PI / 180.0 / 60.0, 2.0));
            //Normal[] NormalX = new Normal[3] { new Normal(mX[0], Math.Sqrt(KX[0, 0])), new Normal(mX[1], Math.Sqrt(KX[1, 1])), new Normal(mX[2], Math.Sqrt(KX[2, 2])) };
            //Normal[] NormalNu = new Normal[3] { new Normal(mNu[0], Math.Sqrt(KNu[0, 0])), new Normal(mNu[1], Math.Sqrt(KNu[1, 1])), new Normal(mNu[2], Math.Sqrt(KNu[2, 2])) }; ;

            //TestEnvironmentStatic testPolar = new TestEnvironmentStatic
            //{
            //    Phi = x => Extensions.cart2sphere(x),
            //    InvPhi = y => Extensions.sphere2cart(y),
            //    W = () => Vector(NormalX[0].Sample(), NormalX[1].Sample(), NormalX[2].Sample()),
            //    Nu = () => Vector(NormalNu[0].Sample(), NormalNu[1].Sample(), NormalNu[2].Sample()),
            //    MX = mX,
            //    KX = KX,
            //    KNu = KNu,
            //    utOptimizationType = UTOptimizationType.ImplicitAlphaBetaKappa
            //};

            //testPolar.Initialize(N, 100, 100);
            //Vector<double> mErr;
            //Matrix<double> KErr;
            //Matrix<double> KErrTh;
            //Vector<double> mErr_inv;
            //Matrix<double> KErr_inv;
            //Matrix<double> KErrTh_inv;
            //Vector<double> mErr_lin;
            //Matrix<double> KErr_lin;
            //Matrix<double> KErrTh_lin;
            //Vector<double> mErr_UT;
            //Matrix<double> KErr_UT;
            //Matrix<double> KErrTh_UT;

            //string fileName_alldata = Path.Combine(Settings.Default.OutputFolder, "test_polar_alldata.txt");
            //testPolar.GenerateBundle(N, out mErr, out KErr, out KErrTh, out mErr_inv, out KErr_inv, out KErrTh_inv, out mErr_lin, out KErr_lin, out KErrTh_lin, out mErr_UT, out KErr_UT, out KErrTh_UT, fileName_alldata);


            //string fileName = Path.Combine(Settings.Default.OutputFolder, "test_polar.txt");
            //using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
            //{
            //    //outputfile.WriteLine($"P = {P}");
            //    outputfile.WriteLine($"mErr = {mErr}");
            //    outputfile.WriteLine($"KErr = {KErr}");
            //    outputfile.WriteLine($"KErrTh = {KErrTh}");

            //    //outputfile.WriteLine($"P_inv = {P_inv}");
            //    outputfile.WriteLine($"mErr_inv = {mErr_inv}");
            //    outputfile.WriteLine($"KErr_inv = {KErr_inv}");
            //    outputfile.WriteLine($"KErrTh_inv = {KErrTh_inv}");

            //    //outputfile.WriteLine($"P_lin = {P_lin}");
            //    outputfile.WriteLine($"mErr_lin = {mErr_lin}");
            //    outputfile.WriteLine($"KErr_lin = {KErr_lin}");
            //    outputfile.WriteLine($"KErrTh_lin = {KErrTh_lin}");

            //    //outputfile.WriteLine($"P_UT = {P_UT}");
            //    outputfile.WriteLine($"mErr_UT = {mErr_UT}");
            //    outputfile.WriteLine($"KErr_UT = {KErr_UT}");
            //    outputfile.WriteLine($"KErrTh_UT = {KErrTh_UT}");

            //    outputfile.Close();
            //}


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

            TestCubicSensor testCubicSensor = new TestCubicSensor();
            testCubicSensor.Initialize(100, 100, true, Settings.Default.OutputFolder);
            testCubicSensor.GenerateBundle(100, Settings.Default.OutputFolder);
            testCubicSensor.GenerateOne(Settings.Default.OutputFolder);
            testCubicSensor.ProcessResults(Settings.Default.OutputFolder, Settings.Default.ScriptsFolder, Settings.Default.LatexFolder);
            testCubicSensor.GenerateReport(Settings.Default.LatexFolder);
         }
    }
}
