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
using MathNetExtensions;

namespace CMNFTest
{
    class Program
    {
        static void Main(string[] args)
        {
            #region sphere
            //int N = 10000;
            //Vector<double> mX = Exts.Vector(30, 40, 100); Matrix<double> KX = Exts.Diag(30 * 30, 30 * 30, 30 * 30);
            //Vector<double> mNu = Exts.Vector(0, 0, 0); Matrix<double> KNu = Exts.Diag(30 * 30, Math.Pow(5 * Math.PI / 180.0, 2.0), Math.Pow(5 * Math.PI / 180.0, 2.0));
            //Normal[] NormalX = new Normal[3] { new Normal(mX[0], Math.Sqrt(KX[0, 0])), new Normal(mX[1], Math.Sqrt(KX[1, 1])), new Normal(mX[2], Math.Sqrt(KX[2, 2])) };
            //Normal[] NormalNu = new Normal[3] { new Normal(mNu[0], Math.Sqrt(KNu[0, 0])), new Normal(mNu[1], Math.Sqrt(KNu[1, 1])), new Normal(mNu[2], Math.Sqrt(KNu[2, 2])) }; ;

            //TestEnvironmentStatic testSphere = new TestEnvironmentStatic
            //{
            //    Phi = x => Utils.cart2sphere(x),
            //    InvPhi = y => Utils.sphere2cart(y),
            //    W = () => Exts.Vector(NormalX[0].Sample(), NormalX[1].Sample(), NormalX[2].Sample()),
            //    Nu = () => Exts.Vector(NormalNu[0].Sample(), NormalNu[1].Sample(), NormalNu[2].Sample()),
            //    MX = mX,
            //    KX = KX,
            //    KNu = KNu
            //};

            //testSphere.Initialize(N, Settings.Default.OutputFolder);
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

            //string fileName_alldata = Path.Combine(Settings.Default.OutputFolder, "test_sphere_close_alldata.txt");
            //testSphere.GenerateBundle(N, out mErr, out KErr, out KErrTh, out mErr_inv, out KErr_inv, out KErrTh_inv, out mErr_lin, out KErr_lin, out KErrTh_lin, out mErr_UT, out KErr_UT, out KErrTh_UT, fileName_alldata);


            //string fileName = Path.Combine(Settings.Default.OutputFolder, "test_sphere_close.txt");
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
            #endregion

            #region polar 
            //Vector<double> mX = Exts.Vector(300, 400); Matrix<double> KX = Exts.Diag(30 * 30, 30 * 30);
            //Vector<double> mNu = Exts.Vector(0, 0); Matrix<double> KNu = Exts.Diag(Math.Pow(5 * Math.PI / 180.0, 2.0), 30 * 30);
            //Normal[] NormalX = new Normal[2] { new Normal(mX[0], Math.Sqrt(KX[0, 0])), new Normal(mX[1], Math.Sqrt(KX[1, 1])) };
            //Normal[] NormalNu = new Normal[2] { new Normal(mNu[0], Math.Sqrt(KNu[0, 0])), new Normal(mNu[1], Math.Sqrt(KNu[1, 1])) }; ;

            ////Console.WriteLine(mX.ToLine());

            //TestEnvironmentStatic testPolar = new TestEnvironmentStatic
            //{
            //    Phi = x => Utils.cart2pol(x),
            //    InvPhi = y => Utils.pol2cart(y),
            //    W = () => Exts.Vector(NormalX[0].Sample(), NormalX[1].Sample()),
            //    Nu = () => Exts.Vector(NormalNu[0].Sample(), NormalNu[1].Sample()),
            //    MX = mX,
            //    KX = KX,
            //    KNu = KNu
            //};

            //int N = 10000;
            //testPolar.Initialize(N, Settings.Default.OutputFolder);
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

            //string fileName_alldata = Path.Combine(Settings.Default.OutputFolder, "test_polar_far_alldata.txt");
            //testPolar.GenerateBundle(N, out mErr, out KErr, out KErrTh, out mErr_inv, out KErr_inv, out KErrTh_inv, out mErr_lin, out KErr_lin, out KErrTh_lin, out mErr_UT, out KErr_UT, out KErrTh_UT, fileName_alldata);


            //string fileName = Path.Combine(Settings.Default.OutputFolder, "test_polar_far.txt");
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
            #endregion

            #region cubic sensor
            //TestCubicSensorScalar testCubicSensor = new TestCubicSensorScalar();
            //testCubicSensor.Initialize(50, 10000, true, Settings.Default.OutputFolder);
            //testCubicSensor.GenerateBundle(10000, Settings.Default.OutputFolder);
            //testCubicSensor.GenerateOne(Settings.Default.OutputFolder);
            //testCubicSensor.ProcessResults(Settings.Default.OutputFolder, Settings.Default.ScriptsFolder, Settings.Default.LatexFolder);
            //testCubicSensor.GenerateReport(Settings.Default.LatexFolder);
            #endregion

            #region inverse proportion good
            //TestInverseProportionGoodScalar testInverseProportion = new TestInverseProportionGoodScalar();
            //testInverseProportion.Initialize(50, 10000, true, Settings.Default.OutputFolder);
            //testInverseProportion.GenerateBundle(10000, Settings.Default.OutputFolder);
            //testInverseProportion.GenerateOne(Settings.Default.OutputFolder);
            //testInverseProportion.ProcessResults(Settings.Default.OutputFolder, Settings.Default.ScriptsFolder, Settings.Default.LatexFolder);
            //testInverseProportion.GenerateReport(Settings.Default.LatexFolder);
            #endregion

            #region inverse proportion bad
            //TestInverseProportionBadScalar testInverseProportion = new TestInverseProportionBadScalar();
            //testInverseProportion.Initialize(50, 10000, true, Settings.Default.OutputFolder);
            //testInverseProportion.GenerateBundle(1000000, Settings.Default.OutputFolder);
            //testInverseProportion.GenerateOne(Settings.Default.OutputFolder);
            //testInverseProportion.ProcessResults(Settings.Default.OutputFolder, Settings.Default.ScriptsFolder, Settings.Default.LatexFolder);
            //testInverseProportion.GenerateReport(Settings.Default.LatexFolder); 
            #endregion

            #region logistic regression
            //TestLogisticModelScalar testLogisticModel = new TestLogisticModelScalar();
            //testLogisticModel.Initialize(50, 1000, true, Settings.Default.OutputFolder);
            //testLogisticModel.GenerateBundle(10000, Settings.Default.OutputFolder, true);
            //testLogisticModel.GenerateOne(Settings.Default.OutputFolder, true);
            //testLogisticModel.ProcessResults(Settings.Default.OutputFolder, Settings.Default.ScriptsFolder, Settings.Default.LatexFolder);
            //testLogisticModel.GenerateReport(Settings.Default.LatexFolder);
            #endregion

            #region logistic regression zero
            //TestLogisticModelZeroScalar testLogisticZeroModel = new TestLogisticModelZeroScalar();
            //testLogisticZeroModel.Initialize(50, 1000, true, Settings.Default.OutputFolder);
            //testLogisticZeroModel.GenerateBundle(10000, Settings.Default.OutputFolder, true);
            //testLogisticZeroModel.GenerateOne(Settings.Default.OutputFolder, true);
            //testLogisticZeroModel.ProcessResults(Settings.Default.OutputFolder, Settings.Default.ScriptsFolder, Settings.Default.LatexFolder);
            //testLogisticZeroModel.GenerateReport(Settings.Default.LatexFolder);
            #endregion

            #region logistic regression uniform noise
            TestLogisticModelUniformNoiseScalar testLogisticUniformNoiseModel = new TestLogisticModelUniformNoiseScalar();
            //testLogisticUniformNoiseModel.Initialize(50, 1000, true, Settings.Default.OutputFolder);
            //testLogisticUniformNoiseModel.GenerateBundle(10000, Settings.Default.OutputFolder, true);
            //testLogisticUniformNoiseModel.GenerateOne(Settings.Default.OutputFolder, true);
            testLogisticUniformNoiseModel.ProcessResults(Settings.Default.OutputFolder, Settings.Default.ScriptsFolder, Settings.Default.LatexFolder);
            testLogisticUniformNoiseModel.GenerateReport(Settings.Default.LatexFolder);
            #endregion

        }
    }
}
