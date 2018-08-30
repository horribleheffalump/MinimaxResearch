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
using CommandLine;

namespace CMNFTest
{
    class Program
    {
        public class Options
        {
            [Option('m', "model", Required = true, HelpText = "Model name, one of: sphere, polar, cubic, invprop-good, invprop-bad, logreg-simple, logreg-zero, logreg-uniform")]
            public string Model { get; set; }

            [Option('n', "samples", Required = false, HelpText = "Number of samples for static models (sphere and polar)")]
            public int N { get; set; }

            [Option('T', "upper-bound", Required = false, HelpText = "The upper bound of the observation interval for dynamic models")]
            public int T { get; set; }

            [Option('t', "train-count", Required = false, HelpText = "Number of trajectoris in the training set for dynamic models")]
            public int TrainCount { get; set; }

            [Option('e', "test-count", Required = false, HelpText = "Number of trajectoris in the test set for dynamic models")]
            public int TestCount { get; set; }

            [Option('U', "UKF", Required = false, Default = true, HelpText = "Do calculate Unscented Kalman Filter")]
            public bool UKF { get; set; }

            [Option('b', "bound", Required = false, HelpText = "Upper bound for the state")]
            public double Bound { get; set; }

            [Option('o', "output-folder", Required = false, HelpText = "Folder to store numeric results")]
            public string OutputFolder { get; set; }

            [Option('p', "plots-folder", Required = false, HelpText = "Folder to store plots")]
            public string PlotsFolder { get; set; }

            [Option('s', "scripts-folder", Required = false, HelpText = "Folder where the python scripts are stored")]
            public string ScriptsFolder { get; set; }

            [Option('q', "templates-folder", Required = false, HelpText = "Folder where the latex templates are stored")]
            public string TemplatesFolder { get; set; }


        }

        static void Main(string[] args)
        {
            CommandLine.Parser.Default.ParseArguments<Options>(args).WithParsed<Options>(opts => Run(opts));
        }

        static void Run(Options o)
        {
            if (string.IsNullOrWhiteSpace(o.OutputFolder))
                o.OutputFolder = Settings.Default.OutputFolder;

            if (string.IsNullOrWhiteSpace(o.PlotsFolder))
                o.PlotsFolder = Settings.Default.LatexFolder;

            if (string.IsNullOrWhiteSpace(o.ScriptsFolder))
                o.ScriptsFolder = Settings.Default.ScriptsFolder;

            if (string.IsNullOrWhiteSpace(o.TemplatesFolder))
                o.TemplatesFolder = Settings.Default.LatexFolder;

            #region sphere
            if (o.Model == "sphere")
            {
                int N = o.N;
                Vector<double> mX = Exts.Vector(30, 40, 100); Matrix<double> KX = Exts.Diag(30 * 30, 30 * 30, 30 * 30);
                Vector<double> mNu = Exts.Vector(0, 0, 0); Matrix<double> KNu = Exts.Diag(30 * 30, Math.Pow(5 * Math.PI / 180.0, 2.0), Math.Pow(5 * Math.PI / 180.0, 2.0));
                Normal[] NormalX = new Normal[3] { new Normal(mX[0], Math.Sqrt(KX[0, 0])), new Normal(mX[1], Math.Sqrt(KX[1, 1])), new Normal(mX[2], Math.Sqrt(KX[2, 2])) };
                Normal[] NormalNu = new Normal[3] { new Normal(mNu[0], Math.Sqrt(KNu[0, 0])), new Normal(mNu[1], Math.Sqrt(KNu[1, 1])), new Normal(mNu[2], Math.Sqrt(KNu[2, 2])) }; ;

                TestEnvironmentStatic testSphere = new TestEnvironmentStatic
                {
                    Phi = x => Utils.cart2sphere(x),
                    InvPhi = y => Utils.sphere2cart(y),
                    W = () => Exts.Vector(NormalX[0].Sample(), NormalX[1].Sample(), NormalX[2].Sample()),
                    Nu = () => Exts.Vector(NormalNu[0].Sample(), NormalNu[1].Sample(), NormalNu[2].Sample()),
                    MX = mX,
                    KX = KX,
                    KNu = KNu
                };

                testSphere.Initialize(N, o.OutputFolder);
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

                string fileName_alldata = Path.Combine(o.OutputFolder, "test_sphere_close_alldata.txt");
                testSphere.GenerateBundle(N, out mErr, out KErr, out KErrTh, out mErr_inv, out KErr_inv, out KErrTh_inv, out mErr_lin, out KErr_lin, out KErrTh_lin, out mErr_UT, out KErr_UT, out KErrTh_UT, fileName_alldata);


                string fileName = Path.Combine(o.OutputFolder, "test_sphere_close.txt");
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
            }
            #endregion

            #region polar 
            if (o.Model == "polar")
            {
                int N = o.N;

                Vector<double> mX = Exts.Vector(300, 400); Matrix<double> KX = Exts.Diag(30 * 30, 30 * 30);
                Vector<double> mNu = Exts.Vector(0, 0); Matrix<double> KNu = Exts.Diag(Math.Pow(5 * Math.PI / 180.0, 2.0), 30 * 30);
                Normal[] NormalX = new Normal[2] { new Normal(mX[0], Math.Sqrt(KX[0, 0])), new Normal(mX[1], Math.Sqrt(KX[1, 1])) };
                Normal[] NormalNu = new Normal[2] { new Normal(mNu[0], Math.Sqrt(KNu[0, 0])), new Normal(mNu[1], Math.Sqrt(KNu[1, 1])) }; ;

                //Console.WriteLine(mX.ToLine());

                TestEnvironmentStatic testPolar = new TestEnvironmentStatic
                {
                    Phi = x => Utils.cart2pol(x),
                    InvPhi = y => Utils.pol2cart(y),
                    W = () => Exts.Vector(NormalX[0].Sample(), NormalX[1].Sample()),
                    Nu = () => Exts.Vector(NormalNu[0].Sample(), NormalNu[1].Sample()),
                    MX = mX,
                    KX = KX,
                    KNu = KNu
                };

                testPolar.Initialize(N, o.OutputFolder);
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

                string fileName_alldata = Path.Combine(o.OutputFolder, "test_polar_far_alldata.txt");
                testPolar.GenerateBundle(N, out mErr, out KErr, out KErrTh, out mErr_inv, out KErr_inv, out KErrTh_inv, out mErr_lin, out KErr_lin, out KErrTh_lin, out mErr_UT, out KErr_UT, out KErrTh_UT, fileName_alldata);


                string fileName = Path.Combine(o.OutputFolder, "test_polar_far.txt");
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
            }
            #endregion

            #region cubic sensor
            if (o.Model == "cubic")
            {
                TestCubicSensorScalar testCubicSensor = new TestCubicSensorScalar();
                testCubicSensor.Initialize(o.T, o.TrainCount, o.UKF, o.OutputFolder);
                testCubicSensor.GenerateBundle(o.TestCount, o.OutputFolder);
                testCubicSensor.GenerateOne(o.OutputFolder);
                testCubicSensor.ProcessResults(o.OutputFolder, o.ScriptsFolder, o.PlotsFolder);
                testCubicSensor.GenerateReport(o.TemplatesFolder, o.PlotsFolder);
            }
            #endregion

            #region inverse proportion good
            if (o.Model == "invprop-good")
            {
                TestInverseProportionGoodScalar testInverseProportion = new TestInverseProportionGoodScalar();
                testInverseProportion.Initialize(o.T, o.TrainCount, o.UKF, o.OutputFolder);
                testInverseProportion.GenerateBundle(o.TestCount, o.OutputFolder);
                testInverseProportion.GenerateOne(o.OutputFolder);
                testInverseProportion.ProcessResults(o.OutputFolder, o.ScriptsFolder, o.PlotsFolder);
                testInverseProportion.GenerateReport(o.TemplatesFolder, o.PlotsFolder);
            }
            #endregion

            #region inverse proportion bad
            if (o.Model == "invprop-bad")
            {
                TestInverseProportionBadScalar testInverseProportion = new TestInverseProportionBadScalar(o.Bound);
                testInverseProportion.Initialize(o.T, o.TrainCount, o.UKF, o.OutputFolder);
                testInverseProportion.GenerateBundle(o.TestCount, o.OutputFolder);
                testInverseProportion.GenerateOne(o.OutputFolder);
                testInverseProportion.ProcessResults(o.OutputFolder, o.ScriptsFolder, o.PlotsFolder);
                testInverseProportion.GenerateReport(o.TemplatesFolder, o.PlotsFolder);
            }
            #endregion

            #region logistic regression
            if (o.Model == "logreg-simple")
            {
                TestLogisticModelScalar testLogisticModel = new TestLogisticModelScalar(o.Bound);
                testLogisticModel.Initialize(o.T, o.TrainCount, o.UKF, o.OutputFolder);
                testLogisticModel.GenerateBundle(o.TestCount, o.OutputFolder, o.UKF);
                testLogisticModel.GenerateOne(o.OutputFolder, o.UKF);
                testLogisticModel.ProcessResults(o.OutputFolder, o.ScriptsFolder, o.PlotsFolder);
                testLogisticModel.GenerateReport(o.TemplatesFolder, o.PlotsFolder);
            }
            #endregion

            #region logistic regression zero
            if (o.Model == "logreg-zero")
            {
                TestLogisticModelZeroScalar testLogisticZeroModel = new TestLogisticModelZeroScalar(o.Bound);
                testLogisticZeroModel.Initialize(o.T, o.TrainCount, o.UKF, o.OutputFolder);
                testLogisticZeroModel.GenerateBundle(o.TestCount, o.OutputFolder, o.UKF);
                testLogisticZeroModel.GenerateOne(o.OutputFolder, o.UKF);
                testLogisticZeroModel.ProcessResults(o.OutputFolder, o.ScriptsFolder, o.PlotsFolder);
                testLogisticZeroModel.GenerateReport(o.TemplatesFolder, o.PlotsFolder);
            }
            #endregion

            #region logistic regression uniform noise
            if (o.Model == "logreg-uniform")
            {
                TestLogisticModelUniformNoiseScalar testLogisticUniformNoiseModel = new TestLogisticModelUniformNoiseScalar();
                testLogisticUniformNoiseModel.Initialize(o.T, o.TrainCount, o.UKF, o.OutputFolder);
                testLogisticUniformNoiseModel.GenerateBundle(o.TestCount, o.OutputFolder, o.UKF);
                testLogisticUniformNoiseModel.GenerateOne(o.OutputFolder, o.UKF);
                testLogisticUniformNoiseModel.ProcessResults(o.OutputFolder, o.ScriptsFolder, o.PlotsFolder);
                testLogisticUniformNoiseModel.GenerateReport(o.TemplatesFolder, o.PlotsFolder);
                #endregion
            }

        }
    }
}
