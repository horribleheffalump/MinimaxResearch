using CMNFvsUTFTest.Properties;
using CommandLine;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using NonlinearSystem;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TestEnvironments;

namespace CMNFvsUTFTest
{
    class Program
    {
        public class Options
        {
            [Option('m', "model", Required = true, HelpText = "Model name, one of: sphere, polar, polartwo, cubic, invprop-good, invprop-bad, logreg-simple, logreg-zero, logreg-uniform, samplereg, switchingobs")]
            public string Model { get; set; }

            [Option('n', "samples", Required = false, HelpText = "Number of samples for static models (sphere and polar)")]
            public int N { get; set; }

            [Option('T', "upper-bound", Required = false, HelpText = "The upper bound of the observation interval for dynamic models")]
            public int T { get; set; }

            [Option('t', "train-count", Required = false, HelpText = "Number of trajectoris in the training set for dynamic models")]
            public int TrainCount { get; set; }

            [Option('e', "test-count", Required = false, HelpText = "Number of trajectoris in the test set for dynamic models")]
            public int TestCount { get; set; }

            [Option('E', "bundle-count", Required = false, Default = 1, HelpText = "Number of bundles of trajectoris in the test set for dynamic models")]
            public int BundleCount { get; set; }

            [Option('U', "UKF", Required = false, Default = false, HelpText = "Do calculate Unscented Kalman Filter")]
            public bool UKF { get; set; }

            [Option('S', "UKF-stepwise", Required = false, Default = false, HelpText = "Do calculate Unscented Kalman Filter with stepwise parameter estimation")]
            public bool UKFStepwise { get; set; }

            [Option('R', "UKF-randomshoot", Required = false, Default = false, HelpText = "Do optimize Unscented Kalman Filter parameters with random shoot method")]
            public bool UKFRandomShoot { get; set; }

            [Option('C', "CMNF", Required = false, Default = false, HelpText = "Do calculate conditionally minimax nonlinear filter")]
            public bool CMNF { get; set; }

            [Option('M', "MCMNF", Required = false, Default = false, HelpText = "Do calculate modified conditionally minimax nonlinear filter")]
            public bool MCMNF { get; set; }

            [Option('l', "MCMNF-train-count", Required = false, Default = 0, HelpText = "Number of test set points for MCMNF wtepwise parameters fitting")]
            public int MCMNFTrainCount { get; set; }

            [Option('b', "bound", Required = false, HelpText = "Upper bound for the state")]
            public double Bound { get; set; }

            [Option('W', "DW", Required = false, Default = 1.0, HelpText = "State noise variation")]
            public double DW { get; set; }

            [Option('N', "DNu", Required = false, Default = 1.0, HelpText = "Observations noise variation")]
            public double DNu { get; set; }

            [Option('o', "output-folder", Required = false, HelpText = "Folder to store numeric results")]
            public string OutputFolder { get; set; }

            [Option('p', "plots-folder", Required = false, HelpText = "Folder to store plots")]
            public string PlotsFolder { get; set; }

            [Option('s', "scripts-folder", Required = false, HelpText = "Folder where the python scripts are stored")]
            public string ScriptsFolder { get; set; }

            [Option('q', "templates-folder", Required = false, HelpText = "Folder where the latex templates are stored")]
            public string TemplatesFolder { get; set; }

            [Option('B', "bulk", Required = false, Default = false, HelpText = "Bulk output of the state process trajectories bundle")]
            public bool Bulk { get; set; }

            [Option('G', "generate-samples", Required = false, Default = 1, HelpText = "Number of single trajectories to generate")]
            public int SamplesCount { get; set; }

            [Option('P', "parallel", Required = false, Default = false, HelpText = "Parallel bundles calculation (on to save time, off to save memory in case of multiple large bundles)")]
            public bool Parallel { get; set; }



        }

        static void Main(string[] args)
        {
            CommandLine.Parser.Default.ParseArguments<Options>(args).WithParsed<Options>(opts => Run(opts, args));
        }

        static void Run(Options o, string[] args)
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

                //Vector<double> mX = Exts.Vector(300, 400); Matrix<double> KX = Exts.Diag(30 * 30, 30 * 30);
                //Vector<double> mX = Exts.Vector(30000, 40000); Matrix<double> KX = Exts.Diag(100 * 100, 100 * 100);
                Vector<double> mX = Exts.Vector(30000, 40000); Matrix<double> KX = Exts.Diag(4500 * 4500, 4500 * 4500);
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

            #region polartwo 
            if (o.Model == "polartwo")
            {
                int N = o.N;
                Vector<double> secondpoint = Exts.Vector(10000, 0);
                Vector<double> mX = Exts.Vector(30000, 40000); Matrix<double> KX = Exts.Diag(4500 * 4500, 4500 * 4500);
                Vector<double> mNu = Exts.Vector(0, 0, 0, 0);
                Matrix<double> KNu = Exts.Diag(Math.Pow(5 * Math.PI / 180.0, 2.0),
                                                30 * 30,
                                                Math.Pow(5 * Math.PI / 180.0, 2.0),
                                                30 * 30);
                Normal[] NormalX = new Normal[2] { new Normal(mX[0], Math.Sqrt(KX[0, 0])), new Normal(mX[1], Math.Sqrt(KX[1, 1])) };
                Normal[] NormalNu = new Normal[4] { new Normal(mNu[0], Math.Sqrt(KNu[0, 0])),
                                                    new Normal(mNu[1], Math.Sqrt(KNu[1, 1])),
                                                    new Normal(mNu[2], Math.Sqrt(KNu[2, 2])),
                                                    new Normal(mNu[3], Math.Sqrt(KNu[3, 3])) };

                //Console.WriteLine(mX.ToLine());

                TestEnvironmentStatic testPolar = new TestEnvironmentStatic
                {
                    Phi = x => Exts.Stack(Utils.cart2pol(x), Utils.cart2pol(x - secondpoint)),
                    InvPhi = y => Exts.Stack(Utils.pol2cart(Exts.Vector(y[0], y[1])), Utils.pol2cart(Exts.Vector(y[2], y[3])) + secondpoint),
                    W = () => Exts.Vector(NormalX[0].Sample(), NormalX[1].Sample()),
                    Nu = () => Exts.Vector(NormalNu[0].Sample(), NormalNu[1].Sample(), NormalNu[2].Sample(), NormalNu[3].Sample()),
                    //Nu = () => Exts.Vector(0,0,0,0),
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

            TestEnvironmentVector testEnv = new TestEnvironmentVector();

            if (o.Model == "cubic")
            {
                testEnv = new TestCubicSensorScalar(o.DW, o.DNu);
            }
            if (o.Model == "invprop-good")
            {
                testEnv = new TestInverseProportionGoodScalar(o.Bound, o.DW, o.DNu);
            }
            if (o.Model == "invprop-bad")
            {
                testEnv = new TestInverseProportionBadScalar(o.Bound, o.DW, o.DNu);
            }
            if (o.Model == "logreg-simple")
            {
                testEnv = new TestLogisticModelScalar(o.Bound, o.DW, o.DNu);
            }
            if (o.Model == "logreg-zero")
            {
                testEnv = new TestLogisticModelZeroScalar(o.Bound, o.DW, o.DNu);
            }
            if (o.Model == "logreg-uniform")
            {
                testEnv = new TestLogisticModelUniformNoiseScalar();
            }
            if (o.Model == "samplereg")
            {
                testEnv = new TestSampledRegression(o.DNu);
            }
            if (o.Model == "switchingobs")
            {
                testEnv = new TestSwitchingObservations(o.DNu);
            }

            List<FilterType> filters = new List<FilterType>();
            if (o.CMNF) filters.Add(FilterType.CMNF);
            if (o.MCMNF) filters.Add(FilterType.MCMNF);
            if (o.UKF)
            {
                if (o.UKFStepwise)
                {
                    if (o.UKFRandomShoot) filters.Add(FilterType.UKFStepwiseRandomShoot);
                    else filters.Add(FilterType.UKFStepwise);
                }
                else
                {
                    if (o.UKFRandomShoot) filters.Add(FilterType.UKFIntegralRandomShoot);
                    else filters.Add(FilterType.UKFIntegral);
                }
            }

            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(Path.Combine(o.OutputFolder, "parameters.txt"), true))
            {
                outputfile.WriteLine($"{DateTime.Now}\t{string.Join(" ", args)}");
                outputfile.Close();
            }

            if (o.Bulk)
                testEnv.GenerateBundleSamples(o.T, o.TrainCount, o.OutputFolder);
            else
            {
                testEnv.Initialize(o.T, o.TrainCount, o.MCMNFTrainCount, o.OutputFolder, filters);
                if (o.BundleCount > 1)
                    testEnv.GenerateBundles(o.BundleCount, o.TestCount, o.OutputFolder, o.Parallel);
                else
                    testEnv.GenerateBundle(o.TestCount, o.OutputFolder);
                if (o.SamplesCount == 1)
                    testEnv.GenerateOne(o.OutputFolder);
                else
                {
                    for (int i = 0; i < o.SamplesCount; i++)
                    {
                        testEnv.GenerateOne(o.OutputFolder, i);
                    }
                }
                testEnv.ProcessResults(o.OutputFolder, o.ScriptsFolder, o.PlotsFolder);
                //testEnv.GenerateReport(o.TemplatesFolder, o.PlotsFolder);
            }
        }
    }
}
