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

namespace TargetTrackingTest
{
    class Program
    {
        public class Options
        {
            [Option('T', "upper-bound", Required = false, HelpText = "The upper bound of the observation interval for dynamic models")]
            public int T { get; set; }

            [Option("h-state", Required = false, Default = 0.01, HelpText = "Discretization step for stochastic differential equation of the state dynamics")]
            public double h_state { get; set; }

            [Option("h-obs", Required = false, Default = 0.01, HelpText = "Intervals between the observations")]
            public double h_obs { get; set; }

            [Option("alpha", Required = false, Default = 0.05, HelpText = "State multiplier in ARMA model for acceleration")]
            public double alpha { get; set; }
            [Option("beta", Required = false, Default = 1.0, HelpText = "Noise multiplier in ARMA model for acceleration")]
            public double beta { get; set; }
            [Option("gamma", Required = false, Default = 0.5, HelpText = "Constant in ARMA model for acceleration")]
            public double gamma { get; set; }

            [Option('t', "train-count", Required = false, HelpText = "Number of trajectoris in the training set for dynamic models")]
            public int TrainCount { get; set; }

            [Option('e', "test-count", Required = false, HelpText = "Number of trajectoris in the test set for dynamic models")]
            public int TestCount { get; set; }

            [Option('N', "bundle-count", Required = false, Default = 1, HelpText = "Number of bundles of trajectoris in the test set for dynamic models")]
            public int BundleCount { get; set; }

            [Option('D', "Dummy", Required = false, Default = false, HelpText = "Do calculate dummuy filter")]
            public bool Dummy { get; set; }

            [Option('E', "EKF", Required = false, Default = false, HelpText = "Do calculate Extended Kalman Filter")]
            public bool EKF { get; set; }

            [Option('U', "UKF", Required = false, Default = false, HelpText = "Do calculate Unscented Kalman Filter, default parameters with no optimization")]
            public bool UKF { get; set; }
            [Option("UKF-file", Required = false, HelpText = "File name to save/load parameters of UKF filter")]
            public string UKFFileName { get; set; }

            [Option("UKF-integral-randomshoot", Required = false, Default = false, HelpText = "Do calculate Unscented Kalman Filter, integral parameter optimization with random shooting")]
            public bool UKFIntegralRandomShoot { get; set; }
            [Option("UKF-integral-randomshoot-file", Required = false, HelpText = "File name to save/load parameters of trained UKF filter")]
            public string UKFIntegralRandomShootFileName { get; set; }

            [Option("UKF-integral-NelderMead", Required = false, Default = false, HelpText = "Do calculate Unscented Kalman Filter, integral parameter optimization with Nelder-Mead method")]
            public bool UKFIntegralNelderMead { get; set; }
            [Option("UKF-integral-NelderMead-file", Required = false, HelpText = "File name to save/load parameters of trained UKF filter")]
            public string UKFIntegralNelderMeadFileName { get; set; }

            [Option("UKF-stepwise-randomshoot", Required = false, Default = false, HelpText = "Do calculate Unscented Kalman Filter, stepwise parameter optimization with random shooting")]
            public bool UKFStepwiseRandomShoot { get; set; }
            [Option("UKF-stepwise-randomshoot-file", Required = false, HelpText = "File name to save/load parameters of trained UKF filter")]
            public string UKFStepwiseRandomShootFileName { get; set; }

            [Option("UKF-stepwise-NelderMead", Required = false, Default = false, HelpText = "Do calculate Unscented Kalman Filter, stepwise parameter optimization with Nelder-Mead method")]
            public bool UKFStepwiseNelderMead { get; set; }
            [Option("UKF-stepwise-NelderMead-file", Required = false, HelpText = "File name to save/load parameters of trained UKF filter")]
            public string UKFStepwiseNelderMeadFileName { get; set; }

            [Option('C', "CMNF", Required = false, Default = false, HelpText = "Do calculate conditionally minimax nonlinear filter")]
            public bool CMNF { get; set; }
            [Option("CMNF-file", Required = false, HelpText = "File name to save/load parameters of trained CMNF filter")]
            public string CMNFFileName { get; set; }

            [Option('B', "BCMNF", Required = false, Default = false, HelpText = "Do calculate another conditionally minimax nonlinear filter")]
            public bool BCMNF { get; set; }
            [Option("BCMNF-file", Required = false, HelpText = "File name to save/load parameters of trained another CMNF filter")]
            public string BCMNFFileName { get; set; }

            [Option('M', "MCMNF", Required = false, Default = false, HelpText = "Do calculate modified conditionally minimax nonlinear filter")]
            public bool MCMNF { get; set; }

            [Option('l', "MCMNF-train-count", Required = false, Default = 0, HelpText = "Number of test set points for MCMNF wtepwise parameters fitting")]
            public int MCMNFTrainCount { get; set; }

            [Option('o', "output-folder", Required = false, HelpText = "Folder to store numeric results")]
            public string OutputFolder { get; set; }

            [Option('s', "scripts-folder", Required = false, HelpText = "Folder where the python scripts are stored")]
            public string ScriptsFolder { get; set; }

            [Option("bulk", Required = false, Default = false, HelpText = "Bulk output of the state process trajectories bundle")]
            public bool Bulk { get; set; }

            [Option('G', "generate-samples", Required = false, Default = 0, HelpText = "Number of single trajectories to generate")]
            public int SamplesCount { get; set; }

            [Option('P', "parallel", Required = false, Default = false, HelpText = "Parallel bundles calculation (on to save time, off to save memory in case of multiple large bundles)")]
            public bool Parallel { get; set; }

            [Option('d', "degree-parallelism", Required = false, Default = 1, HelpText = "Maximum degreee of parallelism")]
            public int ParallelismDegree { get; set; }

            [Option('L', "load", Required = false, Default = false, HelpText = "Load filter parameters from files")]
            public bool Load { get; set; }

            [Option('S', "Save", Required = false, Default = false, HelpText = "Save trained filter parameters to files for future use")]
            public bool Save { get; set; }

            [Option("sift", Required = false, Default = false, HelpText = "Sift diverged estimates before estimate error statistic calculation")]
            public bool Sift { get; set; }

            [Option("sift-bound", Required = false, Default = 1000.0, HelpText = "Sift bound estimate to be considered as diverged: Math.Sqrt(dx*dx + dy*dy) > bound, where dx = x_hat - x, dy = y_hat - y")]
            public double SiftBound { get; set; }

            [Option("skip-filter", Required = false, Default = false, HelpText = "Skip filter calculation step")]
            public bool Skip { get; set; }

            [Option('A', "aggregate", Required = false, Default = false, HelpText = "Aggregate the previously generated data")]
            public bool Aggregate { get; set; }

            [Option("no-bin", Required = false, Default = false, HelpText = "Do not save the calculated data to reuse later")]
            public bool NoBin { get; set; }

            [Option("no-text", Required = false, Default = false, HelpText = "Do not save the calculated data to text file")]
            public bool NoText { get; set; }

            [Option("no-python", Required = false, Default = false, HelpText = "Do not process the generated data")]
            public bool NoPython { get; set; }

        }

        static void Main(string[] args)
        {
            CommandLine.Parser.Default.ParseArguments<Options>(args).WithParsed<Options>(opts => Run(opts, args));
        }

        static void Run(Options o, string[] args)
        {
            if (string.IsNullOrWhiteSpace(o.ScriptsFolder))
                o.ScriptsFolder = "..\\..\\..\\OutputScripts\\";

            Directory.CreateDirectory(o.OutputFolder);
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(Path.Combine(o.OutputFolder, "parameters.txt"), true))
            {
                outputfile.WriteLine($"{DateTime.Now}\t{string.Join(" ", args)}");
                outputfile.Close();
            }

            #region 3d model Ienkaran Arasaratnam, Simon Haykin, and Tom R. Hurd

            //double sigma1 = Math.Sqrt(0.2);
            //double sigma2 = 7.0 * 1e-3;
            //double sigma_r = 50;
            //double sigma_th = 0.1;
            //double sigma_ph = 0.1;


            //Vector<double> mEta = Exts.Vector(1000, 0, 2650, 150, 200, 0, 1.0);
            //Matrix<double> dEta = Exts.Diag(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);

            //Vector<double> mW = Exts.Vector(0, 0, 0, 0, 0, 0, 0);
            //Matrix<double> dW = Exts.Diag(1e1, Math.Pow(sigma1, 2), 1e1, Math.Pow(sigma1, 2), 1e1, Math.Pow(sigma1, 2), Math.Pow(sigma2, 2));

            //Vector<double> mNu = Exts.Vector(0, 0, 0);
            //Matrix<double> dNu = Exts.Diag(Math.Pow(sigma_r, 2), Math.Pow(sigma_th, 2), Math.Pow(sigma_ph, 2));


            //Func<Vector<double>, Vector<double>> Phi1 = (x) => Exts.Vector(x[1], -x[6] * x[3], x[3], x[6] * x[1], x[5], 0, 0);
            //Func<Matrix<double>> Phi2 = () => Exts.Diag(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
            //Func<Vector<double>, Vector<double>> Psi1 = (x) => Utils.cart2sphere(Exts.Vector(x[0], x[2], x[4]));
            //Func<Matrix<double>> Psi2 = () => Exts.Diag(1.0, 1.0, 1.0);

            ////Vector<double> X_R2 = Exts.Vector(0, 0);
            ////Vector<double> mNu = Exts.Vector(0, 0, 0, 0);
            ////Matrix<double> dNu = Exts.Diag(Math.Pow(0.1 * Math.PI / 180, 2), Math.Pow(50, 2), Math.Pow(0.1 * Math.PI / 180, 2), Math.Pow(50, 2));
            ////Func<double, Vector<double>, Vector<double>> Psi1 = (s, x) => Exts.Stack(Utils.cart2pol(Exts.Vector(x[0], x[1]) - X_R1), Utils.cart2pol(Exts.Vector(x[0], x[1]) - X_R2));
            ////Func<double, Vector<double>, Matrix<double>> Psi2 = (s, x) => Exts.Diag(1.0, 1.0, 1.0, 1.0);

            //RandomVector<Normal> NormalW = new RandomVector<Normal>(mW, dW);
            //RandomVector<Normal> NormalNu = new RandomVector<Normal>(mNu, dNu);
            //RandomVector<Normal> NormalEta = new RandomVector<Normal>(mEta, dEta);

            //Func<Vector<double>> W;
            //Func<Vector<double>> Nu;
            //Func<Vector<double>> X0;
            //W = () => NormalW.Sample();
            //Nu = () => NormalNu.Sample();
            ////X0 = () => NormalEta.Sample();
            //X0 = () => mEta;

            //double h_state = 0.01;
            //double h_obs = 0.01;
            //double T = 1.0 + h_state / 2;

            #endregion


            #region 2d model with acceleration

            double Alpha_n = o.alpha; // 0.05;
            double Beta_n = o.beta; //1.0;
            //double Gamma_n = 2.5;
            double Gamma_n = o.gamma; //0.5;


            //Vector<double> mEta = Exts.Vector(30000, 40000, 400, 10 * Math.PI / 180, Gamma_n / Alpha_n);
            Vector<double> mEta = Exts.Vector(0, 25000, 400, 10 * Math.PI / 180, Gamma_n / Alpha_n);
            Matrix<double> dEta = Exts.Diag(Math.Pow(2000, 2), Math.Pow(2000, 2), Math.Pow(115, 2), Math.Pow(15 * Math.PI / 180, 2), Math.Pow(Beta_n, 2) / 2 / Alpha_n);

            Vector<double> mW = Exts.Vector(0, 0, 0, 0, 0);
            Matrix<double> dW = Exts.Diag(0, 0, 0, 0, Math.Pow(Beta_n, 2));

            Vector<double> X_R1 = Exts.Vector(-10000, 10000);


            Func<Vector<double>, Vector<double>> Phi1 = (x) => Exts.Vector(x[2] * Math.Cos(x[3]), x[2] * Math.Sin(x[3]), 0, x[4] / x[2], -Alpha_n * x[4] + Gamma_n);
            Func<Matrix<double>> Phi2 = () => Exts.Diag(0, 0, 0, 0, 1.0);

            //Func<Vector<double>, Vector<double>> Psi1 = (x) => Utils.cart2pol(Exts.Vector(x[0], x[1]) - X_R1);
            //Func<Matrix<double>> Psi2 = () => Exts.Diag(1.0, 1.0);
            Vector<double> X_R2 = Exts.Vector(0, 0);
            Vector<double> mNu = Exts.Vector(0, 0, 0, 0);
            Matrix<double> dNu = Exts.Diag(Math.Pow(0.1 * Math.PI / 180, 2), Math.Pow(50, 2), Math.Pow(0.1 * Math.PI / 180, 2), Math.Pow(50, 2));
            //Matrix<double> dNu = Exts.Diag(Math.Pow(0.01 * Math.PI / 180, 2), Math.Pow(0.01, 2), Math.Pow(0.01 * Math.PI / 180, 2), Math.Pow(0.01, 2));
            Func<Vector<double>, Vector<double>> Psi1 = (x) => Exts.Stack(Utils.cart2pol(Exts.Vector(x[0], x[1]) - X_R1), Utils.cart2pol(Exts.Vector(x[0], x[1]) - X_R2));
            Func<Matrix<double>> Psi2 = () => Exts.Diag(1.0, 1.0, 1.0, 1.0);

            Normal NormalW = new Normal(mW[4], Math.Sqrt(dW[4, 4]));
            RandomVector<Normal> NormalNu = new RandomVector<Normal>(mNu, dNu);
            RandomVector<Normal> NormalEta = new RandomVector<Normal>(mEta, dEta);
            ContinuousUniform UniformEtaV = new ContinuousUniform(200, 600);

            Func<Vector<double>> W;
            Func<Vector<double>> Nu;
            Func<Vector<double>> X0;
            W = () => Exts.Vector(0, 0, 0, 0, NormalW.Sample());
            Nu = () => NormalNu.Sample();
            //W = () => Exts.Vector(0, 0, 0, 0, 0);
            //Nu = () => Exts.Vector(0, 0, 0, 0);
            X0 = () =>
            {
                var x = NormalEta.Sample();
                x[2] = UniformEtaV.Sample();
                return x;
            };
            //X0 = () => mEta;

            //discretize


            double h_state = o.h_state;//0.01;
            double h_obs = o.h_obs; //1.0;
            double T = o.T + h_state / 2;


            //Func<int, Vector<double>, Vector<double>> Phi1_discr = (i, x) => x + h_obs * Phi1(x);
            Func<int, Vector<double>, Vector<double>> Phi1_discr = (i, x) =>
            {
                Vector<double> xx = x;
                for (double s = h_state; s < h_obs + h_state / 2; s += h_state)
                {
                    xx += h_state * Phi1(xx);
                }
                return xx;
            };
            Func<int, Vector<double>, Matrix<double>> Phi2_discr = (i, x) => Math.Sqrt(h_obs) * Phi2();
            Func<int, Vector<double>, Vector<double>> Psi1_discr = (i, x) => Psi1(x);
            Func<int, Vector<double>, Matrix<double>> Psi2_discr = (i, x) => Psi2();

            Func<Vector<double>, Matrix<double>> dPhi = (x) => Matrix<double>.Build.Dense(5, 5, new double[5 * 5] {
                0, 0, 0, 0, 0,
                0, 0, 0, 0, 0,
                Math.Cos(x[3]), Math.Sin(x[3]), 0, -x[4] / Math.Pow(x[2], 2), 0,
                -x[2] * Math.Sin(x[3]), x[2] * Math.Cos(x[3]), 0, 0, 0,
                0, 0, 0, 1.0 / x[2], -Alpha_n
            }); // Column-major order
            Func<Vector<double>, Matrix<double>> dPsi = (x) => Matrix<double>.Build.Dense(4, 5, new double[5 * 4] {
                -x[1] / (x[0]*x[0] + x[1] * x[1]), x[0] / Math.Sqrt(x[0] * x[0] + x[1] * x[1]), -x[1]/ (x[0] * x[0] + x[1] * x[1]), x[0]/ Math.Sqrt(x[0] * x[0] + x[1] * x[1]),
                 x[0] / (x[0]*x[0] + x[1] * x[1]), x[1] / Math.Sqrt(x[0] * x[0] + x[1] * x[1]),  x[0]/ (x[0] * x[0] + x[1] * x[1]), x[1]/ Math.Sqrt(x[0] * x[0] + x[1] * x[1]),
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0
            }); // Column-major order

            Func<int, Vector<double>, Matrix<double>, (Vector<double>, Matrix<double>)> Ricatti = (i, x, P) =>
            {
                Vector<double> xx = x;
                Matrix<double> PP = P;
                for (double s = h_state; s < h_obs + h_state / 2; s += h_state)
                {
                    PP += h_state * (dPhi(xx) * PP + PP * dPhi(xx).Transpose() + Phi2() * dW * Phi2().Transpose());
                    xx += h_state * (Phi1(xx) + Phi2() * mW);
                }
                return (xx, PP);
            };

            Func<int, Vector<double>, Vector<double>, Matrix<double>, (Vector<double>, Matrix<double>)> DummyEstimate = (i, y, x, P) =>
            {
                Vector<double> xx = x;
                Matrix<double> PP = P;
                for (double s = h_state; s < h_obs + h_state / 2; s += h_state)
                {
                    PP += h_state * (dPhi(xx) * PP + PP * dPhi(xx).Transpose() + Phi2() * dW * Phi2().Transpose());
                    xx += h_state * (Phi1(xx) + Phi2() * mW);
                }
                return ((0.5 * (Utils.pol2cart(Exts.Vector(y[0], y[1])) + X_R1 + Utils.pol2cart(Exts.Vector(y[2], y[3])) + X_R2)).Stack(Exts.Vector(xx[2], xx[3], xx[4])), PP);
            };

            #endregion

            int N = (int)(T / h_obs);

            //string folder = "d:\\results\\cont\\";
            //Directory.CreateDirectory(folder);
            string CMNFFileName = Path.Combine(o.OutputFolder, "cmnf.params");
            if (!string.IsNullOrWhiteSpace(o.CMNFFileName)) CMNFFileName = o.CMNFFileName;
            string BCMNFFileName = Path.Combine(o.OutputFolder, "bcmnf.params");
            if (!string.IsNullOrWhiteSpace(o.BCMNFFileName)) BCMNFFileName = o.BCMNFFileName;

            string UKFFileName = Path.Combine(o.OutputFolder, "ukf.params");
            if (!string.IsNullOrWhiteSpace(o.UKFFileName)) UKFFileName = o.UKFFileName;
            string UKFOptStepwiseNMFileName = Path.Combine(o.OutputFolder, "ukfoptstepwiseNM.params");
            if (!string.IsNullOrWhiteSpace(o.UKFStepwiseNelderMeadFileName)) UKFOptStepwiseNMFileName = o.UKFStepwiseNelderMeadFileName;
            string UKFOptIntegralNMFileName = Path.Combine(o.OutputFolder, "ukfoptintegralNM.params");
            if (!string.IsNullOrWhiteSpace(o.UKFIntegralNelderMeadFileName)) UKFOptIntegralNMFileName = o.UKFIntegralNelderMeadFileName;
            string UKFOptStepwiseRandFileName = Path.Combine(o.OutputFolder, "ukfoptstepwiserand.params");
            if (!string.IsNullOrWhiteSpace(o.UKFStepwiseRandomShootFileName)) UKFOptStepwiseRandFileName = o.UKFStepwiseRandomShootFileName;
            string UKFOptIntegralRandFileName = Path.Combine(o.OutputFolder, "ukfoptintegralrand.params");
            if (!string.IsNullOrWhiteSpace(o.UKFIntegralRandomShootFileName)) UKFOptIntegralRandFileName = o.UKFIntegralRandomShootFileName;
            List<(FilterType, string)> filters = new List<(FilterType, string)>();
            if (o.CMNF) filters.Add((FilterType.CMNF, CMNFFileName));
            if (o.BCMNF) filters.Add((FilterType.BCMNF, BCMNFFileName));
            if (o.MCMNF) filters.Add((FilterType.BCMNF, string.Empty));

            if (o.UKF) filters.Add((FilterType.UKFNoOptimization, UKFFileName));
            if (o.UKFStepwiseNelderMead) filters.Add((FilterType.UKFStepwise, UKFOptStepwiseNMFileName));
            if (o.UKFIntegralNelderMead) filters.Add((FilterType.UKFIntegral, UKFOptIntegralNMFileName));
            if (o.UKFStepwiseRandomShoot) filters.Add((FilterType.UKFStepwiseRandomShoot, UKFOptStepwiseRandFileName));
            if (o.UKFIntegralRandomShoot) filters.Add((FilterType.UKFIntegralRandomShoot, UKFOptIntegralRandFileName));
            if (o.EKF) filters.Add((FilterType.EKF, string.Empty));
            if (o.Dummy) filters.Add((FilterType.Dummy, string.Empty));

                TestEnvironmentVector testEnv = new TestEnvironmentVector()
            {
                TestName = "Target tracking",
                TestFileName = "TargetTracking",
                Phi1 = Phi1_discr,
                Phi2 = Phi2_discr,
                Psi1 = Psi1_discr,
                Psi2 = Psi2_discr,
                dPhi = (i, x) => dPhi(x),
                dPsi = (i, x) => dPsi(x),
                Xi = (i, x) => Phi1_discr(i, x) + Phi2_discr(i, x) * mW,
                //Zeta = (i, x, y, k) => (y - Psi1_discr(i, x) - Psi2_discr(i, x) * mNu).Stack(Utils.pol2cart(Exts.Vector(y[0], y[1]))+X_R1).Stack(Utils.pol2cart(Exts.Vector(y[2], y[3]))+X_R2),
                Zeta = (i, x, y, k) => (y - Psi1_discr(i, x) - Psi2_discr(i, x) * mNu),
                Alpha = (i, x) => Phi1_discr(i, x) + Phi2_discr(i, x) * mW,
                Gamma = (i, x, y) => (y).Stack(Utils.pol2cart(Exts.Vector(y[0], y[1])) + X_R1).Stack(Utils.pol2cart(Exts.Vector(y[2], y[3])) + X_R2),
                //Gamma = (i, x, y) => y - Psi1_discr(i, x) - Psi2_discr(i, x) * mNu,
                W = (i) => W(),
                Nu = (i) => Nu(),
                DW = dW,
                DNu = dNu,
                X0 = () => X0(),
                X0Hat = mEta,
                DX0Hat = dEta,
                Predict = Ricatti,
                DummyEstimate = DummyEstimate
            };

            Func<DiscreteVectorModel> ModelGenerator = () =>
            {
                DiscreteVectorModel model = null;

                int n = 0;
                double h_tolerance = h_state / 2.0;
                double t_nextobservation = h_obs;
                Vector<double> State = X0();
                Vector<double> Obs;

                model = new DiscreteVectorModel(Phi1_discr, Phi2_discr, Psi1_discr, Psi2_discr, (i) => W(), (i) => Nu(), X0(), true);
                //for (double s = 0; s < T; s += h_obs)
                //{
                //    model.Step();
                //}
                for (double s = h_state; s < T; s += h_state)
                {
                    if (s > 0)
                    {
                        State = State + h_state * Phi1(State) + Math.Sqrt(h_state) * Phi2() * W();
                    }
                    if (Math.Abs(s - t_nextobservation) < h_tolerance)
                    {
                        Obs = Psi1(State) + Psi2() * Nu();
                        t_nextobservation += h_obs;

                        n++;
                        model.Trajectory.Add(n, new Vector<double>[] { State, Obs });
                    }
                }
                return model;
            };

            //bool doCalculateFilter = true; 
            //testEnv.Initialize(N, 100000, 100, folder, filters, doCalculateFilter, !doCalculateFilter, ModelGenerator);
            //testEnv.Sifter = (x) => Math.Sqrt(x[0] * x[0] + x[1] * x[1]) > 500.0;
            //testEnv.GenerateBundleSamples(50, 1000, "d:\\results\\cont\\");
            //for (int i = 0; i < 10; i++)
            //{
            //    testEnv.GenerateOne("folder", i);
            //}

            //testEnv.GenerateOne(folder);
            //testEnv.RunScript(Path.Combine("..\\..\\..\\OutputScripts\\", "estimate_sample.py"), folder);
            //testEnv.RunScript(Path.Combine("..\\..\\..\\OutputScripts\\", "trajectory.py"), folder);


            //testEnv.GenerateBundles(1, 10, folder, true, 8);
            //testEnv.RunScript(Path.Combine("..\\..\\..\\OutputScripts\\", "estimate_statistics.py"), folder);

            //ContinuousVectorModel model = new ContinuousVectorModel(h_state, h_obs, Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
            ////for (double s = 0; s < T; s += h_state)
            ////{
            ////    model.Step();
            ////}
            //for (int n = 0; n * h_obs < T+h_state; n++)
            //{
            //    model.StepObs();
            //}
            //model.SaveTrajectory("d:\\results\\cont\\state.txt");
            //model.SaveObservations("d:\\results\\cont\\Obs.txt");

            if (o.Bulk)
                testEnv.GenerateBundleSamples(o.T, o.TrainCount, o.OutputFolder);
            else
            {
                testEnv.Initialize(o.T, o.TrainCount, o.MCMNFTrainCount, o.OutputFolder, filters, o.Save, o.Load, ModelGenerator);
                if (o.Sift) testEnv.Sifter = (x) => Math.Sqrt(x[0] * x[0] + x[1] * x[1]) > o.SiftBound;
                if (o.Aggregate)
                {
                    testEnv.Aggregate(o.OutputFolder, o.OutputFolder, !o.NoBin, !o.NoText);
                }
                if (!o.Skip)
                {
                    if (o.SamplesCount == 0)
                    {
                        testEnv.GenerateBundles(o.BundleCount, o.TestCount, o.OutputFolder, o.Parallel, o.ParallelismDegree, !o.NoBin, !o.NoText);
                        if (!o.NoPython)
                        {
                            testEnv.RunScript(Path.Combine(o.ScriptsFolder, "estimate_statistics.py"), o.OutputFolder);
                        }

                    }
                    else
                    {
                        if (o.SamplesCount == 1)
                        {
                            testEnv.GenerateOne(o.OutputFolder);
                            if (!o.NoPython)
                            {
                                testEnv.RunScript(Path.Combine(o.ScriptsFolder, "estimate_sample.py"), o.OutputFolder);
                                testEnv.RunScript(Path.Combine(o.ScriptsFolder, "trajectory.py"), o.OutputFolder);
                            }
                        }
                        else
                        {
                            for (int i = 0; i < o.SamplesCount; i++)
                            {
                                testEnv.GenerateOne(o.OutputFolder, i);
                            }
                        }
                    }
                }
            }

        }
    }
}
