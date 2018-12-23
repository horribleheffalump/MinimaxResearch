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
        static void Main(string[] args)
        {
            #region simple 2d model

            //Vector<double> mEta = Exts.Vector(100000, 200000);
            //Matrix<double> dEta = Exts.Diag(1.0, 1.0);

            //Vector<double> mW = Exts.Vector(0, 0);
            //Matrix<double> dW = Exts.Diag(1e-10, 1e-10);

            //Vector<double> X_R1 = Exts.Vector(20000, 0);
            //Vector<double> mNu = Exts.Vector(0, 0);
            //Matrix<double> dNu = Exts.Diag(Math.Pow(1e-10 * Math.PI / 180, 2), Math.Pow(1e-10, 2));


            //Func<Vector<double>, Vector<double>> Phi1 = (x) => Exts.Vector(x[0], x[1]);
            //Func<Matrix<double>> Phi2 = () => Exts.Diag(1.0, 1.0);
            //Func<Vector<double>, Vector<double>> Psi1 = (x) => Exts.Vector(1,1);
            //Func<Matrix<double>> Psi2 = () => Exts.Diag(1.0, 1.0);


            //RandomVector<Normal> NormalW = new RandomVector<Normal>(mW, dW);
            //RandomVector<Normal> NormalNu = new RandomVector<Normal>(mNu, dNu);


            //Func<Vector<double>> W;
            //Func<Vector<double>> Nu;
            //Func<Vector<double>> X0;
            //W = () => NormalW.Sample();
            //Nu = () => NormalNu.Sample();
            ////X0 = () => NormalEta.Sample();
            //X0 = () => mEta;

            //double h_state = 0.1;
            //double h_obs = 0.1;
            //double T = 1.0 + h_state / 2.0;

            #endregion

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

            //double Alpha_t = 0.1;
            //double Beta_t = 2;
            //double Gamma_t = -0.005;
            //double Alpha_n = 0.05;
            double Alpha_n = 0.05;
            double Beta_n = 5.0;
            double Gamma_n = 2.0;


            Vector<double> mEta = Exts.Vector(30000, 40000, 400, 10 * Math.PI / 180, Gamma_n / Alpha_n);
            Matrix<double> dEta = Exts.Diag(Math.Pow(2000, 2), Math.Pow(2000, 2), Math.Pow(115, 2), Math.Pow(15 * Math.PI / 180, 2), Math.Pow(Beta_n, 2) / 2 / Alpha_n);
            //Matrix<double> dEta = Exts.Diag(0.0, 0.0, 0.0, 0.0, 0.0);

            Vector<double> mW = Exts.Vector(0, 0, 0, 0, 0);
            Matrix<double> dW = Exts.Diag(0, 0, 0, 0, Math.Pow(Beta_n, 2));

            Vector<double> X_R1 = Exts.Vector(20000, 0);
            //Vector<double> mNu = Exts.Vector(0, 0);
            //Matrix<double> dNu = Exts.Diag(Math.Pow(0.1 * Math.PI / 180, 2), Math.Pow(50, 2));


            Func<Vector<double>, Vector<double>> Phi1 = (x) => Exts.Vector(x[2] * Math.Cos(x[3]), x[2] * Math.Sin(x[3]), 0, x[4] / x[2], -Alpha_n * x[4] + Gamma_n);
            Func<Matrix<double>> Phi2 = () => Exts.Diag(0, 0, 0, 0, 1.0);
            //Func<Vector<double>, Vector<double>> Psi1 = (x) => Utils.cart2pol(Exts.Vector(x[0], x[1]) - X_R1);
            //Func<Matrix<double>> Psi2 = () => Exts.Diag(1.0, 1.0);

            Vector<double> X_R2 = Exts.Vector(0, 0);
            Vector<double> mNu = Exts.Vector(0, 0, 0, 0);
            Matrix<double> dNu = Exts.Diag(Math.Pow(0.1 * Math.PI / 180, 2), Math.Pow(50, 2), Math.Pow(0.1 * Math.PI / 180, 2), Math.Pow(50, 2));
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
            X0 = () =>
            {
                var x = NormalEta.Sample();
                x[2] = UniformEtaV.Sample();
                return x;
            };
            //X0 = () => mEta;

            double h_state = 0.1;
            double h_obs = 0.1;
            double T = 20.0 + h_state / 2;

            #endregion


            //discretize

            int N = (int)(T / h_obs);

            string folder = "d:\\results\\cont_EKF\\";

            string CMNFFileName = Path.Combine(folder, "cmnf.params");
            string UKFFileName = Path.Combine(folder, "ukf.params");
            List<(FilterType, string)> filters = new List<(FilterType, string)>
            {
                (FilterType.CMNF, CMNFFileName),
                //(FilterType.MCMNF, string.Empty),
                (FilterType.UKFNoOptimization, UKFFileName),
                (FilterType.EKF, string.Empty)
            };

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
                Zeta = (i, x, y, k) => (y - Psi1_discr(i, x) - Psi2_discr(i, x) * mNu).Stack(Utils.pol2cart(Exts.Vector(y[0], y[1]))).Stack(Utils.pol2cart(Exts.Vector(y[2], y[3]))),
                W = (i) => W(),
                Nu = (i) => Nu(),
                DW = dW,
                DNu = dNu,
                X0 = () => X0(),
                X0Hat = mEta,
                DX0Hat = dEta
            };

            Func<DiscreteVectorModel> ModelGenerator = () =>
            {
                int n = 0;
                double h_tolerance = h_state / 2.0;
                double t_nextobservation = h_obs;
                Vector<double> State = X0();
                Vector<double> Obs;
                DiscreteVectorModel model = new DiscreteVectorModel(Phi1_discr, Phi2_discr, Psi1_discr, Psi2_discr, (i) => W(), (i) => Nu(), X0(), true);
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


            bool doCalculateFilter = true;
            testEnv.Initialize(N, 1000, 100, "d:\\results\\cont_EKF\\", filters, doCalculateFilter, !doCalculateFilter, ModelGenerator);
            //testEnv.GenerateBundleSamples(50, 1000, "d:\\results\\cont\\");
            testEnv.GenerateOne("d:\\results\\cont_EKF\\");
            testEnv.RunScript(Path.Combine("..\\..\\..\\OutputScripts\\", "estimate_sample.py"), folder);
            testEnv.RunScript(Path.Combine("..\\..\\..\\OutputScripts\\", "trajectory.py"), folder);

            
            testEnv.GenerateBundles(1, 1000, "d:\\results\\cont_EKF\\", false, 4);
            testEnv.RunScript(Path.Combine("..\\..\\..\\OutputScripts\\", "estimate_statistics.py"), folder);

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
        }
    }
}
