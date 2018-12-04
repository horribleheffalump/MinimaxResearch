using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using NonlinearSystem;
using System;
using System.Collections.Generic;
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
            //double Alpha_t = 0.1;
            //double Beta_t = 2;
            //double Gamma_t = -0.005;
            double Alpha_n = 0.05;
            double Beta_n = 1;
            double Gamma_n = 0.02;


            Vector<double> mEta = Exts.Vector(100000, 200000, 450, 10 * Math.PI / 180, Gamma_n / Alpha_n);
            Matrix<double> dEta = Exts.Diag(Math.Pow(2000, 2), Math.Pow(2000, 2), Math.Pow(100, 2), Math.Pow(15 * Math.PI / 180, 2), Math.Pow(Beta_n, 2) / 2 / Alpha_n);

            Vector<double> mW = Exts.Vector(0, 0, 0, 0, 0);
            Matrix<double> dW = Exts.Diag(1.0, 1.0, 1.0, 1.0, Math.Pow(Beta_n, 2));

            Vector<double> X_R1 = Exts.Vector(20000, 0);
            Vector<double> mNu = Exts.Vector(0, 0);
            Matrix<double> dNu = Exts.Diag(Math.Pow(0.1 * Math.PI / 180, 2), Math.Pow(50, 2));


            Func<double, Vector<double>, Vector<double>> Phi1 = (s, x) => Exts.Vector(x[2] * Math.Cos(x[3]), x[2] * Math.Sin(x[3]), 0, x[4] / x[2], -Alpha_n * x[4] + Gamma_n);
            Func<double, Matrix<double>> Phi2 = (s) => Exts.Diag(0, 0, 0, 0, 1.0);
            Func<double, Vector<double>, Vector<double>> Psi1 = (s, x) => Utils.cart2pol(Exts.Vector(x[0], x[1]) - X_R1);
            Func<double, Vector<double>, Matrix<double>> Psi2 = (s, x) => Exts.Diag(1.0, 1.0);
            
            //Vector<double> X_R2 = Exts.Vector(0, 0);
            //Vector<double> mNu = Exts.Vector(0, 0, 0, 0);
            //Matrix<double> dNu = Exts.Diag(Math.Pow(0.1 * Math.PI / 180, 2), Math.Pow(50, 2), Math.Pow(0.1 * Math.PI / 180, 2), Math.Pow(50, 2));
            //Func<double, Vector<double>, Vector<double>> Psi1 = (s, x) => Exts.Stack(Utils.cart2pol(Exts.Vector(x[0], x[1]) - X_R1), Utils.cart2pol(Exts.Vector(x[0], x[1]) - X_R2));
            //Func<double, Vector<double>, Matrix<double>> Psi2 = (s, x) => Exts.Diag(1.0, 1.0, 1.0, 1.0);

            Normal NormalW = new Normal(mW[4], Math.Sqrt(dW[4, 4]));
            RandomVector<Normal> NormalNu = new RandomVector<Normal>(mNu, dNu);
            RandomVector<Normal> NormalEta = new RandomVector<Normal>(mEta, dEta);

            Func<Vector<double>> W;
            Func<Vector<double>> Nu;
            Func<Vector<double>> X0;
            W = () => Exts.Vector(0, 0, 0, 0, NormalW.Sample());
            Nu = () => NormalNu.Sample();
            X0 = () => NormalEta.Sample();

            double h_state = 0.01;
            double h_obs = 1.0;
            double T = 600;



            //discretize

            int N = (int)(T / h_obs);

            string CMNFFileName = "d:\\results\\cont\\cmnf.params";
            string UKFFileName = "d:\\results\\cont\\ukf.params";
            List<(FilterType, string)> filters = new List<(FilterType, string)>
            {
                (FilterType.CMNF, CMNFFileName),
                (FilterType.MCMNF, string.Empty)//,
                //(FilterType.UKFIntegralRandomShoot, UKFFileName)
            };

            Func<int, Vector<double>, Vector<double>> Phi1_discr = (i, x) =>
            {
                double tt = i * h_obs;
                Vector<double> xx = x;
                while (tt < (i + 1) * h_obs + h_state)
                {
                    xx += h_state * Phi1(tt, xx);
                    tt += h_state;
                }
                return xx;
            };
            Func<int, Vector<double>, Matrix<double>> Phi2_discr = (i, x) => Math.Sqrt(h_obs) * Phi2(h_obs * i);
            Func<int, Vector<double>, Vector<double>> Psi1_discr = (i, x) => Psi1(h_obs * i, x);
            Func<int, Vector<double>, Matrix<double>> Psi2_discr = (i, x) => Psi2(h_obs * i, x);

            TestEnvironmentVector testEnv = new TestEnvironmentVector()
            {
                TestName = "Target tracking",
                TestFileName = "TargetTracking",
                Phi1 = Phi1_discr,
                Phi2 = Phi2_discr,
                Psi1 = Psi1_discr,
                Psi2 = Psi2_discr,
                Xi = (i, x) => Phi1_discr(i, x) + Phi2_discr(i, x) * mW,
                Zeta = (i, x, y, k) => y - Psi1_discr(i, x) - Psi2_discr(i, x) * mNu,
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
                Vector<double> Obs = Psi1(n, State) + Psi2(n, State) * Nu();
                DiscreteVectorModel model = new DiscreteVectorModel(Phi1_discr, Phi2_discr, Psi1_discr, Psi2_discr, (i) => W(), (i) => Nu(), X0(), true);
                for (double s = 0; s < T; s += h_state)
                {
                    if (s > 0)
                    {
                        State += h_state * Phi1(s, State) + Math.Sqrt(h_state) * Phi2(s) * W();
                    }
                    if (Math.Abs(s - t_nextobservation) < h_tolerance)
                    {
                        Obs = Psi1(s, State) + Psi2(s, State) * Nu();
                        t_nextobservation += h_obs;

                        model.Trajectory.Add(n, new Vector<double>[] { State, Obs });
                        n++;
                    }
                }
                return model;
            };

            bool doCalculateFilter = false;
            testEnv.Initialize(N, 1000, 100000, "d:\\results\\cont\\", filters, doCalculateFilter, !doCalculateFilter, ModelGenerator);
            //testEnv.GenerateBundles(1, 1, "d:\\results\\cont\\", false);
            testEnv.GenerateOne("d:\\results\\cont\\");

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
