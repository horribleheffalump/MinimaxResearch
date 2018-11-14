using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using NonlinearSystem;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TargetTrackingTest
{
    class Program
    {
        static void Main(string[] args)
        {
            Vector<double> mW = Exts.Vector(0); Matrix<double> dW = Exts.Diag(1);
            Vector<double> mNu = Exts.Vector(0); Matrix<double> dNu = Exts.Diag(1);
            Vector<double> mEta = Exts.Vector(0); Matrix<double> dEta = Exts.Diag(1);
            Func<double, Vector<double>, Vector<double>> Phi1 = (s, x) => Exts.Vector(1);
            Func<double, Vector<double>, Matrix<double>> Phi2 = (s, x) => Exts.Diag(1);
            Func<double, Vector<double>, Vector<double>> Psi1 = (s, x) => Exts.Vector(x[0]);
            Func<double, Vector<double>, Matrix<double>> Psi2 = (s, x) => Exts.Diag(1);

            Normal[] NormalW = new Normal[1] { new Normal(mW[0], Math.Sqrt(dW[0, 0])) };
            Normal[] NormalNu = new Normal[1] { new Normal(mNu[0], Math.Sqrt(dNu[0, 0])) };
            Normal[] NormalEta = new Normal[1] { new Normal(mEta[0], Math.Sqrt(dEta[0, 0]))}; ;

            Func<double, Vector<double>> W;
            Func<double, Vector<double>> Nu;
            Func<Vector<double>> X0;
            W = (s) => Exts.Vector(NormalW[0].Sample());
            Nu = (s) => Exts.Vector(NormalNu[0].Sample());
            X0 = () => Exts.Vector(NormalEta[0].Sample());

            double h_state = 0.01;
            double h_obs = 1.0;
            double T = 10;
            ContinuousVectorModel model = new ContinuousVectorModel(h_state, h_obs, Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
            for (double s = 0; s < T; s+=h_state)
            {
                model.Step();
            }

            model.SaveTrajectory("d:\\results\\cont\\state.txt");
            model.SaveObservations("d:\\results\\cont\\Obs.txt");
        }
    }
}
