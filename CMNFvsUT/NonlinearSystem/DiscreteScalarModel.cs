using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;

namespace NonlinearSystem
{
    public class DiscreteScalarModel
    {
        private Func<double, double> Phi1;
        private Func<double, double> Phi2;
        private Func<double, double> Psi;
        private Func<int, double> W;
        private Func<int, double> Nu;
        private int t;
        public double State;
        public double Obs;
        private bool doSave;


      

        public Dictionary<int,double[]> Trajectory;

        public DiscreteScalarModel(Func<double, double> phi1, Func<double, double> phi2, Func<double, double> psi, Func<int, double> w, Func<int, double> nu, double X0, bool saveHistory = false)
        {
            t = 0;
            Phi1 = phi1;
            Phi2 = phi2;
            Psi = psi;
            W = w;
            Nu = nu;

            doSave = saveHistory;

            Trajectory = new Dictionary<int, double[]>();

            State = X0;
            //if (doSave) Trajectory.Add(t, State);
            //Obs = Psi(State) + Nu(t);
        }

        public double Step()
        {
            if (t>0)
            {
                State = Phi1(State) + Phi2(State) * W(t);
            }
            Obs = Psi(State) + Nu(t);
            if (doSave) Trajectory.Add(t, new double[] { State, Obs });
            t++;
            return Obs;
        }

        public void SaveTrajectory(string path)
        {
            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(path))
            {
                foreach (var x in Trajectory.OrderBy(s => s.Key))
                {
                    outputfile.WriteLine(string.Format(provider, "{0} {1}", x.Key, x.Value));
                }
                outputfile.Close();
            }
        }

    }
}
