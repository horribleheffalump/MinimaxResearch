using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NonlinearSystem
{
    public class ContinuousVectorModel
    {
        private Func<double, Vector<double>, Vector<double>> Phi1; // Phi1(t, X)
        private Func<double, Vector<double>, Matrix<double>> Phi2;
        private Func<double, Vector<double>, Vector<double>> Psi1;
        private Func<double, Vector<double>, Matrix<double>> Psi2;
        private Func<double, Vector<double>> W;
        private Func<double, Vector<double>> Nu;
        private double t;
        public Vector<double> State;
        public Vector<double> Obs;
        private bool doSave;
        private double h_state;
        private double h_obs;
        private double h_tolerance;
        private double t_nextobservation;
        private int n; // observation counter


        public Dictionary<double, Vector<double>> Trajectory;
        public Dictionary<double, Vector<double>> Observation;
        public DiscreteVectorModel DiscreteModel;

        NumberFormatInfo provider = new NumberFormatInfo();
       
        public ContinuousVectorModel(
                                       double h_state,
                                       double h_obs,
                                       Func<double, Vector<double>, Vector<double>> Phi1,
                                       Func<double, Vector<double>, Matrix<double>> Phi2,
                                       Func<double, Vector<double>, Vector<double>> Psi1,
                                       Func<double, Vector<double>, Matrix<double>> Psi2,
                                       Func<double, Vector<double>> W,
                                       Func<double, Vector<double>> Nu,
                                       Vector<double> X0,
                                       bool saveHistory = false)
        {
            if (h_obs < h_state)
                throw new ArgumentException("Discretization step for observations should not be lower then the same for the state");

            t = 0;
            this.h_state = h_state;
            this.h_obs = h_obs;
            h_tolerance = h_state / 2.0;
            t_nextobservation = t + h_obs;
            n = 0;
            this.Phi1 = Phi1;
            this.Phi2 = Phi2;
            this.Psi1 = Psi1;
            this.Psi2 = Psi2;
            this.W = W;
            this.Nu = Nu;

            doSave = saveHistory;

            Trajectory = new Dictionary<double, Vector<double>>();
            Observation = new Dictionary<double, Vector<double>>();
            State = X0;
            if (doSave) Trajectory.Add(t, X0);
            DiscreteModel = new DiscreteVectorModel(); // dummy discrete model to save discretized trajectory
            Obs = Psi1(t,State) + Psi2(t, State) * Nu(t);
            if (doSave) Observation.Add(t, Obs);

            provider.NumberDecimalSeparator = ".";
        }

        public Vector<double>[] Step()
        {
            if (t > 0)
            {
                State += h_state * Phi1(t, State) + Math.Sqrt(h_state) * Phi2(t, State) * W(t);
                if (doSave) Trajectory.Add(t, State);
            }
            if (Math.Abs(t - t_nextobservation) < h_tolerance)
            {
                Obs = Psi1(t, State) + Psi2(t, State) * Nu(t);
                t_nextobservation += h_obs;
                if (doSave) Observation.Add(t, Obs);

                DiscreteModel.Trajectory.Add(n, new Vector<double>[] { State, Obs });
                n++;
            }
            else
                Obs = null;
            Vector<double>[] result = new Vector<double>[] { State, Obs };
            t+= h_state;
            return result;
        }

        public void SaveTrajectory(string path)
        {
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(path))
            {
                foreach (var x in Trajectory.OrderBy(s => s.Key))
                {
                    outputfile.WriteLine(string.Format(provider, "{0} {1}", x.Key, string.Join(" ", x.Value.Select(s => s.ToString(provider)).ToArray())));
                }
                outputfile.Close();
            }
        }

        public void SaveObservations(string path)
        {
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(path))
            {
                foreach (var x in Observation.OrderBy(s => s.Key))
                {
                    outputfile.WriteLine(string.Format(provider, "{0} {1}", x.Key, string.Join(" ", x.Value.Select(s => s.ToString(provider)).ToArray())));
                }
                outputfile.Close();
            }
        }


    }
}
