using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;
using MathNet.Numerics.LinearAlgebra;

namespace NonlinearSystem
{
    public class DiscreteVectorModel
    {
        private Func<int, Vector<double>, Vector<double>> Phi1; // Phi1(t, X)
        private Func<int, Vector<double>, Matrix<double>> Phi2;
        private Func<int, Vector<double>, Vector<double>> Psi1;
        private Func<int, Vector<double>, Matrix<double>> Psi2;
        private Func<int, Vector<double>> W;
        private Func<int, Vector<double>> Nu;
        private int t;
        public Vector<double> State;
        public Vector<double> Obs;
        private bool doSave;


      

        public Dictionary<int,Vector<double>[]> Trajectory;

        public DiscreteVectorModel() // dummy constructor
        {
            t = 0;
            Trajectory = new Dictionary<int, Vector<double>[]>();
        }


        public DiscreteVectorModel(Func<int, Vector<double>, Vector<double>> phi1, 
                                   Func<int, Vector<double>, Matrix<double>> phi2,
                                   Func<int, Vector<double>, Vector<double>> psi1,
                                   Func<int, Vector<double>, Matrix<double>> psi2,
                                   Func<int, Vector<double>> w, 
                                   Func<int, Vector<double>> nu, 
                                   Vector<double> X0,
                                   bool saveHistory = false)
        {
            t = 0;
            Phi1 = phi1;
            Phi2 = phi2;
            Psi1 = psi1;
            Psi2 = psi2;
            W = w;
            Nu = nu;

            doSave = saveHistory;

            Trajectory = new Dictionary<int, Vector<double>[]>();

            State = X0;
            //if (doSave) Trajectory.Add(t, State);
            //Obs = Psi(State) + Nu(t);
        }

        public Vector<double>[] Step()
        {
            if (t>0)
            {
                State = Phi1(t, State) + Phi2(t, State) * W(t);
            }
            Obs = Psi1(t,State) + Psi2(t, State) * Nu(t);
            Vector<double>[] result = new Vector<double>[] { State, Obs };
            if (doSave) Trajectory.Add(t, result);
            t++;
            return result;
        }

        public void SaveTrajectory(string path)
        {
            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(path))
            {
                foreach (var x in Trajectory.OrderBy(s => s.Key))
                {
                    outputfile.WriteLine(string.Format(provider, "{0} {1} {2}", x.Key, x.Value[0].ToArray(), x.Value[1].ToArray()));
                }
                outputfile.Close();
            }
        }



    }
}
