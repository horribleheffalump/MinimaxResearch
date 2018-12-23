using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CMNF;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;

namespace TestEnvironments
{
    public class MCMNFWrapper : BasicFilter
    {
        private ModifiedCMNFilter MCMNF;
        public int N;
        public int T;
        public Vector<double> X0Hat;
        public Func<int, Vector<double>, Vector<double>> Xi;
        public Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> Zeta;
        public Func<int, Vector<double>, Vector<double>> Phi1;
        public Func<int, Vector<double>, Matrix<double>> Phi2;
        public Func<int, Vector<double>, Vector<double>> Psi1;
        public Func<int, Vector<double>, Matrix<double>> Psi2;
        public Func<int, Vector<double>> W;
        public Func<int, Vector<double>> Nu;

        public override void Initialize()
        {
            MCMNF = new ModifiedCMNFilter(Xi, Zeta, Phi1, Phi2, Psi1, Psi2, W, Nu);
        }
        public override void InitializeAndTrain()
        {
            Initialize();
        }

        public override (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat)
        {
            return MCMNF.Step(t, y, xHat, kHat, N);
        }
    }
}
