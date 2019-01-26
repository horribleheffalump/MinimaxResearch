using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CMNF;
using System.Xml.Serialization;
using System.IO;
using MathNetExtensions;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;

namespace TestEnvironments.Filters
{
    internal class RCMNFWrapper: BasicFilter
    {
        private ResamplingCMNFilter RCMNF;
        public int N;
        public int T;
        public Func<Vector<double>> X0;
        public Vector<double> X0Hat;
        public Func<int, Vector<double>, Vector<double>> Xi;
        public Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> Zeta;
        public Func<int, Vector<double>, Vector<double>> Phi1;
        public Func<int, Vector<double>, Matrix<double>> Phi2;
        public Func<int, Vector<double>, Vector<double>> Psi1;
        public Func<int, Vector<double>, Matrix<double>> Psi2;
        public Func<int, Vector<double>> W;
        public Func<int, Vector<double>> Nu;
        public Matrix<double> DNu;

        public override void Initialize()
        {
            FilterName = "RCMNF";
            RCMNF = new ResamplingCMNFilter(Xi, Zeta, Phi1, Phi2, Psi1, Psi2, W, Nu, DNu);
            RCMNF.Initialize(N, X0, X0Hat);
        }

        public override void InitializeAndTrain()
        {
            Initialize();
        }

        public override (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat)
        {
            return RCMNF.Step(t, y, xHat, kHat);
        }

    }
}
