using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CMNF;

namespace CMNFTest
{
    public class CMNFWrapper: BasicFilter
    {
        private CMNFilter CMNF;
        public int T;
        public Vector<double> X0Hat;
        public DiscreteVectorModel[] Models;
        public Func<int, Vector<double>, Vector<double>> Xi;
        public Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> Zeta;

        public override void Initialize()
        {
            CMNF = new CMNFilter(Xi, Zeta);
            CMNF.EstimateParameters(Models, X0Hat, T);
        }

        public override (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat)
        {
            return CMNF.Step(t, y, xHat);
        }
    }
}
