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
    public class StaticVectorModel // X = W; Y = Phi(X) + Nu
    {
        private Func<Vector<double>, Vector<double>> Phi;
        private Func<Vector<double>, Vector<double>> InvPhi;
        private Func<Vector<double>> W;
        private Func< Vector<double>> Nu;

        public Vector<double> X;
        public Vector<double> Y;
        public Vector<double> Xinv;
        public Vector<double> YXinv;


        public StaticVectorModel(Func<Vector<double>, Vector<double>> phi,
                                 Func<Vector<double>, Vector<double>> invphi,
                                   Func<Vector<double>> w,
                                   Func<Vector<double>> nu)
        {
            Phi = phi;
            InvPhi = invphi;
            W = w;
            Nu = nu;

            X = W();
            Y = Phi(X) + Nu();
            Xinv = InvPhi(Y);
            YXinv = Y.Stack(Xinv);
            //YXinv = Y.Stack(Y);
        }
    }
}
