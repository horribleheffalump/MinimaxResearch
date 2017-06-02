using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using UKF;

namespace CMNFTest
{
    class TestEnvironmentStatic
    {
        //public StaticVectorModel[] models;

        //public Vector<double>[] X; // = Enumerable.Repeat(0.0, N).Select(i => Vector(NormalX[0].Sample(), NormalX[1].Sample())).ToArray();
        //public Vector<double>[] Y; // = X.Select(x => Extensions.cart2pol(x) + Vector(NormalNu[0].Sample(), NormalNu[1].Sample())).ToArray();
        //public Vector<double>[] Xinv; // = Y.Select(y => Extensions.pol2cart(y)).ToArray();
        //public Vector<double>[] YXinv;// = Y.Stack(Xinv);

        public Func<Vector<double>, Vector<double>> Phi;
        public Func<Vector<double>, Vector<double>> InvPhi;
        public Func<Vector<double>> W;
        public Func<Vector<double>> Nu;

        private NumberFormatInfo provider;

        private Matrix<double> Kxy;
        private Matrix<double> Kyy;
        private Vector<double> My;
        public Matrix<double> P;

        private Matrix<double> Kxy_lin;
        private Matrix<double> Kyy_lin;
        private Vector<double> My_lin;
        public Matrix<double> P_lin;

        public Vector<double> MX;
        public Matrix<double> KX;

        public Matrix<double> KNu;

        private UTStaticEstimate utStaticEstimate;

        public UTOptimizationType utOptimizationType;

        public void Initialize(int n)
        {
            provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            StaticVectorModel[] models = new StaticVectorModel[n];
            Vector<double>[] X = new Vector<double>[n];
            Vector<double>[] Y = new Vector<double>[n];
            //Vector<double>[] Xinv = new Vector<double>[n];
            Vector<double>[] YXinv = new Vector<double>[n];
            for (int i = 0; i < n; i++)
            {
                models[i] = new StaticVectorModel(Phi, InvPhi, W, Nu);
                X[i] = models[i].X;
                Y[i] = models[i].Y;
                //Xinv[i] = models[i].Xinv;
                YXinv[i] = models[i].YXinv;
            }

            Kxy = Extensions.Cov(X, YXinv);
            Kyy = Extensions.Cov(YXinv, YXinv);
            My = YXinv.Average();
            P = Kxy * (Kyy.PseudoInverse());

            Kxy_lin = Extensions.Cov(X, Y);
            Kyy_lin = Extensions.Cov(Y, Y);
            My_lin = Y.Average();
            P_lin = Kxy_lin * (Kyy_lin.PseudoInverse());

            utStaticEstimate = new UTStaticEstimate(utOptimizationType);
            utStaticEstimate.EstimateParameters(100000, 10000, Phi, x => x.Trace(), X, Y, MX, KX, KNu);
        }

        public void GenerateBundle(int n, 
            out Vector<double> mErr, out Matrix<double> KErr, out Matrix<double> KErrTh, 
            out Vector<double> mErr_lin, out Matrix<double> KErr_lin, out Matrix<double> KErrTh_lin,
            out Vector<double> mErr_UT, out Matrix<double> KErr_UT, out Matrix<double> KErrTh_UT
            )
        {
            StaticVectorModel[] models = new StaticVectorModel[n];
            Vector<double>[] X = new Vector<double>[n];
            Vector<double>[] Y = new Vector<double>[n];
            //Vector<double>[] Xinv = new Vector<double>[n];
            Vector<double>[] YXinv = new Vector<double>[n];
            for (int i = 0; i < n; i++)
            {
                models[i] = new StaticVectorModel(Phi, InvPhi, W, Nu);
                X[i] = models[i].X;
                Y[i] = models[i].Y;
                //Xinv[i] = models[i].Xinv;
                YXinv[i] = models[i].YXinv;
            }

            Vector<double>[] Xhat = YXinv.Select(y => MX + P * (y - My)).ToArray();
            Vector<double>[] Err = Xhat.Subtract(X);
            mErr = Err.Average();
            KErr = Extensions.Cov(Err, Err);
            KErrTh = KX - P * Kyy * P.Transpose();

            Vector<double>[] Xhat_lin = Y.Select(y => MX + P_lin * (y - My_lin)).ToArray();
            Vector<double>[] Err_lin = Xhat_lin.Subtract(X);
            mErr_lin = Err_lin.Average();
            KErr_lin = Extensions.Cov(Err_lin, Err_lin);
            KErrTh_lin = KX - P_lin * Kyy_lin * P_lin.Transpose();

            utStaticEstimate.Estimate(Phi, X, Y, MX, KX, KNu, out mErr_UT, out KErr_UT, out KErrTh_UT);

        }


    }
}
