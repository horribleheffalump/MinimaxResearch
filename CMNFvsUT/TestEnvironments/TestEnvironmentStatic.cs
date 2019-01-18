using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using UKF;
using MathNetExtensions;
using System.IO;

namespace TestEnvironments
{
    public class TestEnvironmentStatic
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

        private Matrix<double> Kxx;
        private Matrix<double> Kxy;
        private Matrix<double> Kyy;
        private Vector<double> My;
        public Matrix<double> P;

        private Matrix<double> Kxy_inv;
        private Matrix<double> Kyy_inv;
        private Vector<double> My_inv;
        public Matrix<double> P_inv;

        private Matrix<double> Kxy_lin;
        private Matrix<double> Kyy_lin;
        private Vector<double> My_lin;
        public Matrix<double> P_lin;

        public Vector<double> MX;
        public Matrix<double> KX;

        public Matrix<double> KNu;

        private UTStaticEstimate utStaticEstimate;

        public void Initialize(int n, string outputFolder)
        {
            provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            StaticVectorModel[] models = new StaticVectorModel[n];
            Vector<double>[] X = new Vector<double>[n];
            Vector<double>[] Y = new Vector<double>[n];
            Vector<double>[] Xinv = new Vector<double>[n];
            Vector<double>[] YXinv = new Vector<double>[n];
            for (int i = 0; i < n; i++)
            {
                models[i] = new StaticVectorModel(Phi, InvPhi, W, Nu);
                X[i] = models[i].X;
                Y[i] = models[i].Y;
                Xinv[i] = models[i].Xinv;
                YXinv[i] = models[i].YXinv;
            }

            Kxx = Exts.Cov(X, X);
            Kxy = Exts.Cov(X, YXinv);
            Kyy = Exts.Cov(YXinv, YXinv);
            My = YXinv.Average();
            P = Kxy * (Kyy.PseudoInverse());

            Kxy_inv = Exts.Cov(X, Xinv);
            Kyy_inv = Exts.Cov(Xinv, Xinv);
            My_inv = Xinv.Average();
            P_inv = Kxy_inv * (Kyy_inv.PseudoInverse());

            Kxy_lin = Exts.Cov(X, Y);
            Kyy_lin = Exts.Cov(Y, Y);
            My_lin = Y.Average();
            P_lin = Kxy_lin * (Kyy_lin.PseudoInverse());

            utStaticEstimate = new UTStaticEstimate(UTDefinitionType.ImplicitAlphaBetaKappa, OptimizationMethod.NelderMeed);
            utStaticEstimate.EstimateParameters(Phi, x => x.Trace(), X, Y, MX, KX, KNu, outputFolder);
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(Path.Combine(outputFolder, "UTStaticEstimateParams.txt")))
            {
                outputfile.WriteLine(utStaticEstimate.utParams.ToString());
                outputfile.Close();
            }
        }

        public void GenerateBundle(int n,
            out Vector<double> mErr, out Matrix<double> KErr, out Matrix<double> KErrTh,
            out Vector<double> mErr_inv, out Matrix<double> KErr_inv, out Matrix<double> KErrTh_inv,
            out Vector<double> mErr_lin, out Matrix<double> KErr_lin, out Matrix<double> KErrTh_lin,
            out Vector<double> mErr_UT, out Matrix<double> KErr_UT, out Matrix<double> KErrTh_UT,
            string fileName = null
            )
        {
            StaticVectorModel[] models = new StaticVectorModel[n];
            Vector<double>[] X = new Vector<double>[n];
            Vector<double>[] Y = new Vector<double>[n];
            Vector<double>[] Xinv = new Vector<double>[n];
            Vector<double>[] YXinv = new Vector<double>[n];
            for (int i = 0; i < n; i++)
            {
                models[i] = new StaticVectorModel(Phi, InvPhi, W, Nu);
                X[i] = models[i].X;
                Y[i] = models[i].Y;
                Xinv[i] = models[i].Xinv;
                YXinv[i] = models[i].YXinv;
            }

            Vector<double>[] Xhat = YXinv.Select(y => MX + P * (y - My)).ToArray();
            Vector<double>[] Err = Xhat.Subtract(X);
            mErr = Err.Average();
            KErr = Exts.Cov(Err, Err);
            //KErrTh = KX - P * Kyy * P.Transpose();
            KErrTh = Kxx - P * Kyy * P.Transpose();

            Vector<double>[] Xhat_inv = Xinv.Select(y => MX + P_inv * (y - My_inv)).ToArray();
            Vector<double>[] Err_inv = Xhat_inv.Subtract(X);
            mErr_inv = Err_inv.Average();
            KErr_inv = Exts.Cov(Err_inv, Err_inv);
            //KErrTh_inv = KX - P_inv * Kyy_inv * P_inv.Transpose();
            KErrTh_inv = Kxx - P_inv * Kyy_inv * P_inv.Transpose();

            Vector<double>[] Xhat_lin = Y.Select(y => MX + P_lin * (y - My_lin)).ToArray();
            Vector<double>[] Err_lin = Xhat_lin.Subtract(X);
            mErr_lin = Err_lin.Average();
            KErr_lin = Exts.Cov(Err_lin, Err_lin);
            //KErrTh_lin = KX - P_lin * Kyy_lin * P_lin.Transpose();
            KErrTh_lin = Kxx - P_lin * Kyy_lin * P_lin.Transpose();

            Vector<double>[] Xhat_UT = utStaticEstimate.Estimate(Phi, X, Y, MX, KX, KNu, out mErr_UT, out KErr_UT, out KErrTh_UT);

            if (!string.IsNullOrWhiteSpace(fileName))
            {
                string X_head = string.Join("; ", Enumerable.Range(0, X[0].Count).Select(i => $"X_{i}"));
                string Y_head = string.Join("; ", Enumerable.Range(0, Y[0].Count).Select(i => $"Y_{i}"));
                string Xinv_head = string.Join("; ", Enumerable.Range(0, Xinv[0].Count).Select(i => $"Xinv_{i}"));
                string Xhat_head = string.Join("; ", Enumerable.Range(0, Xhat[0].Count).Select(i => $"Xhat_{i}"));
                string Xhat_inv_head = string.Join("; ", Enumerable.Range(0, Xhat_inv[0].Count).Select(i => $"Xhat_inv_{i}"));
                string Xhat_lin_head = string.Join("; ", Enumerable.Range(0, Xhat_lin[0].Count).Select(i => $"Xhat_lin_{i}"));
                string Xhat_UT_head = string.Join("; ", Enumerable.Range(0, Xhat_lin[0].Count).Select(i => $"Xhat_UT_{i}"));

                using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
                {
                    outputfile.WriteLine($"{X_head}; {Y_head}; {Xinv_head}; {Xhat_head}; {Xhat_inv_head}; {Xhat_lin_head}; {Xhat_UT_head}");
                    for (int i = 0; i < n; i++)
                    {
                        outputfile.WriteLine($"{X[i].ToLine()}; {Y[i].ToLine()}; {Xinv[i].ToLine()}; {Xhat[i].ToLine()}; {Xhat_inv[i].ToLine()}; {Xhat_lin[i].ToLine()}; {Xhat_UT[i].ToLine()}");
                    }
                    outputfile.Close();
                }
            }
        }


    }
}
