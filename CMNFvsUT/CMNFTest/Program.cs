using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using CMNF;
using CMNFTest.Properties;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using UKF;

namespace CMNFTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //int T = 25;
            //int N = 1000;
            //double mW = 0; double DW = 1;
            //double mNu = 0; double DNu = 1;
            //double mEta = 100; double DEta = 100;
            //Func<double, double> phi1 = new Func<double, double>(x => x / (1 + x * x));
            //Func<double, double> phi2 = new Func<double, double>(x => 1.0);
            //Func<double, double> psi = new Func<double, double>(x => Math.Pow(x, 3) + Math.Pow(x, 1));

            //TestEnvironment test1 = new TestEnvironment(T, N, phi1, phi2, psi, mW, mNu, mEta, DW, DNu, DEta, x => phi1(x) + phi2(x) * mW, (x, y) => y - psi(x) - mNu);

            //test1.GenerateOne(Path.Combine(Settings.Default.OutputFolder, "test1_estimate.txt"));
            //test1.GenerateBundle(N, Path.Combine(Settings.Default.OutputFolder, "test1_estimateAvg.txt"));

            int T = 25;
            int N = 1000000;
            double mW = 0; double DW = 1;
            double mNu = 0; double DNu = 1;
            double mEta = 100; double DEta = 100;
            Func<double, double> phi1 = new Func<double, double> (x => 1.0 / ( Math.Sign(x) * Math.Pow(Math.Abs(x), 1.0/3.0)));//(x => Math.Abs(x) / (x * x + 0.01) + 1.0);
            Func<double, double> phi2 = new Func<double, double>(x => 1.0);
            Func<double, double> psi = new Func<double, double>(x => Math.Pow(x, 3));

            TestEnvironment test2 = new TestEnvironment(false, T, N, phi1, phi2, psi, mW, mNu, mEta, DW, DNu, DEta, x => phi1(x) + phi2(x) * mW, (x, y) => y - psi(x) - mNu);

            //test2.GenerateOne(Path.Combine(Settings.Default.OutputFolder, "test2_estimate.txt"));
            test2.GenerateBundle(1000000, Path.Combine(Settings.Default.OutputFolder, "test2_estimateAvg.txt"), false);


        }
    }

    public class TestEnvironment
    {
        private int T;

        private Func<double, double> Phi1;
        private Func<double, double> Phi2;
        private Func<double, double> Psi;
        private Func<int, double> W;
        private Func<int, double> Nu;
        private double DW;
        private double DNu;
        private Func<double> X0;
        private double X0Hat;
        private double DX0Hat;


        private Func<double, double> Xi;
        private Func<double, double, double> Zeta;

        private NumberFormatInfo provider;



        public CMNFilter CMNF;
        public UKFilter UKF;
        public UKFTest UKFtest;


        private void Initialize(bool doCalculateUKF, int t, int n,
            Func<double, double> phi1,
            Func<double, double> phi2,
            Func<double, double> psi,
            Func<int, double> w,
            Func<int, double> nu,
            Func<double> x0,
            double dw,
            double dnu,
            double x0Hat,
            double dx0Hat,
            Func<double, double> xi,
            Func<double, double, double> zeta
            )
        {
            provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            T = t;

            Phi1 = phi1;
            Phi2 = phi2;
            Psi = psi;
            W = w;
            Nu = nu;
            X0 = x0;
            DW = dw;
            DNu = dnu;
            X0Hat = x0Hat;
            DX0Hat = dx0Hat;
            Xi = xi;
            Zeta = zeta;

            DiscreteScalarModel[] models = new DiscreteScalarModel[n];
            for (int i = 0; i < n; i++)
            {
                models[i] = new DiscreteScalarModel(phi1, phi2, psi, W, Nu, X0(), true);
            }

            CMNF = new CMNFilter(xi, zeta);
            CMNF.EstimateParameters(models, X0Hat, T);

            UKF = new UKFilter(1,
                new Func<Vector<double>, Vector<double>>(x => Vector<double>.Build.Dense(1, Phi1(x[0]))),
                new Func<Vector<double>, Vector<double>>(x => Vector<double>.Build.Dense(1, Psi(x[0]))),
                Matrix<double>.Build.Dense(1, 1, DW), Matrix<double>.Build.Dense(1, 1, DNu));
            //

            //UKF.EstimateParameters(models, T, X0Hat, DX0Hat, Path.Combine(Settings.Default.OutputFolder, "test1_optimize_UKF.txt"));
            if (doCalculateUKF)
                UKF.EstimateParametersRandom(models, T, X0Hat, DX0Hat, Path.Combine(Settings.Default.OutputFolder, "optimize_UKF.txt"));



            //UKFtest = new UKFTest(1, 1);
        }

        public TestEnvironment(bool doCalculateUKF, int T, int N,
            Func<double, double> phi1,
            Func<double, double> phi2,
            Func<double, double> psi,
            Func<int, double> w,
            Func<int, double> nu,
            Func<double> x0,
            double dw,
            double dnu,
            double x0Hat,
            double dx0Hat,
            Func<double, double> xi,
            Func<double, double, double> zeta
            )
        {
            Initialize(doCalculateUKF, T, N, phi1, phi2, psi, w, nu, x0, dw, dnu, x0Hat, dx0Hat, xi, zeta);
        }


        public TestEnvironment(bool doCalculateUKF, int T, int N,
           Func<double, double> phi1,
           Func<double, double> phi2,
           Func<double, double> psi,
           double mW,
           double mNu,
           double mEta,
           double DW,
           double DNu,
           double DEta,
           Func<double, double> xi,
           Func<double, double, double> zeta
           )
        {

            Normal NormalW = new Normal(mW, Math.Sqrt(DW));
            Normal NormalNu = new Normal(mNu, Math.Sqrt(DNu));
            Normal NormalEta = new Normal(mEta, Math.Sqrt(DEta));

            Initialize(doCalculateUKF, T, N, phi1, phi2, psi, (x) => NormalW.Sample(), (x) => NormalNu.Sample(), () => NormalEta.Sample(), DW, DNu, mEta, DEta, xi, zeta);

        }

        public void GenerateOne(string fileName)
        {
            DiscreteScalarModel model = new DiscreteScalarModel(Phi1, Phi2, Psi, W, Nu, X0(), true);
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
            {
                double xHat = X0Hat;
                Vector<double> xHatU = Vector<double>.Build.Dense(1, X0Hat);
                Matrix<double> PHatU = Matrix<double>.Build.Dense(1, 1, DX0Hat);

                Matrix<double> xHatU_ = Matrix<double>.Build.Dense(1, 1, X0Hat);
                Matrix<double> PHatU_ = Matrix<double>.Build.Dense(1, 1, DX0Hat);


                //Vector<double> xHatU = Vector<double>.Build.Dense(1, X0Hat);
                //Matrix<double> PHatU = Matrix<double>.Build.Dense(1, 1, DX0Hat);
                //Matrix<double> t1 = UKF.GenerateSigmaPoints(xHatU, PHatU);
                //Matrix<double> t2 = UKFtest.GetSigmaPoints(xHatU.ToColumnMatrix(), PHatU, UKFtest.c);

                //Matrix<double> ups;
                //Vector<double> y;
                //Matrix<double> Py;
                //UKF.UT(UKF.Phi, t1, UKF.Rw, out ups, out y, out Py);
                //Matrix<double>[] o = UKFtest.UnscentedTransform(
                //    new Func<Matrix<double>, Matrix<double>>(x => Matrix<double>.Build.Dense(1, 1, Phi1(x.At(0, 0)))),
                //    t2, UKFtest.Wm, UKFtest.Wc, 1, UKF.Rw);


                for (int t = 0; t < T; t++)
                {
                    double y = model.Step();
                    xHat = CMNF.Step(t, y, xHat);
                    UKF.Step(Vector<double>.Build.Dense(1, y), xHatU, PHatU, out xHatU, out PHatU);

                    //Matrix<double>[] o = UKFtest.Update(
                    //    new Func<Matrix<double>, Matrix<double>>(x => Matrix<double>.Build.Dense(1, 1, Phi1(x[0, 0]))),
                    //    xHatU_, PHatU_,
                    //    new Func<Matrix<double>, Matrix<double>>(x => Matrix<double>.Build.Dense(1, 1, Psi(x[0, 0]))),
                    //    Matrix<double>.Build.Dense(1, 1, y),
                    //    Matrix<double>.Build.Dense(1, 1, DW),
                    //    Matrix<double>.Build.Dense(1, 1, DNu)
                    //    );
                    //xHatU_ = o[0];
                    //PHatU_ = o[1];
                    outputfile.WriteLine(string.Format(provider, "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}",
                        t, model.State, xHat, Math.Abs(model.State - xHat), Math.Sqrt(CMNF.KHat[t]),
                        xHatU[0], Math.Abs(model.State - xHatU[0]), Math.Sqrt(PHatU[0, 0]),
                        0, 0, 0
                        //xHatU_[0,0], Math.Abs(model.State - xHatU_[0, 0]), Math.Sqrt(PHatU_[0, 0])
                        ));
                }
                outputfile.Close();
            }
        }
        public void GenerateBundle(int N, string fileName, bool doCalculateUKF)
        {
            DiscreteScalarModel[] modelsEst = new DiscreteScalarModel[N];
            for (int i = 0; i < N; i++)
            {
                modelsEst[i] = new DiscreteScalarModel(Phi1, Phi2, Psi, W, Nu, X0(), true);
            }
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
            {
                Vector<double> xHat = Vector<double>.Build.Dense(N, X0Hat);

                Vector<double> xHatU = Vector<double>.Build.Dense(1, 0);// = Vector<double>.Build.Dense(N, X0Hat);
                Vector<double> PHatU = Vector<double>.Build.Dense(1, 0);// = Vector<double>.Build.Dense(N, DX0Hat);

                Vector<double> xHatU_ = Vector<double>.Build.Dense(1, 0);// = Vector<double>.Build.Dense(N, X0Hat);
                Vector<double> PHatU_ = Vector<double>.Build.Dense(1, 0);// = Vector<double>.Build.Dense(N, DX0Hat);

                if (doCalculateUKF)
                {
                    xHatU = Vector<double>.Build.Dense(N, X0Hat);
                    PHatU = Vector<double>.Build.Dense(N, DX0Hat);

                    xHatU_ = Vector<double>.Build.Dense(N, X0Hat);
                    PHatU_ = Vector<double>.Build.Dense(N, DX0Hat);

                }

                for (int t = 0; t < T; t++)
                {
                    Vector<double> y = Vector<double>.Build.Dense(N, i => modelsEst[i].Step());
                    Vector<double> x = Vector<double>.Build.Dense(N, i => modelsEst[i].State);
                    xHat = Vector<double>.Build.Dense(N, i => CMNF.Step(t, y[i], xHat[i]));

                    double mx = x.Average();
                    double Dx = ((x - mx).PointwiseMultiply(x - mx)).Average();
                    Vector<double> error = (x - xHat);
                    Vector<double> errorPow2 = (x - xHat).PointwiseMultiply(x - xHat);
                    outputfile.Write(string.Format(provider, "{0} {1} {2} {3} {4} {5} {6}",
                        t, mx, Dx, xHat.Average(), error.Average(), errorPow2.Average(), CMNF.KHat[t]
                        ));

                    if (doCalculateUKF)
                    {
                        for (int i = 0; i < N; i++)
                        {
                            Vector<double> xHatU_i;
                            Matrix<double> PHatU_i;
                            UKF.Step(Vector<double>.Build.Dense(1, y[i]),
                                Vector<double>.Build.Dense(1, xHatU[i]),
                                Matrix<double>.Build.Dense(1, 1, PHatU[i]),
                                out xHatU_i, out PHatU_i);
                            xHatU[i] = xHatU_i[0];
                            PHatU[i] = PHatU_i[0, 0];

                            ////Matrix<double>[] o = UKFtest.Update(
                            ////    new Func<Matrix<double>, Matrix<double>>(a => Matrix<double>.Build.Dense(1, 1, Phi1(a[0, 0]))),
                            ////    Matrix<double>.Build.Dense(1, 1, xHatU_[i]),
                            ////    Matrix<double>.Build.Dense(1, 1, PHatU_[i]),
                            ////    new Func<Matrix<double>, Matrix<double>>(a => Matrix<double>.Build.Dense(1, 1, Psi(a[0, 0]))),
                            ////    Matrix<double>.Build.Dense(1, 1, y[i]),
                            ////    Matrix<double>.Build.Dense(1, 1, DW),
                            ////    Matrix<double>.Build.Dense(1, 1, DNu)
                            ////    );
                            //xHatU_[i] = o[0][0,0];
                            //PHatU_[i] = o[1][0,0];

                        }
                        Vector<double> errorU = (x - xHatU);
                        Vector<double> errorU_ = (x - xHatU_);
                        Vector<double> errorUPow2 = (x - xHatU).PointwiseMultiply(x - xHatU);
                        Vector<double> errorU_Pow2 = (x - xHatU_).PointwiseMultiply(x - xHatU_);
                        outputfile.Write(string.Format(provider, " {0} {1} {2} {3} {4} {5} {6} {7}",
                            xHatU.Average(), errorU.Average(), errorUPow2.Average(), PHatU.Average(),
                            xHatU_.Average(), errorU_.Average(), errorU_Pow2.Average(), PHatU_.Average()
                            ));
                    }
                    else
                    {
                        outputfile.Write(string.Format(provider, " {0} {1} {2} {3} {4} {5} {6} {7}", 0, 0, 0, 0, 0, 0, 0, 0));
                    }
                    outputfile.WriteLine();


                }
                outputfile.Close();
            }




        }
    }
}
