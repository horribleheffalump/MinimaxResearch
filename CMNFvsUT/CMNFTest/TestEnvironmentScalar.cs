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
using MathNetExtensions;

namespace CMNFTest
{
    public class TestEnvironmentScalar
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



        public CMNScalarFilter CMNFScalar;
        public CMNFilter CMNF;
        public UKFilter UKF;


        private void Initialize(int N1, int N2, bool doCalculateUKF, int t, int n,
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

            //DiscreteScalarModel[] models = new DiscreteScalarModel[n];
            //for (int i = 0; i < n; i++)
            //{
            //    models[i] = new DiscreteScalarModel(phi1, phi2, psi, W, Nu, X0(), true);
            //}

            DiscreteVectorModel[] models = new DiscreteVectorModel[n];
            for (int i = 0; i < n; i++)
            {
                models[i] = new DiscreteVectorModel
                    (
                    new Func<int, Vector<double>, Vector<double>>((s, x) => Exts.Vector(Phi1(x[0]))),
                    new Func<int, Vector<double>, Matrix<double>>((s, x) => Exts.Matrix(Phi2(x[0]))),
                    new Func<int, Vector<double>, Vector<double>>((s, x) => Exts.Vector(Psi(x[0]))),
                    new Func<int, Vector<double>, Matrix<double>>((s, x) => Exts.Matrix(1.0)),
                    new Func<int, Vector<double>>((s) => Exts.Vector(W(s))),
                    new Func<int, Vector<double>>((s) => Exts.Vector(Nu(s))),
                    Exts.Vector(X0()),
                    true);
            }
            for (int s = 0; s < T; s++)
            {
                for (int i = 0; i < n; i++)
                {
                    models[i].Step();
                }
            }

            CMNFScalar = new CMNScalarFilter(xi, zeta);
            CMNFScalar.EstimateParameters(models, Exts.Vector(X0Hat), T);

            CMNF = new CMNFilter(
                new Func<int, Vector<double>, Vector<double>>((s, x) => Exts.Vector(xi(x[0]))),
                new Func<int, Vector<double>, Vector<double>, Vector<double>>((s, x, y) => Exts.Vector(zeta(x[0], y[0])))
                );
            CMNF.EstimateParameters(models, Exts.Vector(X0Hat), T);

            //Console.WriteLine("f, F, h, H, K");
            //for (int i = 0; i < T; i++)
            //{
            //    Console.WriteLine($"{Math.Abs(CMNFScalar.fHat[i] - CMNF.fHat[i][0])}, {Math.Abs(CMNFScalar.FHat[i] - CMNF.FHat[i][0,0])}, {Math.Abs(CMNFScalar.hHat[i] - CMNF.hHat[i][0])}, {Math.Abs(CMNFScalar.HHat[i] - CMNF.HHat[i][0, 0])}, {Math.Abs(CMNFScalar.KHat[i] - CMNF.KHat[i][0, 0])}");
            //}

            UKF = new UKFilter();
                //1,
                //new Func<int, Vector<double>, Vector<double>>((s,x) => Exts.Vector(Phi1(x[0]))),
                //new Func<int, Vector<double>, Vector<double>>((s,x) => Exts.Vector(Psi(x[0]))),
                //Exts.Matrix(DW), Exts.Matrix(DNu));
                //

            //UKF.EstimateParameters(models, T, X0Hat, DX0Hat, Path.Combine(Settings.Default.OutputFolder, "test1_optimize_UKF.txt"));
            if (doCalculateUKF)
                UKF.EstimateParameters(N1, N2,
                    (s, x) => Exts.Vector(Phi1(x[0])),
                    (s, x) => Exts.Vector(Psi(x[0])),
                    Exts.Matrix(DW), Exts.Matrix(DNu),
                    x => x.Trace(), T,  
                    models, Exts.Vector(X0Hat), Exts.Matrix(DX0Hat), Settings.Default.OutputFolder);
            //UKFtest = new UKFTest(1, 1);
        }

        public TestEnvironmentScalar(int N1, int N2, bool doCalculateUKF, int T, int N,
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
            Initialize(N1, N2, doCalculateUKF, T, N, phi1, phi2, psi, w, nu, x0, dw, dnu, x0Hat, dx0Hat, xi, zeta);
        }


        public TestEnvironmentScalar(int N1, int N2, bool doCalculateUKF, int T, int N,
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

            Initialize(N1, N2, doCalculateUKF, T, N, phi1, phi2, psi, (x) => NormalW.Sample(), (x) => NormalNu.Sample(), () => NormalEta.Sample(), DW, DNu, mEta, DEta, xi, zeta);

        }

        public void GenerateOne(string fileName)
        {
            DiscreteScalarModel model = new DiscreteScalarModel(Phi1, Phi2, Psi, W, Nu, X0(), true);
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
            {
                double xHat = X0Hat;
                Vector<double> xHatU = Exts.Vector(X0Hat);
                Matrix<double> PHatU = Exts.Matrix(DX0Hat);

                Matrix<double> xHatU_ = Exts.Matrix(X0Hat);
                Matrix<double> PHatU_ = Exts.Matrix(DX0Hat);


                //Vector<double> xHatU = Exts.Vector(X0Hat);
                //Matrix<double> PHatU = Exts.Matrix(DX0Hat);
                //Matrix<double> t1 = UKF.GenerateSigmaPoints(xHatU, PHatU);
                //Matrix<double> t2 = UKFtest.GetSigmaPoints(xHatU.ToColumnMatrix(), PHatU, UKFtest.c);

                //Matrix<double> ups;
                //Vector<double> y;
                //Matrix<double> Py;
                //UKF.UT(UKF.Phi, t1, UKF.Rw, out ups, out y, out Py);
                //Matrix<double>[] o = UKFtest.UnscentedTransform(
                //    new Func<Matrix<double>, Matrix<double>>(x => Exts.Matrix(Phi1(x.At(0, 0)))),
                //    t2, UKFtest.Wm, UKFtest.Wc, 1, UKF.Rw);


                for (int t = 0; t < T; t++)
                {
                    double y = model.Step();
                    xHat = CMNFScalar.Step(t, y, xHat);
                    UKF.Step((s, x) => Exts.Vector(Phi1(x[0])), (s, x) => Exts.Vector(Psi(x[0])), Exts.Matrix(DW), Exts.Matrix(DNu), t, Exts.Vector(y), xHatU, PHatU, out xHatU, out PHatU);

                    //Matrix<double>[] o = UKFtest.Update(
                    //    new Func<Matrix<double>, Matrix<double>>(x => Exts.Matrix(Phi1(x[0, 0]))),
                    //    xHatU_, PHatU_,
                    //    new Func<Matrix<double>, Matrix<double>>(x => Exts.Matrix(Psi(x[0, 0]))),
                    //    Exts.Matrix(y),
                    //    Exts.Matrix(DW),
                    //    Exts.Matrix(DNu)
                    //    );
                    //xHatU_ = o[0];
                    //PHatU_ = o[1];
                    outputfile.WriteLine(string.Format(provider, "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}",
                        t, model.State, xHat, Math.Abs(model.State - xHat), Math.Sqrt(CMNFScalar.KHat[t]),
                        xHatU[0], Math.Abs(model.State - xHatU[0]), Math.Sqrt(PHatU[0, 0]),
                        0, 0, 0
                        //xHatU_[0,0], Math.Abs(model.State - xHatU_[0, 0]), Math.Sqrt(PHatU_[0, 0])
                        ));
                }
                outputfile.Close();
            }
        }
        //public void GenerateBundle(int N, string fileName, bool doCalculateUKF)
        //{
        //    //DiscreteScalarModel[] modelsEst = new DiscreteScalarModel[N];
        //    DiscreteVectorModel[] modelsEst = new DiscreteVectorModel[N];

        //    for (int i = 0; i < N; i++)
        //    {
        //        modelsEst[i] = new DiscreteVectorModel
        //           (
        //            new Func<int, Vector<double>, Vector<double>>((s, x) => Exts.Vector(Phi1(x[0]))),
        //            new Func<int, Vector<double>, Matrix<double>>((s, x) => Exts.Matrix(Phi2(x[0]))),
        //            new Func<int, Vector<double>, Vector<double>>((s, x) => Exts.Vector(Psi(x[0]))),
        //            new Func<int, Vector<double>, Matrix<double>>((s, x) => Exts.Matrix(1.0)),
        //            new Func<int, Vector<double>>((s) => Exts.Vector(W(s))),
        //            new Func<int, Vector<double>>((s) => Exts.Vector(Nu(s))),
        //            Exts.Vector(X0()),
        //            true
        //           );

        //        for (int s = 0; s < T; s++)
        //        {
        //            modelsEst[i].Step();
        //        }
        //    }

        //    using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
        //    {
        //        Vector<double> xHat = Vector<double>.Build.Dense(N, X0Hat);

        //        Vector<double> xHatU = Exts.Vector(0);// = Vector<double>.Build.Dense(N, X0Hat);
        //        Vector<double> PHatU = Exts.Vector(0);// = Vector<double>.Build.Dense(N, DX0Hat);

        //        Vector<double> xHatU_ = Exts.Vector(0);// = Vector<double>.Build.Dense(N, X0Hat);
        //        Vector<double> PHatU_ = Exts.Vector(0);// = Vector<double>.Build.Dense(N, DX0Hat);

        //        if (doCalculateUKF)
        //        {
        //            xHatU = Vector<double>.Build.Dense(N, X0Hat);
        //            PHatU = Vector<double>.Build.Dense(N, DX0Hat);

        //            xHatU_ = Vector<double>.Build.Dense(N, X0Hat);
        //            PHatU_ = Vector<double>.Build.Dense(N, DX0Hat);

        //        }

        //        for (int t = 0; t < T; t++)
        //        {
        //            Vector<double> x = Vector<double>.Build.Dense(N, i => modelsEst[i].Trajectory[t][0][0]);
        //            Vector<double> y = Vector<double>.Build.Dense(N, i => modelsEst[i].Trajectory[t][1][0]);
        //            //xHat = Vector<double>.Build.Dense(N, i => CMNFScalar.Step(t, y[i], xHat[i]));
        //            xHat = Vector<double>.Build.Dense(N, i => CMNF.Step(t, Vector<double>.Build.Dense(1,y[i]), Exts.Vector(xHat[i]))[0]);

        //            double mx = x.Average();
        //            double Dx = ((x - mx).PointwiseMultiply(x - mx)).Average();
        //            Vector<double> error = (x - xHat);
        //            Vector<double> errorPow2 = (x - xHat).PointwiseMultiply(x - xHat);
        //            outputfile.Write(string.Format(provider, "{0} {1} {2} {3} {4} {5} {6}",
        //                t, mx, Dx, xHat.Average(), error.Average(), errorPow2.Average(), CMNFScalar.KHat[t]
        //                ));

        //            if (doCalculateUKF)
        //            {
        //                for (int i = 0; i < N; i++)
        //                {
        //                    Vector<double> xHatU_i;
        //                    Matrix<double> PHatU_i;
        //                    UKF.Step(Exts.Vector(y[i]),
        //                        Exts.Vector(xHatU[i]),
        //                        Exts.Matrix(PHatU[i]),
        //                        out xHatU_i, out PHatU_i);
        //                    xHatU[i] = xHatU_i[0];
        //                    PHatU[i] = PHatU_i[0, 0];

        //                    ////Matrix<double>[] o = UKFtest.Update(
        //                    ////    new Func<Matrix<double>, Matrix<double>>(a => Exts.Matrix(Phi1(a[0, 0]))),
        //                    ////    Exts.Matrix(xHatU_[i]),
        //                    ////    Exts.Matrix(PHatU_[i]),
        //                    ////    new Func<Matrix<double>, Matrix<double>>(a => Exts.Matrix(Psi(a[0, 0]))),
        //                    ////    Exts.Matrix(y[i]),
        //                    ////    Exts.Matrix(DW),
        //                    ////    Exts.Matrix(DNu)
        //                    ////    );
        //                    //xHatU_[i] = o[0][0,0];
        //                    //PHatU_[i] = o[1][0,0];

        //                }
        //                Vector<double> errorU = (x - xHatU);
        //                Vector<double> errorU_ = (x - xHatU_);
        //                Vector<double> errorUPow2 = (x - xHatU).PointwiseMultiply(x - xHatU);
        //                Vector<double> errorU_Pow2 = (x - xHatU_).PointwiseMultiply(x - xHatU_);
        //                outputfile.Write(string.Format(provider, " {0} {1} {2} {3} {4} {5} {6} {7}",
        //                    xHatU.Average(), errorU.Average(), errorUPow2.Average(), PHatU.Average(),
        //                    xHatU_.Average(), errorU_.Average(), errorU_Pow2.Average(), PHatU_.Average()
        //                    ));
        //            }
        //            else
        //            {
        //                outputfile.Write(string.Format(provider, " {0} {1} {2} {3} {4} {5} {6} {7}", 0, 0, 0, 0, 0, 0, 0, 0));
        //            }
        //            outputfile.WriteLine();


        //        }
        //        outputfile.Close();
        //    }




        //}

        public void GenerateBundle(int N, string fileName, bool doCalculateUKF)
        {
            //DiscreteScalarModel[] modelsEst = new DiscreteScalarModel[N];
            DiscreteVectorModel[] modelsEst = new DiscreteVectorModel[N];

            for (int i = 0; i < N; i++)
            {
                modelsEst[i] = new DiscreteVectorModel
                   (
                    new Func<int, Vector<double>, Vector<double>>((s, x) => Exts.Vector(Phi1(x[0]))),
                    new Func<int, Vector<double>, Matrix<double>>((s, x) => Exts.Matrix(Phi2(x[0]))),
                    new Func<int, Vector<double>, Vector<double>>((s, x) => Exts.Vector(Psi(x[0]))),
                    new Func<int, Vector<double>, Matrix<double>>((s, x) => Exts.Matrix(1.0)),
                    new Func<int, Vector<double>>((s) => Exts.Vector(W(s))),
                    new Func<int, Vector<double>>((s) => Exts.Vector(Nu(s))),
                    Exts.Vector(X0()),
                    true
                   );

                for (int s = 0; s < T; s++)
                {
                    modelsEst[i].Step();
                }
            }

            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
            {
                Vector<double> xHat = Vector<double>.Build.Dense(N, X0Hat);

                Vector<double> xHatU = Exts.Vector(0);// = Vector<double>.Build.Dense(N, X0Hat);
                Vector<double> PHatU = Exts.Vector(0);// = Vector<double>.Build.Dense(N, DX0Hat);

                Vector<double> xHatU_ = Exts.Vector(0);// = Vector<double>.Build.Dense(N, X0Hat);
                Vector<double> PHatU_ = Exts.Vector(0);// = Vector<double>.Build.Dense(N, DX0Hat);

                if (doCalculateUKF)
                {
                    xHatU = Vector<double>.Build.Dense(N, X0Hat);
                    PHatU = Vector<double>.Build.Dense(N, DX0Hat);

                    xHatU_ = Vector<double>.Build.Dense(N, X0Hat);
                    PHatU_ = Vector<double>.Build.Dense(N, DX0Hat);

                }

                for (int t = 0; t < T; t++)
                {
                    Vector<double> x = Vector<double>.Build.Dense(N, i => modelsEst[i].Trajectory[t][0][0]);
                    Vector<double> y = Vector<double>.Build.Dense(N, i => modelsEst[i].Trajectory[t][1][0]);
                    //xHat = Vector<double>.Build.Dense(N, i => CMNFScalar.Step(t, y[i], xHat[i]));
                    xHat = Vector<double>.Build.Dense(N, i => CMNF.Step(t, Exts.Vector(y[i]), Exts.Vector(xHat[i]))[0]);

                    double mx = x.Average();
                    double Dx = ((x - mx).PointwiseMultiply(x - mx)).Average();
                    Vector<double> error = (x - xHat);
                    Vector<double> errorPow2 = (x - xHat).PointwiseMultiply(x - xHat);
                    outputfile.Write(string.Format(provider, "{0} {1} {2} {3} {4} {5} {6}",
                        t, mx, Dx, xHat.Average(), error.Average(), errorPow2.Average(), CMNFScalar.KHat[t]
                        ));

                    if (doCalculateUKF)
                    {
                        for (int i = 0; i < N; i++)
                        {
                            Vector<double> xHatU_i;
                            Matrix<double> PHatU_i;
                            UKF.Step(
                                (s, _x) => Exts.Vector(Phi1(_x[0])),
                                (s, _x) => Exts.Vector(Psi(_x[0])),
                                Exts.Matrix(DW), Exts.Matrix(DNu),
                                t,
                                Exts.Vector(y[i]),
                                Exts.Vector(xHatU[i]),
                                Exts.Matrix(PHatU[i]),
                                out xHatU_i, out PHatU_i);
                            xHatU[i] = xHatU_i[0];
                            PHatU[i] = PHatU_i[0, 0];

                            ////Matrix<double>[] o = UKFtest.Update(
                            ////    new Func<Matrix<double>, Matrix<double>>(a => Exts.Matrix(Phi1(a[0, 0]))),
                            ////    Exts.Matrix(xHatU_[i]),
                            ////    Exts.Matrix(PHatU_[i]),
                            ////    new Func<Matrix<double>, Matrix<double>>(a => Exts.Matrix(Psi(a[0, 0]))),
                            ////    Exts.Matrix(y[i]),
                            ////    Exts.Matrix(DW),
                            ////    Exts.Matrix(DNu)
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
