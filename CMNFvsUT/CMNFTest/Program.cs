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

namespace CMNFTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //int T = 100;
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

            int T = 100;
            int N = 1000;
            double mW = 0; double DW = 10;
            double mNu = 0; double DNu = 1;
            double mEta = 100; double DEta = 100;
            Func<double, double> phi1 = new Func<double, double>(x => Math.Abs(x) / (x * x) + 1);
            Func<double, double> phi2 = new Func<double, double>(x => 1.0);
            Func<double, double> psi = new Func<double, double>(x => x);

            TestEnvironment test1 = new TestEnvironment(T, N, phi1, phi2, psi, mW, mNu, mEta, DW, DNu, DEta, x => phi1(x) + phi2(x) * mW, (x, y) => y - psi(x) - mNu);

            test1.GenerateOne(Path.Combine(Settings.Default.OutputFolder, "test2_estimate.txt"));
            test1.GenerateBundle(N, Path.Combine(Settings.Default.OutputFolder, "test2_estimateAvg.txt"));


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
        private Func<double> X0;
        private double X0Hat;


        private Func<double, double> Xi;
        private Func<double, double, double> Zeta;

        private NumberFormatInfo provider;



        public CMNFilter CMNF;


        private void Initialize(int t, int n,
            Func<double, double> phi1,
            Func<double, double> phi2,
            Func<double, double> psi,
            Func<int, double> w,
            Func<int, double> nu,
            Func<double> x0,
            double x0Hat,
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
            X0Hat = x0Hat;
            Xi = xi;
            Zeta = zeta;        

            DiscreteScalarModel[] models = new DiscreteScalarModel[n];
            for (int i = 0; i < n; i++)
            {
                models[i] = new DiscreteScalarModel(phi1, phi2, psi, W, Nu, X0(), true);
            }

            CMNF = new CMNFilter(xi, zeta);
            CMNF.EstimateParameters(models, X0Hat, T);
        }

        public TestEnvironment(int T, int N,
            Func<double, double> phi1,
            Func<double, double> phi2,
            Func<double, double> psi,
            Func<int, double> w,
            Func<int, double> nu,
            Func<double> x0,
            double x0Hat,
            Func<double, double> xi,
            Func<double, double, double> zeta
            )
        {
            Initialize(T, N, phi1, phi2, psi, w, nu, x0, x0Hat, xi, zeta);
        }


        public TestEnvironment(int T, int N,
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

            Normal NormalW = new Normal(mW, DW);
            Normal NormalNu = new Normal(mNu, DNu);
            Normal NormalEta = new Normal(mEta, DEta);

            Initialize(T, N, phi1, phi2, psi, (x) => NormalW.Sample(), (x) => NormalNu.Sample(), () => NormalEta.Sample(), mEta, xi, zeta);

        }

        public void GenerateOne(string fileName)
        {
            DiscreteScalarModel model = new DiscreteScalarModel(Phi1, Phi2, Psi, W, Nu, X0(), true);
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
            {
                double xHat = X0Hat;
                for (int t = 0; t < T; t++)
                {
                    double y = model.Step();
                    xHat = CMNF.Step(t, y, xHat);
                    outputfile.WriteLine(string.Format(provider, "{0} {1} {2} {3} {4}", t, model.State, xHat, Math.Abs(model.State - xHat), CMNF.KHat[t]));

                }
                outputfile.Close();
            }
        }

        public void GenerateBundle(int N, string fileName)
        {
            DiscreteScalarModel[] modelsEst = new DiscreteScalarModel[N];
            for (int i = 0; i < N; i++)
            {
                modelsEst[i] = new DiscreteScalarModel(Phi1, Phi2, Psi, W, Nu, X0(), true);
            }

            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
            {
                Vector<double> xHat = Vector<double>.Build.Dense(N, X0Hat);
                for (int t = 0; t < T; t++)
                {
                    Vector<double> y = Vector<double>.Build.Dense(N, i => modelsEst[i].Step());
                    Vector<double> x = Vector<double>.Build.Dense(N, i => modelsEst[i].State);
                    xHat = Vector<double>.Build.Dense(N, i => CMNF.Step(t, y[i], xHat[i]));
                    Vector<double> error = (x - xHat).PointwiseAbs();
                    outputfile.WriteLine(string.Format(provider, "{0} {1} {2} {3} {4}", t, x.Average(), xHat.Average(), error.Average(), CMNF.KHat[t]));
                }
                outputfile.Close();
            }


        }
    }
}
