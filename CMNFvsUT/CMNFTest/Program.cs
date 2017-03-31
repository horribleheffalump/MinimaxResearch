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
            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            int T = 1000;

            double mW = 0;
            double mNu = 0;
            double mEta = 100;
            double DW = 1;
            double DNu = 1;
            double DEta = 100;

            Normal NormalW = new Normal(mW, DW);
            Normal NormalNu = new Normal(mNu, DNu);
            Normal NormalEta = new Normal(mEta, DEta);


            Func<double, double> phi1 = new Func<double, double>(x => x / (1 + x * x));
            Func<double, double> phi2 = new Func<double, double>(x => 1.0);
            Func<double, double> psi = new Func<double, double>(x => Math.Pow(x, 3) + Math.Pow(x, 1));
            Func<int, double> W = new Func<int, double>(t => NormalW.Sample());
            Func<int, double> Nu = new Func<int, double>(t => NormalNu.Sample());

            Func<double, double> xi = new Func<double, double>(x => phi1(x) + phi2(x) * mW);
            Func<double, double, double> zeta = new Func<double, double, double>((x, y) => y - psi(x) - mNu);

            //filter adjustment
            int N = 1000; // bundle size

            DiscreteScalarModel[] models = new DiscreteScalarModel[N];
            for (int i = 0; i < N; i++)
            {
                models[i] = new DiscreteScalarModel(phi1, phi2, psi, W, Nu, NormalEta.Sample(), true);
            }

            CMNFilter cmnf = new CMNFilter(xi, zeta);
            cmnf.EstimateParameters(models, mEta, T);



            //bundle
            DiscreteScalarModel[] modelsEst = new DiscreteScalarModel[N];
            for (int i = 0; i < N; i++)
            {
                modelsEst[i] = new DiscreteScalarModel(phi1, phi2, psi, W, Nu, NormalEta.Sample(), true);
            }

            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(Path.Combine(Settings.Default.OutputFolder, "estimateAverage.txt")))
            {
                Vector<double> xHat = Vector<double>.Build.Dense(N, mEta);
                for (int t = 0; t < T; t++)
                {
                    Vector<double> y = Vector<double>.Build.Dense(N, i => models[i].Step());
                    Vector<double> x = Vector<double>.Build.Dense(N, i => models[i].State);
                    xHat = Vector<double>.Build.Dense(N, i => cmnf.Step(t, y[i], xHat[i]));
                    Vector<double> error = (x - xHat).PointwiseAbs();
                    outputfile.WriteLine(string.Format(provider, "{0} {1} {2} {3} {4}", t, x.Average(), xHat.Average(), error.Average(), cmnf.KHat[t]));
                }
                outputfile.Close();
            }


            //one trajevtory
            DiscreteScalarModel model = new DiscreteScalarModel(phi1, phi2, psi, W, Nu, NormalEta.Sample(), true);
            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(Path.Combine(Settings.Default.OutputFolder, "estimate.txt")))
            {
                double xHat = mEta;
                for (int t = 0; t < T; t++)
                {
                    double y = model.Step();
                    xHat = cmnf.Step(t, y, xHat);
                    outputfile.WriteLine(string.Format(provider, "{0} {1} {2} {3} {4}", t, model.State, xHat, Math.Abs(model.State - xHat), cmnf.KHat[t]));

                }
                outputfile.Close();
            }
            //model.SaveTrajectory(Path.Combine(Settings.Default.OutputFolder, "state.txt"));


        }
    }
}
