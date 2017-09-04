﻿
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
using PythonInteract;

namespace CMNFTest
{
    /// <summary>
    /// Test environment for CMN and UT filters comparison on discrete stochastic dynamic systems
    /// </summary>
    public class TestEnvironmentVector
    {
        public string TestName;
        private int T;

        public Func<int, Vector<double>, Vector<double>> Phi1;
        public Func<int, Vector<double>, Matrix<double>> Phi2;
        public Func<int, Vector<double>, Vector<double>> Psi;

        public string Phi1_latex;
        public string Phi2_latex;
        public string Psi_latex;

        public Func<int, Vector<double>> W;
        public Func<int, Vector<double>> Nu;
        public Matrix<double> DW;
        public Matrix<double> DNu;
        public Func<Vector<double>> X0;
        public Vector<double> X0Hat;
        public Matrix<double> DX0Hat;

        public Func<int, Vector<double>, Vector<double>> Xi;
        public Func<int, Vector<double>, Vector<double>, Vector<double>> Zeta;

        private NumberFormatInfo provider;

        public CMNFilter CMNF;
        public UKFilter UKF;

        /// <summary>
        /// Initializes the test environment by calculating the statistics for CMN and UT filters
        /// </summary>
        /// <param name="doCalculateUKF"></param>
        /// <param name="t">time interval right margin (number of steps)</param>
        /// <param name="n">number of trajectories</param>
        public void Initialize(bool doCalculateUKF, int t, int n)
        {
            provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            T = t;

            DiscreteVectorModel[] models = new DiscreteVectorModel[n];
            for (int i = 0; i < n; i++)
            {
                models[i] = new DiscreteVectorModel(Phi1, Phi2, Psi, new Func<int, Vector<double>, Matrix<double>>((s, x) => Matrix<double>.Build.Dense(1, 1, 1.0)), W, Nu, X0(), true);
                for (int s = 0; s < T; s++)
                {
                    models[i].Step();
                }
            }


            CMNF = new CMNFilter(Xi, Zeta);
            CMNF.EstimateParameters(models, X0Hat, T);

            UKF = new UKFilter(X0().Count, Phi1, Psi, DW, DNu);

            if (doCalculateUKF)
                UKF.EstimateParametersRandom(models, T, X0Hat, DX0Hat, Path.Combine(Settings.Default.OutputFolder, "optimize_UKF.txt"));
        }

        //public TestEnvironmentVector()//bool doCalculateUKF, int T, int N)
        //{
        //    //Initialize(doCalculateUKF, T, N);
        //}


        //public void GenerateOne(string fileName)
        //{
        //    DiscreteScalarModel model = new DiscreteScalarModel(Phi1, Phi2, Psi, W, Nu, X0(), true);
        //    using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName))
        //    {
        //        double xHat = X0Hat;
        //        Vector<double> xHatU = Vector<double>.Build.Dense(1, X0Hat);
        //        Matrix<double> PHatU = Matrix<double>.Build.Dense(1, 1, DX0Hat);

        //        Matrix<double> xHatU_ = Matrix<double>.Build.Dense(1, 1, X0Hat);
        //        Matrix<double> PHatU_ = Matrix<double>.Build.Dense(1, 1, DX0Hat);


        //        for (int t = 0; t < T; t++)
        //        {
        //            double y = model.Step();
        //            xHat = CMNFScalar.Step(t, y, xHat);
        //            UKF.Step(Vector<double>.Build.Dense(1, y), xHatU, PHatU, out xHatU, out PHatU);

        //            outputfile.WriteLine(string.Format(provider, "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}",
        //                t, model.State, xHat, Math.Abs(model.State - xHat), Math.Sqrt(CMNFScalar.KHat[t]),
        //                xHatU[0], Math.Abs(model.State - xHatU[0]), Math.Sqrt(PHatU[0, 0]),
        //                0, 0, 0
        //                ));
        //        }
        //        outputfile.Close();
        //    }
        //}


        private void GenerateBundle(int n, string fileName, bool doCalculateUKF)
        {

            //DiscreteScalarModel[] modelsEst = new DiscreteScalarModel[N];
            DiscreteVectorModel[] modelsEst = new DiscreteVectorModel[n];
            int dimX = X0().Count;

            for (int i = 0; i < n; i++)
            {
                modelsEst[i] = new DiscreteVectorModel(Phi1, Phi2, Psi, new Func<int, Vector<double>, Matrix<double>>((s, x) => Matrix<double>.Build.Dense(1, 1, 1.0)), W, Nu, X0(), true);
                for (int s = 0; s < T; s++)
                {
                    modelsEst[i].Step();
                }
            }

            Vector<double>[] xHat = Enumerable.Repeat(X0Hat, n).ToArray();

            Vector<double>[] xHatU = Enumerable.Repeat(X0Hat, n).ToArray();
            Matrix<double>[] PHatU = Enumerable.Repeat(DX0Hat, n).ToArray();

            for (int k = 0; k < dimX; k++)
            {
                using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString())))
                {
                    outputfile.Close();
                }
            }

            for (int t = 0; t < T; t++)
            {
                Vector<double>[] x = new Vector<double>[n];
                Vector<double>[] y = new Vector<double>[n];

                for (int i = 0; i < n; i++)
                {
                    x[i] = modelsEst[i].Trajectory[t][0];
                    y[i] = modelsEst[i].Trajectory[t][1];
                    xHat[i] = CMNF.Step(t, y[i], xHat[i]);
                }

                Vector<double> mx = x.Average();
                Matrix<double> Dx = Extensions.Cov(x, x);

                Vector<double> mxHat = xHat.Average();

                Vector<double> mError = (x.Subtract(xHat)).Average();
                Matrix<double> DError = Extensions.Cov(x.Subtract(xHat), x.Subtract(xHat));


                Vector<double> mErrorU = Vector<double>.Build.Dense(dimX, 0);
                Matrix<double> DErrorU = Matrix<double>.Build.Dense(dimX, dimX, 0);

                Vector<double> mxHatU = Vector<double>.Build.Dense(dimX, 0);
                Matrix<double> mPHatU = Matrix<double>.Build.Dense(dimX, dimX, 0);

                if (doCalculateUKF)
                {
                    for (int i = 0; i < n; i++)
                    {
                        Vector<double> _xHatU;
                        Matrix<double> _PHatU;
                        UKF.Step(y[i], xHatU[i], PHatU[i], out _xHatU, out _PHatU);
                        xHatU[i] = _xHatU;
                        PHatU[i] = _PHatU;
                    }

                    mErrorU = (x.Subtract(xHatU)).Average();
                    DErrorU = Extensions.Cov(x.Subtract(xHatU), x.Subtract(xHatU));

                    mxHatU = xHatU.Average();
                    mPHatU = PHatU.Average();
                }



                for (int k = 0; k < dimX; k++)
                {
                    using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString()), true))
                    {
                        outputfile.Write(string.Format(provider, "{0} {1} {2} {3} {4} {5} {6}",
                            t, mx[k], Dx[k, k], mxHat[k], mError[k], DError[k, k], CMNF.KHat[t][k, k]
                            ));

                        outputfile.Write(string.Format(provider, " {0} {1} {2} {3}",
                            mxHatU[k], mErrorU[k], DErrorU[k, k], mPHatU[k,k]
                            ));
                        outputfile.WriteLine();
                        outputfile.Close();
                    }
                }
            }
        }

        /// <summary>
        /// Runs the test: generates a bundle of trajectories, applies the CMN and UT filters, calculates statistics for estimate errors and saves it to files
        /// </summary> 
        /// <param name="n">Number of trajectories</param>
        /// <param name="outputPath">Output files name template ({0} - number of state vector component)</param>
        public void Run(int n, string outputPath)
        {
            GenerateBundle(n, outputPath, true);
            //test1.GenerateBundle(n, Path.Combine(Settings.Default.OutputFolder, "test1_estimateAvg_{0}.txt"), true);
        }

        /// <summary>
        /// Runs a python script to process the calculated test data
        /// </summary>
        /// <param name="scriptName">Python script file path</param>
        /// <param name="scriptParamsTemplates">Script parameters templates array ({0} - number of state vector component)</param>
        public void ProcessResults(string scriptName, string[] scriptParamsTemplates)
        {
            for (int i = 0; i < X0Hat.Count; i++)
            {
                Python.RunScript(scriptName, scriptParamsTemplates.Select(s => s.Replace("{0}", i.ToString())).ToArray());
            }
            //Python.RunScript(Path.Combine(Settings.Default.ScriptsFolder, "estimate.py"), new string[] { Settings.Default.OutputFolder, "test1_estimateAvg_0.txt", "test1_estimateAvg_0.pdf" });
            //Python.RunScript(Path.Combine(Settings.Default.ScriptsFolder, "estimate.py"), new string[] { Settings.Default.OutputFolder, "test1_estimateAvg_1.txt", "test1_estimateAvg_1.pdf" });

        }

        public void GenerateReport(string templateFileName, string outputFileName)
        {
            Dictionary<string, string> replacements = new Dictionary<string, string>();
            replacements.Add("%Title%", TestName);
            replacements.Add("%phi1%", Phi1_latex);
            replacements.Add("%phi2%", Phi2_latex);
            replacements.Add("%psi%", Psi_latex);
            replacements.Add("%m_w%", "0");
            replacements.Add("%D_w%", "1");
            replacements.Add("%m_nu%", "0");
            replacements.Add("%D_nu%", "1");
            replacements.Add("%m_eta%", "0");
            replacements.Add("%D_eta%", "1");

            string template = File.ReadAllText(templateFileName, Encoding.Default);
            foreach (var pair in replacements)
            {
                template = template.Replace(pair.Key, pair.Value);
            }
            File.WriteAllText(outputFileName, template, Encoding.Default);
        }


    }
}
