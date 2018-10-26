
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using CMNF;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using UKF;
using PythonInteract;
using CMNFTest.Properties;
using MathNetExtensions;

namespace CMNFTest
{
    /// <summary>
    /// Test environment for CMN and UT filters comparison on discrete stochastic dynamic systems
    /// </summary>
    public class TestEnvironmentVector
    {
        public string TestName;
        public string TestFileName;

        private int T;

        public Func<int, Vector<double>, Vector<double>> Phi1;
        public Func<int, Vector<double>, Matrix<double>> Phi2;
        public Func<int, Vector<double>, Vector<double>> Psi1;
        public Func<int, Vector<double>, Matrix<double>> Psi2;

        public string[] Phi1_latex;
        public string[][] Phi2_latex;
        public string[] Psi1_latex;
        public string[][] Psi2_latex;

        public string P_W;
        public string P_Nu;
        public string P_Eta;


        public Func<int, Vector<double>> W;
        public Func<int, Vector<double>> Nu;

        public Vector<double> MW;
        public Matrix<double> DW;
        public Vector<double> MNu;
        public Matrix<double> DNu;
        public Func<Vector<double>> X0;
        public Vector<double> X0Hat;
        public Matrix<double> DX0Hat;

        public Func<int, Vector<double>, Vector<double>> Xi;
        public Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> Zeta;

        private NumberFormatInfo provider;

        public CMNFilter CMNF;
        public ModifiedCMNFilter MCMNF;
        public UKFilter UKF;

        public bool useSimpleModel = true;

        private void HandleNulls()
        {
            if (Phi2 == null)
                Phi2 = new Func<int, Vector<double>, Matrix<double>>((s, x) => Matrix<double>.Build.Dense(1, 1, 1.0));
            else
                useSimpleModel = false;

            if (Psi2 == null)
                Psi2 = new Func<int, Vector<double>, Matrix<double>>((s, x) => Matrix<double>.Build.Dense(1, 1, 1.0));
            else
                useSimpleModel = false;

            if (MW == null)
                MW = W(0) * 0.0;
            else
                useSimpleModel = false;

            if (MNu == null)
                MNu = Nu(0) * 0.0;
            else
                useSimpleModel = false;
        }

        /// <summary>
        /// Initializes the test environment by calculating the statistics for CMN and UT filters
        /// </summary>
        /// <param name="t">time interval right margin (number of steps)</param>
        /// <param name="n">number of trajectories</param>
        /// <param name="doCalculateUKF"></param>
        /// <param name="doCalculateUKFStepwise"></param>
        /// <param name="outputFolder"></param>
        public void Initialize(int t, int n, bool doCalculateUKF, bool doCalculateUKFStepwise, bool doOptimizeWithRandomShoot, string outputFolder)
        {
            Console.WriteLine("Init");
            provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            HandleNulls();

            T = t;

            DiscreteVectorModel[] models = new DiscreteVectorModel[n];
            for (int i = 0; i < n; i++)
            {
                if (i % 1000 == 0) // inform every 1000-th trajectory
                    Console.WriteLine($"model {i}");
                //Console.WriteLine($"model {i}");
                models[i] = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
                for (int s = 0; s < T; s++)
                {
                    models[i].Step();
                }
            }


            CMNF = new CMNFilter(Xi, Zeta);
            CMNF.EstimateParameters(models, X0Hat, T);

            MCMNF = new ModifiedCMNFilter(Xi, Zeta, Phi1, Phi2, Psi1, Psi2, W, Nu, DX0Hat);

            if (doOptimizeWithRandomShoot)
                UKF = new UKFilter(UTDefinitionType.ImplicitAlphaBetaKappa, OptimizationMethod.RandomShoot);
            else
                UKF = new UKFilter(UTDefinitionType.ImplicitAlphaBetaKappa, OptimizationMethod.NelderMeed);
    

            if (doCalculateUKFStepwise)
            {
                Console.WriteLine($"UKF estimate parameters stepwise");
                //UKF.EstimateParametersStepwise(Phi1, Phi2, Psi1, Psi2, MW, DW, MNu, DNu, x => x.Trace(), T, models, X0Hat, DX0Hat, outputFolder);
                UKF.EstimateParametersStepwise(Phi1, Phi2, Psi1, Psi2, MW, DW, MNu, DNu, x => x[0,0], T, models, X0Hat, DX0Hat, outputFolder);
            }
            else
            {
                Console.WriteLine($"UKF estimate parameters");
                if (doCalculateUKF)
                {
                    if (useSimpleModel)
                        //UKF.EstimateParameters(Phi1, Psi1, DW, DNu, x => x.Trace(), T, models, X0Hat, DX0Hat, outputFolder);
                        UKF.EstimateParameters(Phi1, Psi1, DW, DNu, x => x[0,0], T, models, X0Hat, DX0Hat, outputFolder);
                    else
                        //UKF.EstimateParameters(Phi1, Phi2, Psi1, Psi2, MW, DW, MNu, DNu, x => x.Trace(), T, models, X0Hat, DX0Hat, outputFolder);
                        UKF.EstimateParameters(Phi1, Phi2, Psi1, Psi2, MW, DW, MNu, DNu, x => x[0,0], T, models, X0Hat, DX0Hat, outputFolder);
                }
            }
        }

        public TestEnvironmentVector()
        {
        }


        /// <summary>
        /// Generates a bundle of trajectories and saves the state dynamics to files
        /// </summary>
        /// <param name="t">time interval right margin (number of steps)</param>
        /// <param name="n">number of trajectories</param>
        public void GenerateBundleSamples(int t, int n, string outputFolder)
        {
            HandleNulls();

            provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            T = t;

            DiscreteVectorModel[] models = new DiscreteVectorModel[n];
            for (int i = 0; i < n; i++)
            {
                if (i % 1000 == 0) // inform every 1000-th trajectory
                    Console.WriteLine($"model {i}");
                //Console.WriteLine($"model {i}");
                models[i] = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
                for (int s = 0; s < T; s++)
                {
                    models[i].Step();
                }
            }

            for (int k = 0; k < models[0].State.Count; k++)
            {
                string fileName = Path.Combine(outputFolder, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeBulk));
                using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString())))
                {
                    for (int s = 0; s < T; s++)
                    {
                        outputfile.Write(string.Format(provider, "{0} ", s));
                    }
                    outputfile.WriteLine();
                    for (int i = 0; i < n; i++)
                    {
                        for (int s = 0; s < T; s++)
                        {
                            outputfile.Write(string.Format(provider, "{0} ", models[i].Trajectory[s][0][k]));
                        }
                        outputfile.WriteLine();
                    }

                    //outputfile.Write(string.Format(provider, "{0} {1} {2} {3} {4} {5} {6}",
                    //    t, mx[k], Dx[k, k], mxHat[k], mError[k], DError[k, k], CMNF.KHat[t][k, k]
                    //    ));

                    //outputfile.Write(string.Format(provider, " {0} {1} {2} {3}",
                    //    mxHatU[k], mErrorU[k], DErrorU[k, k], mPHatU[k, k]
                    //    ));
                    outputfile.WriteLine();
                    outputfile.Close();
                }
            }


        }

        public void GenerateOne(string folderName, bool doCalculateUKF = true, int? n = null)
        {
            string fileName_state = "";
            if (n == null)
                fileName_state = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeOneState));
            else
                fileName_state = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeOneState + "_" + n.ToString()));

            string fileName_obs = "";
            if (n == null)
                fileName_obs = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeOneObs));
            else
                fileName_obs = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeOneObs + "_" + n.ToString()));

            //int dimX = X0().Count;

            DiscreteVectorModel modelEst = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
            for (int s = 0; s < T; s++)
            {
                modelEst.Step();
            }

            int dimX = modelEst.Trajectory[0][0].Count;
            int dimY = modelEst.Trajectory[0][1].Count;

            Vector<double> xHat = X0Hat;
            Vector<double> xHatU = X0Hat;
            Matrix<double> PHatU = DX0Hat;

            for (int k = 0; k < dimX; k++)
            {
                using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName_state.Replace("{0}", k.ToString())))
                {
                    outputfile.Close();
                }
            }
            for (int k = 0; k < dimY; k++)
            {
                using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName_obs.Replace("{0}", k.ToString())))
                {
                    outputfile.Close();
                }
            }

            for (int t = 0; t < T; t++)
            {
                var x = modelEst.Trajectory[t][0];
                var y = modelEst.Trajectory[t][1];
                xHat = CMNF.Step(t, y, xHat);

                Vector<double> mError = x - xHat;

                if (doCalculateUKF)
                {
                    if (useSimpleModel)
                    {
                        UKF.Step(Phi1, Psi1, DW, DNu, t, y, xHatU, PHatU, out Vector<double> _xHatU, out Matrix<double> _PHatU);
                        xHatU = _xHatU;
                        PHatU = _PHatU;
                    }
                    else
                    {
                        UKF.Step(Phi1, Phi2, Psi1, Psi2, MW, DW, MNu, DNu, t, y, xHatU, PHatU, out Vector<double> _xHatU, out Matrix<double> _PHatU);
                        xHatU = _xHatU;
                        PHatU = _PHatU;
                    }
                }
                Vector<double> mErrorU = x - xHatU;

                for (int k = 0; k < dimX; k++)
                {
                    using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName_state.Replace("{0}", k.ToString()), true))
                    {
                        //outputfile.write(string.format(provider, "{0} {1} {2} {3} {4}",
                        //    t, x[k], y[k], merror[k], merroru[k]
                        //    ));
                        outputfile.Write(string.Format(provider, "{0} {1} {2} {3}",
                            t, x[k], mError[k], mErrorU[k]
                            ));
                        outputfile.WriteLine();
                        outputfile.Close();
                    }
                }
                for (int k = 0; k < dimY; k++)
                {
                    using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName_obs.Replace("{0}", k.ToString()), true))
                    {
                        outputfile.Write(string.Format(provider, "{0} {1}",
                            t, y[k]
                            ));
                        outputfile.WriteLine();
                        outputfile.Close();
                    }
                }
            }
        }

        /// <summary>
        /// Generates a bundle of trajectories, applies the CMN and UK filters, calculates statistics for estimate errors and saves it to files
        /// </summary> 
        /// <param name="n">Number of trajectories</param>
        /// <param name="folderName">Output folder name</param>
        /// <param name="doCalculateUKF">(optional, default = true) if true, UKF and CMNF estimates are calculated, if false - only CMNF </param>
        public void GenerateBundle(int n, string folderName, bool doCalculateUKF = true)
        {
            Console.WriteLine($"GenerateBundle");
            string fileName = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeMany));
            //DiscreteScalarModel[] modelsEst = new DiscreteScalarModel[N];
            DiscreteVectorModel[] modelsEst = new DiscreteVectorModel[n];
            int dimX = X0().Count;

            for (int i = 0; i < n; i++)
            {
                if (i % 1000 == 0) // inform every 1000-th trajectory
                    Console.WriteLine($"model {i}");
                modelsEst[i] = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
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

            Console.WriteLine($"calculate estimates");
            for (int t = 0; t < T; t++)
            {
                Console.WriteLine($"t={t}");
                Vector<double>[] x = new Vector<double>[n];
                Vector<double>[] y = new Vector<double>[n];

                for (int i = 0; i < n; i++)
                {
                    x[i] = modelsEst[i].Trajectory[t][0];
                    y[i] = modelsEst[i].Trajectory[t][1];
                    xHat[i] = CMNF.Step(t, y[i], xHat[i]);
                }

                Vector<double> mx = x.Average();
                Matrix<double> Dx = Exts.Cov(x, x);

                Vector<double> mxHat = xHat.Average();

                Vector<double> mError = (x.Subtract(xHat)).Average();
                Matrix<double> DError = Exts.Cov(x.Subtract(xHat), x.Subtract(xHat));


                Vector<double> mErrorU = Vector<double>.Build.Dense(dimX, 0);
                Matrix<double> DErrorU = Matrix<double>.Build.Dense(dimX, dimX, 0);

                Vector<double> mxHatU = Vector<double>.Build.Dense(dimX, 0);
                Matrix<double> mPHatU = Matrix<double>.Build.Dense(dimX, dimX, 0);

                if (doCalculateUKF)
                {
                    for (int i = 0; i < n; i++)
                    {
                        if (useSimpleModel)
                        {
                            UKF.Step(Phi1, Psi1, DW, DNu, t, y[i], xHatU[i], PHatU[i], out Vector<double> _xHatU, out Matrix<double> _PHatU);
                            xHatU[i] = _xHatU;
                            PHatU[i] = _PHatU;
                        }
                        else
                        {
                            UKF.Step(Phi1, Phi2, Psi1, Psi2, MW, DW, MNu, DNu, t, y[i], xHatU[i], PHatU[i], out Vector<double> _xHatU, out Matrix<double> _PHatU);
                            xHatU[i] = _xHatU;
                            PHatU[i] = _PHatU;
                        }
                    }

                    mErrorU = (x.Subtract(xHatU)).Average();
                    DErrorU = Exts.Cov(x.Subtract(xHatU), x.Subtract(xHatU));

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
                            mxHatU[k], mErrorU[k], DErrorU[k, k], mPHatU[k, k]
                            ));
                        outputfile.WriteLine();
                        outputfile.Close();
                    }
                }
            }
        }


        /// <summary>
        /// Generates N bundles of trajectories, applies the CMN and UK filters. The statistics for estimate errors are calculated by taking average on each bundle, and then on the whole set of bundles. 
        /// Same as GenerateBundle but for larger numbers of trajectories. 
        /// </summary> 
        /// <param name="N">Number of bundles</param>
        /// <param name="n">Number of trajectories</param>
        /// <param name="folderName">Output folder name</param>
        /// <param name="doCalculateUKF">(optional, default = true) if true, UKF and CMNF estimates are calculated, if false - only CMNF </param>
        public void GenerateBundles(int N, int n, string folderName, bool doCalculateUKF = true)
        {
            Console.WriteLine($"GenerateBundles");
            string fileName = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeMany));

            Vector<double>[,] mx = new Vector<double>[N, T];
            Matrix<double>[,] Dx = new Matrix<double>[N, T];

            Vector<double>[,] mxHat = new Vector<double>[N, T];
            Vector<double>[,] mError = new Vector<double>[N, T];
            Matrix<double>[,] DError = new Matrix<double>[N, T];

            Vector<double>[,] mxHatM = new Vector<double>[N, T];
            Matrix<double>[,] mPHatM = new Matrix<double>[N, T];
            Vector<double>[,] mErrorM = new Vector<double>[N, T];
            Matrix<double>[,] DErrorM = new Matrix<double>[N, T];

            Vector<double>[,] mxHatU = new Vector<double>[N, T];
            Matrix<double>[,] mPHatU = new Matrix<double>[N, T];
            Vector<double>[,] mErrorU = new Vector<double>[N, T];
            Matrix<double>[,] DErrorU = new Matrix<double>[N, T];


            int dimX = X0().Count;

            for (int m = 0; m < N; m++)
            {
                Console.WriteLine($"GenerateBundle {m}");
                DiscreteVectorModel[] modelsEst = new DiscreteVectorModel[n];


                for (int i = 0; i < n; i++)
                {
                    if (i % 1000 == 0) // inform every 1000-th trajectory
                        Console.WriteLine($"model {i}");
                    modelsEst[i] = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
                    for (int s = 0; s < T; s++)
                    {
                        modelsEst[i].Step();
                    }
                }

                Vector<double>[] xHat = Enumerable.Repeat(X0Hat, n).ToArray();
                Vector<double>[] xHatM = Enumerable.Repeat(X0Hat, n).ToArray();
                Matrix<double>[] PHatM = Enumerable.Repeat(DX0Hat, n).ToArray();

                Vector<double>[] xHatU = Enumerable.Repeat(X0Hat, n).ToArray();
                Matrix<double>[] PHatU = Enumerable.Repeat(DX0Hat, n).ToArray();

                //Vector<double>[] xHatU2 = Enumerable.Repeat(X0Hat, n).ToArray();
                //Matrix<double>[] PHatU2 = Enumerable.Repeat(DX0Hat, n).ToArray();

                Console.WriteLine($"calculate estimates");
                for (int t = 0; t < T; t++)
                {
                    Console.WriteLine($"t={t}");
                    Vector<double>[] x = new Vector<double>[n];
                    Vector<double>[] y = new Vector<double>[n];

                    for (int i = 0; i < n; i++)
                    {
                        x[i] = modelsEst[i].Trajectory[t][0];
                        y[i] = modelsEst[i].Trajectory[t][1];
                        xHat[i] = CMNF.Step(t, y[i], xHat[i]);
                        (xHatM[i], PHatM[i]) = MCMNF.Step(t, y[i], xHatM[i], PHatM[i]);
                    }

                    mx[m, t] = x.Average();
                    Dx[m, t] = Exts.Cov(x, x);

                    mxHat[m, t] = xHat.Average();
                    mError[m, t] = (x.Subtract(xHat)).Average();
                    DError[m, t] = Exts.Cov(x.Subtract(xHat), x.Subtract(xHat));

                    mxHatM[m, t] = xHatM.Average();
                    mPHatM[m, t] = PHatM.Average();
                    mErrorM[m, t] = (x.Subtract(xHatM)).Average();
                    DErrorM[m, t] = Exts.Cov(x.Subtract(xHatM), x.Subtract(xHatM));


                    mErrorU[m, t] = Vector<double>.Build.Dense(dimX, 0);
                    DErrorU[m, t] = Matrix<double>.Build.Dense(dimX, dimX, 0);

                    mxHatU[m, t] = Vector<double>.Build.Dense(dimX, 0);
                    mPHatU[m, t] = Matrix<double>.Build.Dense(dimX, dimX, 0);


                    if (doCalculateUKF)
                    {
                        for (int i = 0; i < n; i++)
                        {
                            if (useSimpleModel)
                            {
                                UKF.Step(Phi1, Psi1, DW, DNu, t, y[i], xHatU[i], PHatU[i], out Vector<double> _xHatU, out Matrix<double> _PHatU);
                                xHatU[i] = _xHatU;
                                PHatU[i] = _PHatU;
                            }
                            else
                            {
                                UKF.Step(Phi1, Phi2, Psi1, Psi2, MW, DW, MNu, DNu, t, y[i], xHatU[i], PHatU[i], out Vector<double> _xHatU, out Matrix<double> _PHatU);
                                xHatU[i] = _xHatU;
                                PHatU[i] = _PHatU;
                            }
                        }

                        mErrorU[m, t] = (x.Subtract(xHatU)).Average();
                        DErrorU[m, t] = Exts.Cov(x.Subtract(xHatU), x.Subtract(xHatU));

                        mxHatU[m, t] = xHatU.Average();
                        mPHatU[m, t] = PHatU.Average();
                    }
                }
            }
            for (int k = 0; k < dimX; k++)
            {
                using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString())))
                {
                    outputfile.Close();
                }
            }

            for (int t = 0; t < T; t++)
            {
                for (int k = 0; k < dimX; k++)
                {
                    using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString()), true))
                    {
                        outputfile.Write(string.Format(provider, "{0} {1} {2} {3} {4} {5} {6}",
                            t, mx.Average(axis: 1)[t][k], Dx.Average(axis: 1)[t][k, k], mxHat.Average(axis: 1)[t][k], mError.Average(axis: 1)[t][k], DError.Average(axis: 1)[t][k, k], CMNF.KHat[t][k, k]
                            ));
                        outputfile.Write(string.Format(provider, " {0} {1} {2} {3}",
                            mxHatM.Average(axis: 1)[t][k], mErrorM.Average(axis: 1)[t][k], DErrorM.Average(axis: 1)[t][k, k], mPHatM.Average(axis: 1)[t][k, k]
                            ));
                        outputfile.Write(string.Format(provider, " {0} {1} {2} {3}",
                            mxHatU.Average(axis: 1)[t][k], mErrorU.Average(axis: 1)[t][k], DErrorU.Average(axis: 1)[t][k, k], mPHatU.Average(axis: 1)[t][k, k]
                            ));
                        outputfile.WriteLine();
                        outputfile.Close();
                    }
                }
            }
        }


        public void ProcessResults(string dataFolder, string scriptsFolder, string outputFolder)
        {
            Console.WriteLine("Running scripts");

            string[] scriptNamesOne = new string[] { "process_sample", "estimate_sample" };
            string[] scriptNamesMany = new string[] { "process_statistics", "estimate_statistics" };

            string fileNameOne_state = Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeOneState);
            string fileNameOne_obs = Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeOneObs);
            string fileNameMany = Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeMany);

            string scriptOutputFileNameTemplate = Resources.OutputPictureNameTemplate.Replace("{name}", TestFileName);

            foreach (string s in scriptNamesOne)
            {
                Console.WriteLine($"Running {s}");
                RunScript(
                        Path.Combine(scriptsFolder, s + ".py"),
                        new string[] {
                                                Path.Combine(dataFolder, fileNameOne_state),
                                                Path.Combine(dataFolder, fileNameOne_obs),
                                                Path.Combine(outputFolder, scriptOutputFileNameTemplate.Replace("{script}", s))
                                    });
            }
            foreach (string s in scriptNamesMany)
            {
                Console.WriteLine($"Running {s}");
                RunScript(
                        Path.Combine(scriptsFolder, s + ".py"),
                        new string[] {
                                                Path.Combine(dataFolder, fileNameMany),
                                                Path.Combine(outputFolder, scriptOutputFileNameTemplate.Replace("{script}", s))
                                    });
            }
        }

        /// <summary>
        /// Runs a python script to process the calculated test data
        /// </summary>
        /// <param name="scriptName">Python script file path</param>
        /// <param name="scriptParamsTemplates">Script parameters templates array ({0} - number of state vector component)</param>
        private void RunScript(string scriptName, string[] scriptParamsTemplates)
        {
            for (int i = 0; i < X0Hat.Count; i++)
            {
                Python.RunScript(scriptName, scriptParamsTemplates.Select(s => s.Replace("{0}", i.ToString())).ToArray());
            }
            //Python.RunScript(Path.Combine(Settings.Default.ScriptsFolder, "estimate.py"), new string[] { Settings.Default.OutputFolder, "test1_estimateAvg_0.txt", "test1_estimateAvg_0.pdf" });
            //Python.RunScript(Path.Combine(Settings.Default.ScriptsFolder, "estimate.py"), new string[] { Settings.Default.OutputFolder, "test1_estimateAvg_1.txt", "test1_estimateAvg_1.pdf" });

        }

        private string ProcessPics(string picFileNameTemplate, string caption)
        {
            StringBuilder procPics = new StringBuilder();
            for (int i = 0; i < X0Hat.Count; i++)
            {
                string pic = Resources.LatexPictureTemplte.Replace("%file%", picFileNameTemplate.Replace("{0}", i.ToString()));
                pic = pic.Replace("%caption%", caption + (X0Hat.Count > 1 ? $" Компонента {i + 1}." : ""));
                procPics.AppendLine(pic);
                if (i != X0Hat.Count - 1)
                    procPics.AppendLine(@"\vspace{2em}");
            }

            return procPics.ToString();

        }

        public void GenerateReport(string templateFolderName, string folderName)
        {
            Dictionary<string, string> replacements = new Dictionary<string, string>();
            replacements.Add("%Title%", TestName);
            replacements.Add("%phi1%", Phi1_latex == null ? "" : Phi1_latex.ToLatex());
            replacements.Add("%phi2%", Phi2_latex == null ? "" : Phi2_latex.ToLatex());
            replacements.Add("%psi1%", Psi1_latex == null ? "" : Psi1_latex.ToLatex());
            replacements.Add("%psi2%", Psi2_latex == null ? "" : Psi2_latex.ToLatex());
            replacements.Add("%P_w%", P_W);
            replacements.Add("%P_nu%", P_Nu);
            replacements.Add("%P_eta%", P_Eta);

            string picTemplate = Resources.OutputPictureNameTemplate.Replace("{name}", TestFileName);

            replacements.Add("%figs_process%", ProcessPics(picTemplate.Replace("{script}", "process_statistics"), "Статистика процесса."));
            replacements.Add("%figs_filter%", ProcessPics(picTemplate.Replace("{script}", "estimate_statistics"), "Результаты фильтрации."));
            replacements.Add("%figs_process_sample%", ProcessPics(picTemplate.Replace("{script}", "process_sample"), "Пример траектории процесса и наблюдений."));
            replacements.Add("%figs_estimate_sample%", ProcessPics(picTemplate.Replace("{script}", "estimate_sample"), "Пример оценки УМНФ и UKF."));


            string template = File.ReadAllText(Path.Combine(templateFolderName, "modelling_dynamictemplate.tex"), Encoding.Default);
            foreach (var pair in replacements)
            {
                template = template.Replace(pair.Key, pair.Value);
            }
            Directory.CreateDirectory(folderName);
            File.WriteAllText(Path.Combine(folderName, TestFileName + ".tex"), template, Encoding.Default);
        }


    }
}
