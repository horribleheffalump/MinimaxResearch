
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
using MathNetExtensions;
using System.Threading.Tasks;
using TestEnvironments.Properties;
using System.Threading;
using System.Runtime.Serialization.Formatters.Binary;
using System.Runtime.Serialization;

namespace TestEnvironments
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

        //public string[] Phi1_latex;
        //public string[][] Phi2_latex;
        //public string[] Psi1_latex;
        //public string[][] Psi2_latex;

        //public string P_W;
        //public string P_Nu;
        //public string P_Eta;


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

        public Func<DiscreteVectorModel> ModelGenerator;
        //private NumberFormatInfo provider;

        public BasicFilter[] Filters;

        private void HandleNulls()
        {
            if (Phi2 == null)
                Phi2 = new Func<int, Vector<double>, Matrix<double>>((s, x) => Matrix<double>.Build.Dense(1, 1, 1.0));
            if (Psi2 == null)
                Psi2 = new Func<int, Vector<double>, Matrix<double>>((s, x) => Matrix<double>.Build.Dense(1, 1, 1.0));
            if (MW == null)
                MW = W(0) * 0.0;
            if (MNu == null)
                MNu = Nu(0) * 0.0;
        }

        /// <summary>
        /// Initializes the test environment by calculating the statistics for CMN and UT filters
        /// </summary>
        /// <param name="t">time interval right margin (number of steps)</param>
        /// <param name="n">number of trajectories</param>
        /// <param name="doCalculateUKF"></param>
        /// <param name="doCalculateUKFStepwise"></param>
        /// <param name="outputFolder"></param>
        //public void Initialize(int t, int n, bool doCalculateUKF, bool doCalculateUKFStepwise, bool doOptimizeWithRandomShoot, string outputFolder)
        public void Initialize(int t, int n, int nMCMNF, string outputFolder, List<(FilterType type, string fileName)> filters, bool save = false, bool load = false, Func<DiscreteVectorModel> ModelGenerator = null)
        {
            if (ModelGenerator == null)
            {
                this.ModelGenerator = () =>
                       {
                           DiscreteVectorModel model = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
                           for (int s = 1; s <= T; s++)
                           {
                               model.Step();
                           }
                           return model;
                       };
            }
            else
            {
                this.ModelGenerator = ModelGenerator;
            }
            Console.WriteLine("Init");
            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";
            HandleNulls();

            T = t;

            // generate models for filters parameters fitting/optimization if we are not importing them from file
            DiscreteVectorModel[] models = new DiscreteVectorModel[0];
            if (!load)
            {
                models = new DiscreteVectorModel[n];
                for (int i = 0; i < n; i++)
                {
                    if (i % 1000 == 0) // inform every 1000-th trajectory
                        Console.WriteLine($"model {i}");
                    models[i] = this.ModelGenerator();
                    //models[i] = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
                    //for (int s = 1; s < T; s++)
                    //{
                    //    models[i].Step();
                    //}
                }
            }
            Filters = new BasicFilter[filters.Count()];
            for (int j = 0; j < filters.Count(); j++)
            {
                if (filters[j].type == FilterType.CMNF)
                {
                    CMNFWrapper CMNF = new CMNFWrapper
                    {
                        FileName = filters[j].fileName,
                        T = T,
                        X0Hat = X0Hat,
                        Models = models,
                        Xi = Xi,
                        Zeta = Zeta
                    };
                    Filters[j] = CMNF;
                }
                if (filters[j].type == FilterType.MCMNF)
                {
                    int nTrain = nMCMNF;
                    if (nTrain == 0)
                    {
                        Console.WriteLine($"MCMNF stepwise train set size is not provided, using the integral train set size of {n} instead");
                        nTrain = n;
                    }
                    MCMNFWrapper MCMNF = new MCMNFWrapper
                    {
                        FileName = filters[j].fileName,
                        T = T,
                        N = nTrain,
                        X0Hat = X0Hat,
                        Models = models,
                        Xi = Xi,
                        Zeta = Zeta,
                        Phi1 = Phi1,
                        Phi2 = Phi2,
                        Psi1 = Psi1,
                        Psi2 = Psi2,
                        W = W,
                        Nu = Nu
                    };
                    Filters[j] = MCMNF;
                }
                if (new[] { FilterType.UKFNoOptimization, FilterType.UKFIntegral, FilterType.UKFIntegralRandomShoot, FilterType.UKFStepwise, FilterType.UKFStepwiseRandomShoot }.Contains(filters[j].type))
                {
                    UKFWrapper UKF = new UKFWrapper
                    {
                        FileName = filters[j].fileName,
                        T = T,
                        X0Hat = X0Hat,
                        Models = models,
                        Phi1 = Phi1,
                        Phi2 = Phi2,
                        Psi1 = Psi1,
                        Psi2 = Psi2,
                        MW = MW,
                        DW = DW,
                        MNu = MNu,
                        DNu = DNu,
                        Crit = x => x.Trace(),//x[0, 0],
                        DX0Hat = DX0Hat,
                        outputFolder = outputFolder
                    };

                    if (filters[j].type == FilterType.UKFNoOptimization)
                    {
                        UKF.doOptimize = false;
                        UKF.doCalculateUKFStepwise = false;
                    }
                    else
                        UKF.doOptimize = true;

                    if (new[] { FilterType.UKFIntegralRandomShoot, FilterType.UKFStepwiseRandomShoot }.Contains(filters[j].type))
                        UKF.doOptimizeWithRandomShoot = true;
                    else
                        UKF.doOptimizeWithRandomShoot = false;

                    if (new[] { FilterType.UKFStepwise, FilterType.UKFStepwiseRandomShoot }.Contains(filters[j].type))
                        UKF.doCalculateUKFStepwise = true;
                    else
                        UKF.doCalculateUKFStepwise = false;

                    Filters[j] = UKF;
                }
            }


            foreach (var f in Filters)
            {
                if (load)
                {
                    f.Initialize();
                    f.LoadParams();
                }
                else
                {
                    f.InitializeAndTrain();
                }
                if (save)
                    f.SaveParams();
            }
            models = null;
            //GC.Collect();
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

            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            T = t;

            DiscreteVectorModel[] models = new DiscreteVectorModel[n];
            for (int i = 0; i < n; i++)
            {
                if (i % 1000 == 0) // inform every 1000-th trajectory
                    Console.WriteLine($"model {i}");
                //Console.WriteLine($"model {i}");
                models[i] = ModelGenerator();
                //models[i] = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
                //for (int s = 1; s < T; s++)
                //{
                //    models[i].Step();
                //}
            }

            for (int k = 0; k < models[0].State.Count; k++)
            {
                string fileName = Path.Combine(outputFolder, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeBulk));
                using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString())))
                {
                    for (int s = 1; s < T; s++)
                    {
                        outputfile.Write(string.Format(provider, "{0} ", s));
                    }
                    outputfile.WriteLine();
                    for (int i = 0; i < n; i++)
                    {
                        for (int s = 1; s < T; s++)
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

        public void GenerateOne(string folderName, int? n = null)
        {
            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";
            string fileName_state = "";
            string fileName_obs = "";
            if (n == null)
            {
                fileName_state = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeOneState));
                fileName_obs = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeOneObs));
            }
            else
            {
                fileName_state = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeOneState + "_" + n.ToString()));
                fileName_obs = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeOneObs + "_" + n.ToString()));
            }

            DiscreteVectorModel modelEst = ModelGenerator();
            //DiscreteVectorModel modelEst = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
            //for (int s = 1; s < T; s++)
            //{
            //    modelEst.Step();
            //}

            int dimX = modelEst.Trajectory[1][0].Count;
            int dimY = modelEst.Trajectory[1][1].Count;

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

            Vector<double>[] xHat = new Vector<double>[Filters.Count()];
            Matrix<double>[] KHat = new Matrix<double>[Filters.Count()];
            Vector<double>[] mError = new Vector<double>[Filters.Count()];
            for (int j = 0; j < Filters.Count(); j++)
            {
                xHat[j] = X0Hat;
                KHat[j] = DX0Hat;
            }

            for (int t = 1; t < T; t++) // start from 1 because for 0 we do not have observations
            {
                var x = modelEst.Trajectory[t][0];
                var y = modelEst.Trajectory[t][1];

                for (int j = 0; j < Filters.Count(); j++)
                {
                    (xHat[j], KHat[j]) = Filters[j].Step(t, y, xHat[j], KHat[j]);
                    mError[j] = x - xHat[j];
                }

                for (int k = 0; k < dimX; k++)
                {
                    using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName_state.Replace("{0}", k.ToString()), true))
                    {
                        outputfile.Write(string.Format(provider, "{0} {1}",
                            t, x[k]
                            ));
                        for (int j = 0; j < Filters.Count(); j++)
                        {
                            outputfile.Write(string.Format(provider, " {0}",
                                mError[j][k]
                                ));
                        }
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

        ///// <summary>
        ///// Generates a bundle of trajectories, applies the filters, calculates statistics for estimate errors and saves it to files
        ///// </summary> 
        ///// <param name="n">Number of trajectories</param>
        ///// <param name="folderName">Output folder name</param>
        //public void GenerateBundle(int n, string folderName)
        //{
        //    Console.WriteLine($"GenerateBundle");
        //    string fileName = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeMany));
        //    //DiscreteScalarModel[] modelsEst = new DiscreteScalarModel[N];
        //    DiscreteVectorModel[] modelsEst = new DiscreteVectorModel[n];
        //    int dimX = X0().Count;

        //    for (int i = 0; i < n; i++)
        //    {
        //        if (i % 1000 == 0) // inform every 1000-th trajectory
        //            Console.WriteLine($"model {i}");
        //        modelsEst[i] = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
        //        for (int s = 1; s < T; s++)
        //        {
        //            modelsEst[i].Step();
        //        }
        //    }

        //    for (int k = 0; k < dimX; k++)
        //    {
        //        using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString())))
        //        {
        //            outputfile.Close();
        //        }
        //    }

        //    Vector<double>[][] xHat = new Vector<double>[Filters.Count()][];
        //    Matrix<double>[][] KHat = new Matrix<double>[Filters.Count()][];
        //    Vector<double>[] mxHat = new Vector<double>[Filters.Count()];
        //    Matrix<double>[] mKHat = new Matrix<double>[Filters.Count()];
        //    Vector<double>[] mError = new Vector<double>[Filters.Count()];
        //    Matrix<double>[] DError = new Matrix<double>[Filters.Count()];
        //    for (int j = 0; j < Filters.Count(); j++)
        //    {
        //        xHat[j] = Enumerable.Repeat(X0Hat, n).ToArray();
        //        KHat[j] = Enumerable.Repeat(DX0Hat, n).ToArray();
        //    }
        //    Console.WriteLine($"calculate estimates");
        //    for (int t = 1; t < T; t++)
        //    {
        //        Console.WriteLine($"t={t}");
        //        Vector<double>[] x = new Vector<double>[n];
        //        Vector<double>[] y = new Vector<double>[n];

        //        for (int j = 0; j < Filters.Count(); j++)
        //        {
        //            for (int i = 0; i < n; i++)
        //            {
        //                x[i] = modelsEst[i].Trajectory[t][0];
        //                y[i] = modelsEst[i].Trajectory[t][1];
        //                (xHat[j][i], KHat[j][i]) = Filters[j].Step(t, y[i], xHat[j][i], KHat[j][i]);

        //            }
        //            mxHat[j] = xHat[j].Average();
        //            mKHat[j] = KHat[j].Average();
        //            mError[j] = (x.Subtract(xHat[j])).Average();
        //            DError[j] = Exts.Cov(x.Subtract(xHat[j]), x.Subtract(xHat[j]));
        //        }

        //        Vector<double> mx = x.Average();
        //        Matrix<double> Dx = Exts.Cov(x, x);


        //        for (int k = 0; k < dimX; k++)
        //        {
        //            using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString()), true))
        //            {
        //                outputfile.Write(string.Format(provider, "{0} {1} {2}",
        //                    t, mx[k], Dx[k, k]
        //                    ));
        //                for (int j = 0; j < Filters.Count(); j++)
        //                {
        //                    outputfile.Write(string.Format(provider, " {0} {1} {2} {3}",
        //                    mxHat[j][k], mError[j][k], DError[j][k, k], mKHat[j][k, k]
        //                    ));
        //                }
        //                outputfile.WriteLine();
        //                outputfile.Close();
        //            }
        //        }
        //    }
        //}


        /// <summary>
        /// Generates N bundles of trajectories, applies the CMN and UK filters. The statistics for estimate errors are calculated by taking average on each bundle, and then on the whole set of bundles. 
        /// Same as GenerateBundle but for larger numbers of trajectories. 
        /// </summary> 
        /// <param name="N">Number of bundles</param>
        /// <param name="n">Number of trajectories</param>
        /// <param name="folderName">Output folder name</param>
        /// <param name="doCalculateUKF">(optional, default = true) if true, UKF and CMNF estimates are calculated, if false - only CMNF </param>
        public void GenerateBundles(int N, int n, string folderName, bool parallel = false, int parallel_degree = 1, bool doSaveBin = true, bool doSaveText = true)
        {
            Console.WriteLine($"GenerateBundles");
            string fileName = Path.Combine(folderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeMany));

            List<ProcessInfo> processInfos = new List<ProcessInfo>();

            int dimX = X0().Count;

            if (parallel)
            {
                AsyncCalculatorPlanner acp = new MathNetExtensions.AsyncCalculatorPlanner(N, parallel_degree, () => CalculateBundle(n, new ProcessInfo(T, Filters.Count(), X0Hat, DX0Hat)));
                processInfos = acp.DoCalculate().Select(x => x as ProcessInfo).ToList();
            }
            else
            {
                for (int m = 0; m < N; m++)
                {
                    Console.WriteLine($"GenerateBundle {m}");
                    ProcessInfo p = new ProcessInfo(T, Filters.Count(), X0Hat, DX0Hat);
                    p = CalculateBundle(n, p);
                    processInfos.Add(p);
                }
            }

            ProcessInfo totalProcessInfo;
            if (processInfos.Count > 1)
            {
                totalProcessInfo = new ProcessInfo(processInfos.ToArray());
            }
            else
            {
                totalProcessInfo = processInfos[0];
            }
            if (doSaveText)
                totalProcessInfo.SaveToText(fileName);
            if (doSaveBin)
            {
                string fileNameBin = Path.Combine(folderName, Resources.OutputFileNameBinTemplate.Replace("{name}", TestFileName).Replace("{rnd}", Guid.NewGuid().ToString()));
                totalProcessInfo.SaveToFile(fileNameBin);
            }
        }


        /// <summary>
        /// Imports and aggregates the data from prevoiously generated bundles of trajectories.
        /// The statistics for estimate errors are calculated by taking average the whole set of bundles. 
        /// </summary> 
        /// <param name="folderName">Output folder name</param>
        public void Aggregate(string inputFolderName, string outputFolderName, bool doSaveBin = true, bool doSaveText = true)
        {
            Console.WriteLine($"ImportBundles");
            List<ProcessInfo> pinfos = new List<ProcessInfo>();
            foreach (var f in Directory.EnumerateFiles(inputFolderName))
            {
                string result = $"Importing {f}: ";
                try
                {
                    ProcessInfo p = ProcessInfo.LoadFromFile(f);
                    pinfos.Add(p);
                    result += $"succeded, trajectories count: {p.Count}";
                }
                catch (Exception e)
                {
                    result += $"failed {e.Message}";
                }
                Console.WriteLine(result);
            }
            if (pinfos.Count > 0)
            {
                ProcessInfo totalProcessInfo = new ProcessInfo(pinfos.ToArray());
                Console.WriteLine($"Aggregated data. Trajectories count: {totalProcessInfo.Count}");

                if (doSaveText)
                {
                    string fileName = Path.Combine(outputFolderName, Resources.OutputFileNameTemplate.Replace("{name}", TestFileName).Replace("{type}", Resources.OutputTypeMany));
                    totalProcessInfo.SaveToText(fileName);
                }
                if (doSaveBin)
                {
                    string fileNameBin = Path.Combine(outputFolderName, Resources.OutputFileNameBinTemplate.Replace("{name}", TestFileName).Replace("{rnd}", "total"));
                    totalProcessInfo.SaveToFile(fileNameBin);
                }
            }
            else
                Console.WriteLine("No data found");
        }



        //private void CalculateBundle(int n, int m, Vector<double>[,] mx, Matrix<double>[,] Dx, Vector<double>[][,] mxHat, Matrix<double>[][,] mKHat, Vector<double>[][,] mError, Matrix<double>[][,] DError)
        //{
        //    Console.WriteLine($"GenerateBundle {m}");
        //    DiscreteVectorModel[] modelsEst = new DiscreteVectorModel[n];
        //    for (int i = 0; i < n; i++)
        //    {
        //        if (i % 1000 == 0) // inform every 1000-th trajectory
        //            Console.WriteLine($"model {i}");
        //        modelsEst[i] = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
        //        for (int s = 1; s < T; s++)
        //        {
        //            modelsEst[i].Step();
        //        }
        //    }

        //    Vector<double>[][] xHat = new Vector<double>[Filters.Count()][];
        //    Matrix<double>[][] KHat = new Matrix<double>[Filters.Count()][];
        //    for (int j = 0; j < Filters.Count(); j++)
        //    {
        //        xHat[j] = Enumerable.Repeat(X0Hat, n).ToArray();
        //        KHat[j] = Enumerable.Repeat(DX0Hat, n).ToArray();
        //    }
        //    Console.WriteLine($"calculate estimates");
        //    for (int t = 1; t < T; t++)
        //    {
        //        Console.WriteLine($"t={t}");
        //        Vector<double>[] x = new Vector<double>[n];
        //        Vector<double>[] y = new Vector<double>[n];
        //        for (int j = 0; j < Filters.Count(); j++)
        //        {
        //            for (int i = 0; i < n; i++)
        //            {
        //                x[i] = modelsEst[i].Trajectory[t][0];
        //                y[i] = modelsEst[i].Trajectory[t][1];
        //                (xHat[j][i], KHat[j][i]) = Filters[j].Step(t, y[i], xHat[j][i], KHat[j][i]);
        //            }
        //            mxHat[j][m, t] = xHat[j].Average();
        //            mKHat[j][m, t] = KHat[j].Average();
        //            mError[j][m, t] = (x.Subtract(xHat[j])).Average();
        //            DError[j][m, t] = Exts.Cov(x.Subtract(xHat[j]), x.Subtract(xHat[j]));
        //        }
        //        mx[m, t] = x.Average();
        //        Dx[m, t] = Exts.Cov(x, x);
        //    }

        //}

        private ProcessInfo CalculateBundle(int n, ProcessInfo processInfo)
        {
            //Console.WriteLine($"GenerateBundle {m}");
            DiscreteVectorModel[] modelsEst = new DiscreteVectorModel[n];
            for (int i = 0; i < n; i++)
            {
                //if (i % 1000 == 0) // inform every 1000-th trajectory
                //    Console.WriteLine($"model {i}");
                modelsEst[i] = ModelGenerator();
                //modelsEst[i] = new DiscreteVectorModel(Phi1, Phi2, Psi1, Psi2, W, Nu, X0(), true);
                //for (int s = 1; s < T; s++)
                //{
                //    modelsEst[i].Step();
                //}
            }

            Vector<double>[][] xHat = new Vector<double>[Filters.Count()][];
            Matrix<double>[][] KHat = new Matrix<double>[Filters.Count()][];
            for (int j = 0; j < Filters.Count(); j++)
            {
                xHat[j] = Exts.ArrayOf(X0Hat, n);
                KHat[j] = Exts.ArrayOf(DX0Hat, n);
            }
            Console.WriteLine($"calculate estimates");
            for (int t = 1; t < T; t++)
            {
                Console.WriteLine($"t={t}");
                Vector<double>[] x = new Vector<double>[n];
                Vector<double>[] y = new Vector<double>[n];
                for (int j = 0; j < Filters.Count(); j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        x[i] = modelsEst[i].Trajectory[t][0];
                        y[i] = modelsEst[i].Trajectory[t][1];
                        (xHat[j][i], KHat[j][i]) = Filters[j].Step(t, y[i], xHat[j][i], KHat[j][i]);
                    }
                    processInfo.FilterQualityInfos[j].mxHat[t] = xHat[j].Average().ToColumnMatrix();
                    processInfo.FilterQualityInfos[j].mKHat[t] = KHat[j].Average();
                    processInfo.FilterQualityInfos[j].mError[t] = (x.Subtract(xHat[j])).Average().ToColumnMatrix();
                    processInfo.FilterQualityInfos[j].DError[t] = Exts.Cov(x.Subtract(xHat[j]), x.Subtract(xHat[j]));
                }
                processInfo.mx[t] = x.Average().ToColumnMatrix();
                processInfo.Dx[t] = Exts.Cov(x, x);
            }
            processInfo.Count = n;
            return processInfo;
        }

        public void ProcessResults(string dataFolder, string scriptsFolder, string outputFolder)
        {
            Console.WriteLine("Running scripts");

            string[] scriptNamesOne = new string[] { "process_sample", "estimate_sample" };
            string[] scriptNamesMany = new string[] { "process_statistics", "estimate_statistics", "estimate_statistics_single" };

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

        //public void GenerateReport(string templateFolderName, string folderName)
        //{
        //    Dictionary<string, string> replacements = new Dictionary<string, string>();
        //    replacements.Add("%Title%", TestName);
        //    replacements.Add("%phi1%", Phi1_latex == null ? "" : Phi1_latex.ToLatex());
        //    replacements.Add("%phi2%", Phi2_latex == null ? "" : Phi2_latex.ToLatex());
        //    replacements.Add("%psi1%", Psi1_latex == null ? "" : Psi1_latex.ToLatex());
        //    replacements.Add("%psi2%", Psi2_latex == null ? "" : Psi2_latex.ToLatex());
        //    replacements.Add("%P_w%", P_W);
        //    replacements.Add("%P_nu%", P_Nu);
        //    replacements.Add("%P_eta%", P_Eta);

        //    string picTemplate = Resources.OutputPictureNameTemplate.Replace("{name}", TestFileName);

        //    replacements.Add("%figs_process%", ProcessPics(picTemplate.Replace("{script}", "process_statistics"), "Статистика процесса."));
        //    replacements.Add("%figs_filter%", ProcessPics(picTemplate.Replace("{script}", "estimate_statistics"), "Результаты фильтрации."));
        //    replacements.Add("%figs_process_sample%", ProcessPics(picTemplate.Replace("{script}", "process_sample"), "Пример траектории процесса и наблюдений."));
        //    replacements.Add("%figs_estimate_sample%", ProcessPics(picTemplate.Replace("{script}", "estimate_sample"), "Пример оценки УМНФ и UKF."));


        //    string template = File.ReadAllText(Path.Combine(templateFolderName, "modelling_dynamictemplate.tex"), Encoding.Default);
        //    foreach (var pair in replacements)
        //    {
        //        template = template.Replace(pair.Key, pair.Value);
        //    }
        //    Directory.CreateDirectory(folderName);
        //    File.WriteAllText(Path.Combine(folderName, TestFileName + ".tex"), template, Encoding.Default);
        //}


    }

    public enum FilterType { CMNF, MCMNF, UKFNoOptimization, UKFIntegral, UKFIntegralRandomShoot, UKFStepwise, UKFStepwiseRandomShoot };

    [Serializable]
    public class FilterQualityInfo
    {
        public Matrix<double>[] mError;
        public Matrix<double>[] DError;
        public Matrix<double>[] mxHat;
        public Matrix<double>[] mKHat;

        public FilterQualityInfo(int T, Vector<double> X0Hat, Matrix<double> DX0Hat)
        {
            mError = Exts.ZerosArrayOfShape(X0Hat.ToColumnMatrix(), T);
            DError = Exts.ZerosArrayOfShape(DX0Hat, T);
            mxHat = Exts.ZerosArrayOfShape(X0Hat.ToColumnMatrix(), T);
            mKHat = Exts.ZerosArrayOfShape(DX0Hat, T);
        }

        public FilterQualityInfo(int T, Matrix<double> X0Hat, Matrix<double> DX0Hat) : this(T, X0Hat.Column(0), DX0Hat) { }

    }

    [Serializable]
    public class ProcessInfo
    {
        public int Count;
        public Matrix<double>[] mx;
        public Matrix<double>[] Dx;
        public FilterQualityInfo[] FilterQualityInfos;
        public ProcessInfo(int T, int FiltersCount, Vector<double> X0Hat, Matrix<double> DX0Hat)
        {
            Count = 0;
            mx = Exts.ZerosArrayOfShape(X0Hat.ToColumnMatrix(), T);
            Dx = Exts.ZerosArrayOfShape(DX0Hat, T);
            FilterQualityInfos = new FilterQualityInfo[FiltersCount];
            for (int j = 0; j < FiltersCount; j++)
            {
                FilterQualityInfos[j] = new FilterQualityInfo(T, X0Hat, DX0Hat);
            }


        }

        public static ProcessInfo LoadFromFile(string fileName)
        {
            ProcessInfo p;
            IFormatter formatter = new BinaryFormatter();
            using (Stream stream = new FileStream(fileName, FileMode.Open))
            {
                p = (ProcessInfo)formatter.Deserialize(stream);
                stream.Close();
            }
            return p;
        }

        public ProcessInfo(ProcessInfo[] infos)
        {
            Count = (int)infos.Select(i => (double)i.Count).Sum();
            double[] TrCounts = infos.Select(i => (double)i.Count).ToArray();
            mx = Exts.Average(infos.Select(i => i.mx).ToArray(), TrCounts);
            Dx = Exts.Average(infos.Select(i => i.Dx).ToArray(), TrCounts);
            FilterQualityInfos = new FilterQualityInfo[infos[0].FilterQualityInfos.Count()];
            for (int j = 0; j < FilterQualityInfos.Count(); j++)
            {
                FilterQualityInfos[j] = new FilterQualityInfo(infos[0].FilterQualityInfos[j].mError.Length, infos[0].FilterQualityInfos[j].mError[0], infos[0].FilterQualityInfos[j].DError[0]);
                FilterQualityInfos[j].mError = Exts.Average(infos.Select(i => i.FilterQualityInfos[j].mError).ToArray(), TrCounts);
                FilterQualityInfos[j].DError = Exts.Average(infos.Select(i => i.FilterQualityInfos[j].DError).ToArray(), TrCounts);
                FilterQualityInfos[j].mxHat = Exts.Average(infos.Select(i => i.FilterQualityInfos[j].mxHat).ToArray(), TrCounts);
                FilterQualityInfos[j].mKHat = Exts.Average(infos.Select(i => i.FilterQualityInfos[j].mKHat).ToArray(), TrCounts);
            }
        }

        public void SaveToFile(string fileName)
        {
            IFormatter formatter = new BinaryFormatter();
            using (Stream stream = new FileStream(fileName, FileMode.Create))
            {
                formatter.Serialize(stream, this);
                stream.Close();
            }
        }

        public void SaveToText(string fileName)
        {
            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            for (int k = 0; k < mx[0].RowCount; k++)
            {
                using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString())))
                {
                    outputfile.Close();
                }
            }
            for (int t = 0; t < mx.Length; t++)
            {
                for (int k = 0; k < mx[0].RowCount; k++)
                {
                    using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString()), true))
                    {
                        outputfile.Write(string.Format(provider, "{0} {1} {2}",
                            t, mx[t][k, 0], Dx[t][k, k]
                            ));
                        for (int j = 0; j < FilterQualityInfos.Count(); j++)
                        {
                            outputfile.Write(string.Format(provider, " {0} {1} {2} {3}",
                            FilterQualityInfos[j].mxHat[t][k, 0], FilterQualityInfos[j].mError[t][k, 0], FilterQualityInfos[j].DError[t][k, k], FilterQualityInfos[j].mKHat[t][k, k]
                            ));
                        }
                        outputfile.WriteLine();
                        outputfile.Close();
                    }
                }
            }
        }
    }
}