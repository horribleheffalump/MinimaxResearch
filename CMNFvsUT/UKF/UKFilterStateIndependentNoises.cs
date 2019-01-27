using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using MathNetExtensions;
using System.IO;
using MathNet.Numerics.Optimization;

namespace UKF
{
    /// <summary>
    /// In UKFilterStateIndependentNoises.cs all methods from UKFilter.cs are duplicated for the case whern Phi_2 and Psi_2 are ident
    /// </summary>
    public partial class UKFilter
    {
        public void EstimateParameters(Func<int, Vector<double>, Vector<double>> Phi,
                                             Func<int, Vector<double>, Vector<double>> Psi,
                                             Matrix<double> Rw,
                                             Matrix<double> Rnu,
                                             Func<Matrix<double>, double> Crit,
                                             int T,
                                             DiscreteVectorModel[] models,
                                             Vector<double> xhat0,
                                             Matrix<double> DX0Hat,
                                             string outputFolder)
        {
            (_, utParamsForecast, utParamsCorrection) = UTParmsOptimize(optimizationMethod, utDefinitionType, Phi, Psi, Rw, Rnu, Crit, T, models, xhat0, DX0Hat, outputFolder);
        }

        static (double, UTParams, UTParams) UTParmsOptimize(OptimizationMethod method,
                                             UTDefinitionType type,
                                             Func<int, Vector<double>, Vector<double>> Phi,
                                             Func<int, Vector<double>, Vector<double>> Psi,
                                             Matrix<double> Rw,
                                             Matrix<double> Rnu,
                                             Func<Matrix<double>, double> Crit,
                                             int T,
                                             DiscreteVectorModel[] models,
                                             Vector<double> xhat0,
                                             Matrix<double> DX0Hat,
                                             string outputFolder)
        {
            (int n, Vector<double> lowerBound, Vector<double> upperBound, Vector<double> initialGuess, string filename) = DefineOptimizationParameters(type, xhat0, string.IsNullOrWhiteSpace(outputFolder) ? null : Path.Combine(outputFolder, "UT_optimization_{type}.txt"));
            double min = double.MaxValue;
            Vector<double> argmin = initialGuess;

            switch (method)
            {
                case OptimizationMethod.RandomShoot:
                    var OptimumRandom = RandomOptimizer.Minimize((x) => CalculateSampleCriterion(Phi, Psi, Rw, Rnu, Crit, x, T, models, xhat0, DX0Hat), Exts.Stack(lowerBound, lowerBound), Exts.Stack(upperBound, upperBound), 100, 100, filename);
                    min = OptimumRandom.min;
                    argmin = OptimumRandom.argmin;
                    break;
                case OptimizationMethod.NelderMeed:
                    NelderMeadSimplex optimizer = new NelderMeadSimplex(1e-3, 100);
                    var objective = ObjectiveFunction.Value((x) => CalculateSampleCriterion(Phi, Psi, Rw, Rnu, Crit, x, T, models, xhat0, DX0Hat));
                    try
                    {
                        var optimumNM = optimizer.FindMinimum(objective, Exts.Stack(initialGuess, initialGuess));
                        min = optimumNM.FunctionInfoAtMinimum.Value;
                        argmin = optimumNM.MinimizingPoint;
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine($"Optimizer faild, using the initail guess ({e.Message})");
                        argmin = Exts.Stack(initialGuess, initialGuess);
                    }

                    break;

            }
            return (min, new UTParams(xhat0.Count, argmin.Take(n).ToArray()), new UTParams(xhat0.Count, argmin.Skip(n).Take(n).ToArray()));
        }

        public static double CalculateSampleCriterion(Func<int, Vector<double>, Vector<double>> Phi,
                                     Func<int, Vector<double>, Vector<double>> Psi,
                                     Matrix<double> Rw,
                                     Matrix<double> Rnu,
                                     Func<Matrix<double>, double> Crit,
                                     Vector<double> P,
                                     int T,
                                     DiscreteVectorModel[] models,
                                     Vector<double> xhat0,
                                     Matrix<double> DX0Hat
                                     )
        {
            (UTParams p1, UTParams p2) = SampleVectorToUTParams(P, xhat0.Count);
            return CalculateCriterionValue(Phi, Psi, Rw, Rnu, Crit, p1, p2, T, models, xhat0, DX0Hat);
        }

          public static double CalculateCriterionValue(Func<int, Vector<double>, Vector<double>> Phi,
                                                     Func<int, Vector<double>, Vector<double>> Psi,
                                                     Matrix<double> Rw,
                                                     Matrix<double> Rnu,
                                                     Func<Matrix<double>, double> Crit,
                                                     UTParams p1,
                                                     UTParams p2,
                                                     int T,
                                                     DiscreteVectorModel[] models,
                                                     Vector<double> xhat0,
                                                     Matrix<double> DX0Hat
                                                     )
        {
            double crit = 0;
            try
            {
                int N = models.Count();

                Vector<double>[] xHatU = models.Select(x => xhat0).ToArray();
                Matrix<double>[] PHatU = models.Select(x => DX0Hat).ToArray();
                //Vector<double> PHatU = Vector<double>.Build.Dense(N, DX0Hat);



                for (int t = 1; t < T; t++)
                {
                    for (int i = 0; i < N; i++)
                    {
                        (xHatU[i], PHatU[i]) = Step(Phi, Psi, Rw, Rnu, p1, p2, t, models[i].Trajectory[t][1], xHatU[i], PHatU[i]);
                    }

                    Vector<double>[] states = models.Select(x => (x.Trajectory[t][0])).ToArray();
                    Matrix<double> errorUPow2 = Exts.Cov(states.Subtract(xHatU), states.Subtract(xHatU));

                    crit = Crit(errorUPow2);
                }
            }
            catch { crit = double.MaxValue; }
            return crit;
        }


   
        public static (Vector<double>, Matrix<double>) Step(Func<int, Vector<double>, Vector<double>> Phi,
                                Func<int, Vector<double>, Vector<double>> Psi,
                                Matrix<double> Rw,
                                Matrix<double> Rnu,
                                UTParams p1,
                                UTParams p2,
                                int t,
                                Vector<double> y,
                                Vector<double> xHat_,
                                Matrix<double> P_)
        {
            UnscentedTransform.Transform(x => Phi(t, x), xHat_, P_, Rw, p1, out Vector<double> Xtilde, out _, out Matrix<double> Ptilde);
            UnscentedTransform.Transform(x => Psi(t, x), Xtilde, Ptilde, Rnu, p2, out Vector<double> Ytilde, out Matrix<double> PXY, out Matrix<double> PYtilde);
            Matrix<double> K = PXY * PYtilde.Inverse();
            return (Xtilde + K * (y - Ytilde), Ptilde - K * PYtilde * K.Transpose());
        }
     
        public (Vector<double>, Matrix<double>) Step(Func<int, Vector<double>, Vector<double>> Phi,
                                Func<int, Vector<double>, Vector<double>> Psi,
                                Matrix<double> Rw,
                                Matrix<double> Rnu,
                                int t,
                                Vector<double> y,
                                Vector<double> xHat_,
                                Matrix<double> P_)
        {
            return Step(Phi, Psi, Rw, Rnu, utParamsForecast, utParamsCorrection, t, y, xHat_, P_);
        }
    }
}

