using System; 
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using MathNetExtensions;
using System.IO;
using MathNet.Numerics.Optimization;

namespace UKF
{
    public enum OptimizationMethod { RandomShoot, NelderMeed}

    public class UKFilter
    {
        public UTParams utParamsForecast; // parameters of the unscented transform on forecast phase
        public UTParams utParamsCorrection; // parameters of the unscented transform on correction phase
        UTDefinitionType utDefinitionType;
        OptimizationMethod optimizationMethod;

        public UKFilter(UTDefinitionType type = UTDefinitionType.ImplicitAlpha, OptimizationMethod method = OptimizationMethod.NelderMeed)
        {
            utDefinitionType = type;
            optimizationMethod = method;
        }


        public void EstimateParameters(int N1, int N2,
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
            (_, utParamsForecast, utParamsCorrection) = UTParmsOptimize(optimizationMethod, utDefinitionType, Phi, Psi, Rw, Rnu, Crit, T, models, xhat0, DX0Hat, outputFolder);
        }

        static (double, UTParams, UTParams) UTParmsOptimize(OptimizationMethod method, UTDefinitionType type,
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

            int n = 0;
            string filename = Path.Combine(outputFolder, "UT_optimization_{type}.txt");
            Vector<double> lowerBound;
            Vector<double> upperBound;
            Vector<double> initialGuess;
            switch (type)
            {
                case UTDefinitionType.ImplicitAlpha:
                    n = 1;
                    lowerBound = Exts.Vector(1 - 2/xhat0.Count);
                    upperBound = Exts.Vector(1);
                    initialGuess = Exts.Vector(0.5);
                    filename = filename.Replace("{type}", "ImplicitAlpha");
                    break;
                case UTDefinitionType.ImplicitAlphaBetaKappa:
                    n = 3;
                    lowerBound = Exts.Vector(0,0,0);
                    upperBound = Exts.Vector(5,5,5);
                    initialGuess = Exts.Vector(1, 2, 1);
                    filename = filename.Replace("{type}", "ImplicitABK"); break;
                case UTDefinitionType.Explicit:
                    n = 4;
                    lowerBound = Exts.Vector(-10, -10, -10, -10);
                    upperBound = Exts.Vector(10, 10, 10, 10);
                    initialGuess = Exts.Vector((new UTParams(xhat0.Count, 1, 2, 1)).Params);
                    filename = filename.Replace("{type}", "Explicit"); break;
                default:
                    n = 0;
                    lowerBound = null;
                    upperBound = null;
                    initialGuess = null;
                    break;
            }

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
                    var optimumNM = optimizer.FindMinimum(objective, Exts.Stack(initialGuess, initialGuess));
                    min = optimumNM.FunctionInfoAtMinimum.Value;
                    argmin = optimumNM.MinimizingPoint;
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
            int n = P.Count();
            if (n % 2 != 0)
                new ArgumentException("Distribution number must be even", "distribution");

            int L = xhat0.Count;

            double[] s1 = P.Take(n / 2).ToArray();
            UTParams p1 = new UTParams(L, s1);

            double[] s2 = P.Skip(n / 2).ToArray();
            UTParams p2 = new UTParams(L, s2);

            double crit = CalculateCriterionValue(Phi, Psi, Rw, Rnu, Crit, p1, p2, T, models, xhat0, DX0Hat);

            return crit;
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



                for (int t = 0; t < T; t++)
                {
                    for (int i = 0; i < N; i++)
                    {
                        Step(Phi, Psi, Rw, Rnu, p1, p2, t, models[i].Trajectory[t][1], xHatU[i], PHatU[i], out Vector<double> xHatU_i, out Matrix<double> PHatU_i);
                        xHatU[i] = xHatU_i;
                        PHatU[i] = PHatU_i;
                    }

                    Vector<double>[] states = models.Select(x => (x.Trajectory[t][0])).ToArray();
                    Matrix<double> errorUPow2 = Exts.Cov(states.Subtract(xHatU), states.Subtract(xHatU));

                    crit = Crit(errorUPow2);
                }
            }
            catch { crit = double.MaxValue; }
            return crit;
        }

        public static void Step(Func<int, Vector<double>, Vector<double>> Phi, 
                                Func<int, Vector<double>, Vector<double>> Psi,
                                Matrix<double> Rw, 
                                Matrix<double> Rnu,
                                UTParams p1,
                                UTParams p2,
                                int t, 
                                Vector<double> y, 
                                Vector<double> xHat_, 
                                Matrix<double> P_, 
                                out Vector<double> xHat, 
                                out Matrix<double> PHat)
        {
            UnscentedTransform.Transform(x => Phi(t,x), xHat_, P_, Rw, p1, out Vector<double> Xtilde, out _, out Matrix<double>Ptilde);

            UnscentedTransform.Transform(x => Psi(t,x), Xtilde, Ptilde, Rnu, p2, out Vector<double> Ytilde, out Matrix<double> PXY, out Matrix<double> PYtilde);

            //Matrix<double> PXY = wc[0] * (Xi.Column(0) - Xtilde).ToColumnMatrix() * (Upsilon.Column(0) - Ytilde).ToRowMatrix();
            //for (int i = 1; i < 2 * L + 1; i++)
            //{
            //    PXY = PXY + wc[i] * (Xi.Column(i) - Xtilde).ToColumnMatrix() * (Upsilon.Column(i) - Ytilde).ToRowMatrix();
            //}

            Matrix<double> K = PXY * PYtilde.Inverse();

            xHat = Xtilde + K * (y - Ytilde);
            PHat = Ptilde - K * PYtilde * K.Transpose();

            t++;
        }

        public void Step(       Func<int, Vector<double>, Vector<double>> Phi,
                                Func<int, Vector<double>, Vector<double>> Psi,
                                Matrix<double> Rw,
                                Matrix<double> Rnu,
                                int t,
                                Vector<double> y,
                                Vector<double> xHat_,
                                Matrix<double> P_,
                                out Vector<double> xHat,
                                out Matrix<double> PHat)
        {
            Step(Phi, Psi, Rw, Rnu, utParamsForecast, utParamsCorrection, t, y, xHat_, P_, out xHat, out PHat);
        }


    }
}


//public void UT(Func<Vector<double>, Vector<double>> f, Matrix<double> Xi, Matrix<double> R, Vector<double> wm, Vector<double> wc, out Matrix<double> Upsilon, out Vector<double> y, out Matrix<double> Py)
//{
//    Upsilon = f(Xi.Column(0)).ToColumnMatrix();
//    for (int i = 1; i < 2 * L + 1; i++)
//    {
//        Upsilon = Upsilon.Append(f(Xi.Column(i)).ToColumnMatrix());
//    }

//    y = wm[0] * Upsilon.Column(0);
//    for (int i = 1; i < 2 * L + 1; i++)
//    {
//        y = y + wm[i] * Upsilon.Column(i);
//    }

//    Py = wc[0] * (Upsilon.Column(0) - y).ToColumnMatrix() * (Upsilon.Column(0) - y).ToRowMatrix();
//    for (int i = 1; i < 2 * L + 1; i++)
//    {
//        Py = Py + wc[i] * (Upsilon.Column(i) - y).ToColumnMatrix() * (Upsilon.Column(i) - y).ToRowMatrix();
//    }
//    Py = Py + R;

//}
