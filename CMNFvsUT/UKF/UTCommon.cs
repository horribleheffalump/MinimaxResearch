using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace UKF
{
    /// <summary>
    /// Static sigma points generator
    /// </summary>
    public static class SigmaPoints
    {
        /// <summary>
        /// <para>Generates sigma-points given the mean x, covarianve P and spread parameter lambda</para> 
        /// <para>- Xi_0 = x</para>
        /// <para>- Xi_i = x + sqrt((L+lambda)P)_i, i = 1,...,L</para>
        /// <para>- Xi_i = x - sqrt((L+lambda)P)_{i-L}, i = L+1,...,2L </para>
        /// <para>L - dimention of x, sqrt(P) - Cholesky decomposition </para>
        /// </summary>
        /// <param name="x">Mean</param>
        /// <param name="P">Covariance</param>
        /// <param name="lambda">Spread parameter</param>
        /// <returns>Matrix of sigma points</returns>
        public static Matrix<double> Generate(Vector<double> x, Matrix<double> P, double lambda)
        {
            int L = x.Count;

            Matrix<double> Sqrt = (Math.Sqrt(lambda + L)) * P.Cholesky().Factor.Transpose();

            Matrix<double> Xi = x.ToColumnMatrix();
            for (int i = 0; i < L; i++)
            {
                Xi = Xi.Append((x + Sqrt.Column(i)).ToColumnMatrix());
            }
            for (int i = 0; i < L; i++)
            {
                Xi = Xi.Append((x - Sqrt.Column(i)).ToColumnMatrix());
            }

            return Xi;
        }
    }

    /// <summary>
    /// The unscented transform
    /// </summary>
    public static class UnscentedTransform
    {
        /// <summary>
        /// The unscented transform for y = Phi(x) + nu
        /// </summary>
        /// <param name="Phi">Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu</param>
        /// <param name="mX">Mean of the transformed random variable</param>
        /// <param name="dX">Cov of the transformed random variable</param>
        /// <param name="dNu">Cov of the additive random variable</param>
        /// <param name="p">Parameters of the unscented transform</param>
        /// <param name="y">Returns: approximated mean of the transformed variable</param>
        /// <param name="Kxy">Returns: approximated cross-covariance of the initial and the transformed variable</param>
        /// <param name="Kyy">Returns: approximated covariance of the transormed variable</param>
        public static void Transform(Func<Vector<double>, Vector<double>> Phi, Vector<double> mX, Matrix<double> dX, Matrix<double> dNu, UTParams p, out Vector<double> y, out Matrix<double> Kxy, out Matrix<double> Kyy)
        {
            int L = mX.Count;

            Matrix<double> Xi = SigmaPoints.Generate(mX, dX, p.Lambda);

            Matrix<double> Upsilon = Phi(Xi.Column(0)).ToColumnMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {
                Upsilon = Upsilon.Append(Phi(Xi.Column(i)).ToColumnMatrix());
            }

            y = p.Wm[0] * Upsilon.Column(0);
            for (int i = 1; i < 2 * L + 1; i++)
            {
                y = y + p.Wm[i] * Upsilon.Column(i);
            }


            Kxy = p.Wc[0] * (Xi.Column(0) - mX).ToColumnMatrix() * (Upsilon.Column(0) - y).ToRowMatrix();
            Kyy = p.Wc[0] * (Upsilon.Column(0) - y).ToColumnMatrix() * (Upsilon.Column(0) - y).ToRowMatrix();

            Matrix<double> Z = Xi.Stack(Upsilon);
            Vector<double> MZ = mX.Stack(y);
            Matrix<double> PFull = p.Wc[0] * (Z.Column(0) - MZ).ToColumnMatrix() * (Z.Column(0) - MZ).ToRowMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {
                PFull = PFull + p.Wc[i] * (Z.Column(i) - MZ).ToColumnMatrix() * (Z.Column(i) - MZ).ToRowMatrix();
                Kxy = Kxy + p.Wc[i] * (Xi.Column(i) - mX).ToColumnMatrix() * (Upsilon.Column(i) - y).ToRowMatrix();
                Kyy = Kyy + p.Wc[i] * (Upsilon.Column(i) - y).ToColumnMatrix() * (Upsilon.Column(i) - y).ToRowMatrix();
            }
            Kyy = Kyy + dNu;
            //Py = Py + R;
            //return PFull;
        }
    }

    /// <summary>
    /// Unscented transform parameters optimization method.
    /// </summary>
    public enum OptimizationMethod { RandomShoot, NelderMeed }

    /// <summary>
    /// Unscented transform parameters definition type.
    /// </summary>
    public enum UTDefinitionType { ImplicitAlpha, ImplicitAlphaBetaKappa, Explicit }

    /// <summary>
    /// <para>Parameters of the unscented transform:</para>
    /// <para>- Lambda - scaling parameter</para>
    /// <para>- Wm - weights for sample mean</para>
    /// <para>- Wc - weights for sample covariance</para>
    /// <para>Three ways to define:</para>
    /// <para>- explicitly,</para>
    /// <para>- implicitly with one parameter alpha0,</para>
    /// <para>- implicitly with three parameters alpha, beta, kappa</para>
    /// </summary>
    [Serializable]
    public class UTParams
    {
        public double Lambda;
        public double[] Wm;
        public double[] Wc;
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="L">Dimension of the transformed random variable</param>
        /// <param name="p">Input params to define the unscented transform params. The appropriate method to define UT params is chosen depending on the size of the array.</param>
        public UTParams(int L, params double[] p)
        {
            int n = p.Count();
            switch (n)
            {
                case 1: SetUTParams(L, p[0]); break;
                case 3: SetUTParams(L, p[0], p[1], p[2]); break;
                case 4: SetUTParams(L, p[0], p[1], p[2], p[3]); break;
            }
        }

        /// <summary>
        /// The parameters of the unscented transform are defined by a single parameter alpha0:
        /// </summary>
        /// <param name="L">Dimension of the transformed random variable</param>
        /// <param name="alpha0">Alpha0 - weight of the central points for both sample mean and cov </param>
        public void SetUTParams(int L, double alpha0)
        {
            Lambda = 2.0 / (1 - alpha0) - L;
            //Wm = Vector<double>.Build.Dense(2 * L + 1, (1.0 - alpha0) / 4.0);
            Wm = Enumerable.Repeat((1.0 - alpha0) / 4.0, 2 * L + 1).ToArray();
            Wm[0] = alpha0;
            //Wc = Vector<double>.Build.Dense(2 * L + 1, (1.0 - alpha0) / 4.0);
            Wc = Enumerable.Repeat((1.0 - alpha0) / 4.0, 2 * L + 1).ToArray(); 
            Wc[0] = alpha0;
        }

        /// <summary>
        /// The parameters of the unscented transform are defined by three parameters: alpha, beta, kappa
        /// </summary>
        /// <param name="L">Dimension of the transformed random variable</param>
        /// <param name="alpha">Alpha - determines the spread of the sigma points around the mean of the transformed random variable (small positive value 0 \leq alpha \leq 10e-4)</param>
        /// <param name="beta">Beta - is used to incorporate prior knowledge of the distribution of the transformed random variable (for Gaussian b = 2 is optimal)</param>
        /// <param name="kappa">Kappa - is a secondary scaling parameter</param>
        public void SetUTParams(int L, double alpha, double beta, double kappa)
        {
            Lambda = Math.Pow(alpha, 2.0) * (L + kappa) - L;
            //Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wm = Enumerable.Repeat(0.5 / (Lambda + L), 2 * L + 1).ToArray();
            Wm[0] = Lambda / (Lambda + L);

            //Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wc = Enumerable.Repeat(0.5 / (Lambda + L), 2 * L + 1).ToArray();
            Wc[0] = Lambda / (Lambda + L) + 1.0 - Math.Pow(alpha, 2.0) + beta;
        }

        /// <summary>
        /// Explicit definition of the unscented transform parameters
        /// </summary>
        /// <param name="L">Dimension of the transformed random variable</param>
        /// <param name="lambda">Scaling parameter</param>
        /// <param name="wm0">Central point weight for the sample mean</param>
        /// <param name="wc0">Central point weight for the sample cov</param>
        /// <param name="wi">Non-central points weight for sample mean and cov</param>
        public void SetUTParams(int L, double lambda, double wm0, double wc0, double wi)
        {
            Lambda = lambda;
            //Wm = Vector<double>.Build.Dense(2 * L + 1, wi);
            Wm = Enumerable.Repeat(wi, 2 * L + 1).ToArray();
            Wm[0] = wm0;

            //Wc = Vector<double>.Build.Dense(2 * L + 1, wi);
            Wc = Enumerable.Repeat(wi, 2 * L + 1).ToArray();
            Wc[0] = wc0;
        }

        /// <summary>
        /// Get unscented transformation paramters as array
        /// </summary>
        public double[] Params
        {
            get
            {
                return new double[4] { Lambda, Wm[0], Wc[0], Wm[1] };
            }
        }
    }

}
