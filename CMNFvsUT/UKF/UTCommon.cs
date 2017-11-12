using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
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
        /// <para>L - dimention of x, sqrt(P) - Chlesky decomposition </para>
        /// </summary>
        /// <param name="x">Mean</param>
        /// <param name="P">Covariance</param>
        /// <param name="lambda">Spread parameter</param>
        /// <returns></returns>
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
    /// Unscented transform parameters optimization type. 
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
    public class UTParams
    {
        public double Lambda;
        public Vector<double> Wm;
        public Vector<double> Wc;
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
            Wm = Vector<double>.Build.Dense(2 * L + 1, (1.0 - alpha0) / 4.0);
            Wm[0] = alpha0;
            Wc = Vector<double>.Build.Dense(2 * L + 1, (1.0 - alpha0) / 4.0);
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
            Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wm[0] = Lambda / (Lambda + L);

            Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
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
            Wm = Vector<double>.Build.Dense(2 * L + 1, wi);
            Wm[0] = wm0;

            Wc = Vector<double>.Build.Dense(2 * L + 1, wi);
            Wc[0] = wc0;
        }
    }

}
