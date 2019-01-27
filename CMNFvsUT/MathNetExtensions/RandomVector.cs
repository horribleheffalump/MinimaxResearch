using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNetExtensions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathNetExtensions
{
    /// <summary>
    /// Random vector with given expextation and covariance matrix generator.
    /// </summary>
    /// <typeparam name="T">Continuous distribution</typeparam>
    public class RandomVector<T> where T: IContinuousDistribution, new()
    {
        private T[] distrs;
        private Matrix<double> Sigma;
        private Vector<double> M;
        public RandomVector(Vector<double> M, Matrix<double> Cov)
        {
            this.M = M;
            Svd<double> Cov_svd = Cov.Svd();
            Sigma = Cov_svd.U * Cov_svd.W.PointwiseSqrt() * Cov_svd.VT;
            //if (Cov.FrobeniusNorm() == 0)
            //{
            //    Sigma = Cov;
            //}
            //else
            //{
            //    Sigma = Cov.Cholesky().Factor;
            //}
            distrs = new T[M.Count];
            for (int i = 0; i < M.Count; i++)
            {
                distrs[i] = new T();
            }
        }

        public Vector<double> Sample()
        {
            double[] s = new double[distrs.Length];
            for (int i = 0; i < distrs.Length; i++)
            {
                s[i] =  (distrs[i].Sample() - distrs[i].Mean) / distrs[i].StdDev;
            }
            return M + Sigma * Exts.Vector(s);
        }
    }
}
