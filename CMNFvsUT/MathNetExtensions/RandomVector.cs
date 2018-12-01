using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathNetExtensions
{
    public class RandomVector<T> where T: IContinuousDistribution, new()
    {
        private T[] distrs;
        private Matrix<double> Sigma;
        private Vector<double> M;
        public RandomVector(Vector<double> M, Matrix<double> Cov)
        {
            this.M = M;
            Sigma = Cov.Cholesky().Factor;
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
