using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestEnvironments
{
    [Serializable]
    class SingleFilterInfo
    {
        public bool Marked = false;
        public int T;
        public string FilterName;
        public Matrix<double>[] Err; // Filtering errors
        public Matrix<double>[] KHat; // Filtering error covariances

        public SingleFilterInfo(int T, string Filter)
        {
            this.T = T;
            this.FilterName = Filter;
            Err = new Matrix<double>[T];
            KHat = new Matrix<double>[T];
        }
    }
}
