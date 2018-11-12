using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestEnvironments
{
    public abstract class BasicFilter
    {
        public abstract void Initialize();
        public abstract (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat);
    }
}
