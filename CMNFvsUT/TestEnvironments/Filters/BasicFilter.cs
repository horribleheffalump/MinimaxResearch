using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestEnvironments
{
    [Serializable]
    public abstract class BasicFilter
    {
        public string FilterName;
        public string FileName;
        public abstract void Initialize();
        public abstract void InitializeAndTrain();
        public abstract (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat);

        public virtual void SaveParams() { }
        public virtual void LoadParams() { }
    }
}
