using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestEnvironments
{
    public class DummyFilter : BasicFilter
    {

        public string outputFolder;
        public Func<int, Vector<double>, Vector<double>, Matrix<double>, (Vector<double>, Matrix<double>)> StepFunction;

        public override void Initialize()
        {
            FilterName = "Dummy";
        }

        public override void InitializeAndTrain()
        {
            Initialize();
        }

        public override (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat)
        {
            return StepFunction(t, y, xHat, kHat);
        }
    }
}
