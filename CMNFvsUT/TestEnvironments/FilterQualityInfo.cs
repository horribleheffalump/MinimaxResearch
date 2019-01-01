using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestEnvironments
{
    public enum FilterType { CMNF, MCMNF, RCMNF, UKFNoOptimization, UKFIntegral, UKFIntegralRandomShoot, UKFStepwise, UKFStepwiseRandomShoot, EKF };

    [Serializable]
    public class FilterQualityInfo
    {
        public int Count;
        public string FilterName;
        public Matrix<double>[] mError;
        public Matrix<double>[] DError;
        public Matrix<double>[] mxHat;
        public Matrix<double>[] mKHat;

        public FilterQualityInfo(string FilterName, int T, Vector<double> X0Hat, Matrix<double> DX0Hat)
        {
            Count = 0;
            this.FilterName = FilterName;
            mError = Exts.ZerosArrayOfShape(X0Hat.ToColumnMatrix(), T);
            DError = Exts.ZerosArrayOfShape(DX0Hat, T);
            mxHat = Exts.ZerosArrayOfShape(X0Hat.ToColumnMatrix(), T);
            mKHat = Exts.ZerosArrayOfShape(DX0Hat, T);
        }

        public FilterQualityInfo(string FilterName, int T, Matrix<double> X0Hat, Matrix<double> DX0Hat) : this(FilterName, T, X0Hat.Column(0), DX0Hat) { }

    }
}
