using System;
using System.Collections.Generic;
using System.Globalization;
using MathNet.Numerics.LinearAlgebra;
using System.Text;
using System.Linq;

namespace MathNetExtensions
{
    public static class Exts
    {
        public static string ToLine(this Vector<double> x)
        {
            NumberFormatInfo provider = new NumberFormatInfo()
            {
                NumberDecimalSeparator = ","
            };
            return string.Join("; ", x.ToArray().Select(elem => String.Format(provider, "{0}", elem)));
        }

        public static Vector<double> Average(this Vector<double>[] x)
        {
            Vector<double> mx = x[0];
            for (int i = 1; i < x.Length; i++)
            {
                mx = mx + x[i];
            }
            mx = mx / x.Length;
            return mx;
        }

        public static Matrix<double> Average(this Matrix<double>[] x)
        {
            Matrix<double> mx = x[0];
            for (int i = 1; i < x.Length; i++)
            {
                mx = mx + x[i];
            }
            mx = mx / x.Length;
            return mx;
        }

        public static Vector<double>[] Subtract(this Vector<double>[] v1, Vector<double>[] v2)
        {
            Vector<double>[] result = new Vector<double>[v1.Length];
            for (int i = 0; i < v1.Length; i++)
            {
                result[i] = v1[i] - v2[i];
            }
            return result;
        }

        public static Vector<double> Stack(this Vector<double> v1, Vector<double> v2)
        {
            Vector<double> result = Vector<double>.Build.Dense(v1.Count + v2.Count);
            for (int i = 0; i < v1.Count; i++)
            {
                result[i] = v1[i];
            }
            for (int i = 0; i < v2.Count; i++)
            {
                result[i + v1.Count] = v2[i];
            }
            return result;
        }

        public static Vector<double>[] Stack(this Vector<double>[] v1, Vector<double>[] v2)
        {
            Vector<double>[] result = new Vector<double>[v1.Length];
            for (int i = 0; i < v1.Length; i++)
            {
                result[i] = v1[i].Stack(v2[i]);
            }
            return result;
        }


        public static Matrix<double> Cov(Vector<double>[] x, Vector<double>[] y)
        {
            Vector<double> mx = x.Average();
            Vector<double> my = y.Average();
            //for (int i = 0; i < x.Length; i++)
            //{
            //    mx = mx + x[i];
            //    my = my + y[i];
            //}
            //mx = mx / x.Length;
            //my = my / y.Length;

            Matrix<double> result = (x[0] - mx).ToColumnMatrix() * (y[0] - my).ToRowMatrix();
            for (int i = 1; i < x.Length; i++)
            {
                result = result + (x[i] - mx).ToColumnMatrix() * (y[i] - my).ToRowMatrix();
            }
            return result / (x.Length - 1.0);
        }
        public static Vector<double> Vector(params double[] val)
        {
            return Vector<double>.Build.Dense(val);
        }
        public static Matrix<double> Matrix(double val)
        {
            return Matrix<double>.Build.Dense(1, 1, val);
        }
        public static Matrix<double> Diag(params double[] val)
        {
            return Matrix<double>.Build.DenseDiagonal(val.Length, val.Length, (i) => val[i]);
        }
    }
}
