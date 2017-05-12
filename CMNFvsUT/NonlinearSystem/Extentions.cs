using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;

namespace NonlinearSystem
{
    public static class Extentions
    {
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

    }
}
