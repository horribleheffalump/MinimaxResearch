using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;

namespace NonlinearSystem
{
    public static class Extensions
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

        public static Vector<double> Stack(this Vector<double> v1, Vector<double> v2)
        {
            Vector<double> result = Vector<double>.Build.Dense(v1.Count + v2.Count);
            for (int i = 0; i < v1.Count; i++)
            {
                result[i] = v1[i];
            }
            for (int i = 0; i < v2.Count; i++)
            {
                result[i+v1.Count] = v2[i];
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

        public static Vector<double> cart2pol(Vector<double> x)
        {
            if (x.Count != 2) throw new ArgumentException();
            return Vector<double>.Build.DenseOfArray(new double[] {Math.Atan2(x[1], x[0]), Math.Sqrt(x[0]* x[0] + x[1] * x[1])});
        }
        public static Vector<double> pol2cart(Vector<double> p)
        {
            if (p.Count != 2) throw new ArgumentException();
            return Vector<double>.Build.DenseOfArray(new double[] { p[1] * Math.Cos(p[0]), p[1] * Math.Sin(p[0]) });
        }

        public static Vector<double> cart2sphere(Vector<double> x)
        {
            if (x.Count != 3) throw new ArgumentException();
            double r = Math.Sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
            double theta = Math.Acos(x[2] / r);
            double phi = Math.Atan2(x[1], x[0]);
            return Vector<double>.Build.DenseOfArray(new double[] {r , theta, phi });
        }
        public static Vector<double> sphere2cart(Vector<double> p)
        {
            if (p.Count != 3) throw new ArgumentException();
            double r = p[0];
            double theta = p[1];
            double phi = p[2];
            return Vector<double>.Build.DenseOfArray(new double[] { r * Math.Sin(theta) * Math.Cos(phi), r * Math.Sin(theta) * Math.Sin(phi), r * Math.Cos(theta)});
        }

    }
}
