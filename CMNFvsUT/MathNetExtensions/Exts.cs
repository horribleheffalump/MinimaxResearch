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
        public static T[] GetRow<T>(this T[,] x, int row)
        {
            var width = x.GetLength(0);
            var height = x.GetLength(1);

            if (row >= height)
                throw new IndexOutOfRangeException("Row Index Out of Range");
            // Ensures the row requested is within the range of the 2-d array


            var returnRow = new T[width];
            for (var i = 0; i < width; i++)
                returnRow[i] = x[i, row];

            return returnRow;
        }
        public static T[] GetCol<T>(this T[,] x, int col)
        {
            var width = x.GetLength(0);
            var height = x.GetLength(1);

            if (col >= width)
                throw new IndexOutOfRangeException("Column Index Out of Range");
            // Ensures the column requested is within the range of the 2-d array


            var returnCol = new T[height];
            for (var j = 0; j < height; j++)
                returnCol[j] = x[col, j];

            return returnCol;
        }


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

        public static Matrix<double>[] Average(this Matrix<double>[,] x, int axis = 0)
        {
            var width = x.GetLength(0);
            var height = x.GetLength(1);

            Matrix<double>[] mx = new Matrix<double>[x.GetLength(axis)];

            if (axis == 0)
            {
                mx = x.GetRow(0);
                for (int i = 1; i < height; i++)
                {
                    Matrix<double>[] r = x.GetRow(i);
                    for (int j = 1; j < width; j++)
                    {
                        mx[j] = mx[j] + r[j];
                    }
                }
                for (int j = 1; j < width; j++)
                {
                    mx[j] = mx[j] / height;
                }
            }
            else
            {
                mx = x.GetCol(0);
                for (int j = 1; j < width; j++)
                {
                    Matrix<double>[] c = x.GetCol(j);
                    for (int i = 1; i < height; i++)
                    {
                        mx[i] = mx[i] + c[i];
                    }
                }
                for (int i = 1; i < height; i++)
                {
                    mx[i] = mx[i] / width;
                }
            }
            return mx;
        }

        public static Matrix<double>[] Average(this Matrix<double>[][] x, double[] weights = null)
        {
            double[] coeffs = Enumerable.Repeat(1.0 / x.Length, x.Length).ToArray();
            if (weights != null)
            {
                double sum = weights.Sum();
                coeffs = weights.Select(e => e / sum).ToArray();
            }
            Matrix<double>[] mx = new Matrix<double>[x[0].Length];
            for (int i = 0; i < x[0].Length; i++)
            {
                mx[i] = ZeroOfShape(x[0][i]);
                for (int j = 1; j < x.Length; j++)
                {
                    mx[i] = mx[i] + x[j][i] * coeffs[j];
                }
            }
            return mx;
        }

        public static Vector<double>[] Average(this Vector<double>[,] x, int axis = 0) // TODO: the same code for Matrix and Vector - NOT COOL!!!
        {
            var width = x.GetLength(0);
            var height = x.GetLength(1);

            Vector<double>[] mx = new Vector<double>[x.GetLength(axis)];

            if (axis == 0)
            {
                mx = x.GetRow(0);
                for (int i = 1; i < height; i++)
                {
                    Vector<double>[] r = x.GetRow(i);
                    for (int j = 1; j < width; j++)
                    {
                        mx[j] = mx[j] + r[j];
                    }
                }
                for (int j = 1; j < width; j++)
                {
                    mx[j] = mx[j] / height;
                }
            }
            else
            {
                mx = x.GetCol(0);
                for (int j = 1; j < width; j++)
                {
                    Vector<double>[] c = x.GetCol(j);
                    for (int i = 1; i < height; i++)
                    {
                        mx[i] = mx[i] + c[i];
                    }
                }
                for (int i = 1; i < height; i++)
                {
                    mx[i] = mx[i] / width;
                }
            }
            return mx;
        }

        public static Vector<double>[] Average(this Vector<double>[][] x, double[] weights = null)
        {
            double[] coeffs = Enumerable.Repeat(1.0 / x.Length, x.Length).ToArray();
            if (weights != null)
            {
                double sum = weights.Sum();
                coeffs = weights.Select(e => e / sum).ToArray();
            }
            Vector<double>[] mx = new Vector<double>[x[0].Length];
            for (int i = 0; i < x[0].Length; i++)
            {
                mx[i] = ZeroOfShape(x[0][i]);
                for (int j = 1; j < x.Length; j++)
                {
                    mx[i] = mx[i] + x[j][i] * coeffs[j];
                }
            }
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
        public static Vector<double>[] ArrayOf(Vector<double> v, int n)
        {
            Vector<double>[] result = new Vector<double>[n];
            for (int i = 0; i < n; i++)
            {
                result[i] = Vector<double>.Build.DenseOfVector(v);
            }
            return result;
        }

        public static Vector<double>[] ZerosArrayOfShape(Vector<double> v, int n)
        {
            return ArrayOf(Vector<double>.Build.Dense(v.Count), n);
        }

        public static Vector<double> ZeroOfShape(Vector<double> v)
        {
            return Vector<double>.Build.Dense(v.Count);
        }

        public static Matrix<double>[] ArrayOf(Matrix<double> m, int n)
        {
            Matrix<double>[] result = new Matrix<double>[n];
            for (int i = 0; i < n; i++)
            {
                result[i] = Matrix<double>.Build.DenseOfMatrix(m);
            }
            return result;
        }

        public static Matrix<double>[] ZerosArrayOfShape(Matrix<double> m, int n)
        {
            return ArrayOf(Matrix<double>.Build.Dense(m.RowCount, m.ColumnCount), n);
        }

        public static Matrix<double> ZeroOfShape(Matrix<double> m)
        {
            return Matrix<double>.Build.Dense(m.RowCount, m.ColumnCount);
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

        public static Vector<double>[] Add(this Vector<double>[] v1, Vector<double>[] v2)
        {
            Vector<double>[] result = new Vector<double>[v1.Length];
            for (int i = 0; i < v1.Length; i++)
            {
                result[i] = v1[i] + v2[i];
            }
            return result;
        }

        public static Vector<double>[] Subtract(this Vector<double>[] v1, Vector<double> v2)
        {
            Vector<double>[] result = new Vector<double>[v1.Length];
            for (int i = 0; i < v1.Length; i++)
            {
                result[i] = v1[i] - v2;
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
