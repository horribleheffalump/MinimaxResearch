using System;
using System.Globalization;
using System.Text;
using MathNet.Numerics.LinearAlgebra;


namespace NonlinearSystem
{
    public static class Utils
    {
        public static Vector<double> cart2pol(Vector<double> x)
        {
            if (x.Count != 2) throw new ArgumentException();
            return Vector<double>.Build.DenseOfArray(new double[] { Math.Atan2(x[1], x[0]), Math.Sqrt(x[0] * x[0] + x[1] * x[1]) });
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
            return Vector<double>.Build.DenseOfArray(new double[] { r, theta, phi });
        }
        public static Vector<double> sphere2cart(Vector<double> p)
        {
            if (p.Count != 3) throw new ArgumentException();
            double r = p[0];
            double theta = p[1];
            double phi = p[2];
            return Vector<double>.Build.DenseOfArray(new double[] { r * Math.Sin(theta) * Math.Cos(phi), r * Math.Sin(theta) * Math.Sin(phi), r * Math.Cos(theta) });
        }

        public static string ToLatex(this Matrix<double> x, string format = null)
        {
            StringBuilder result = new StringBuilder();
            if (x.ColumnCount == 1 && x.RowCount == 1)
                result.Append($"{x[0, 0]} ");
            else
            {
                result.AppendLine(@"\left(\begin{array}{" + new String('c', x.ColumnCount) + @"}");
                for (int i = 0; i < x.RowCount; i++)
                {
                    for (int j = 0; j < x.ColumnCount; j++)
                    {
                        if (string.IsNullOrEmpty(format))
                            result.Append($"{x[i, j]} ");
                        else
                            result.Append(x[i, j].ToString(format));
                        if (j < x.ColumnCount - 1) result.Append("& ");
                    }
                    result.AppendLine(@"\\");
                }
                result.AppendLine(@"\end{array}\right)");
            }
            return result.ToString();
        }

        public static string ToMatlab(this Matrix<double> x)
        {
            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";
            StringBuilder result = new StringBuilder();
            result.Append("[");
            for (int i = 0; i < x.RowCount; i++)
            {
                for (int j = 0; j < x.ColumnCount; j++)
                {
                    result.Append(string.Format(provider, "{0}", x[i,j]));
                    if (j < x.ColumnCount - 1) result.Append(", ");
                }
                if (i < x.RowCount - 1) result.Append("; ");
            }
            result.AppendLine("]");
            return result.ToString();
        }


        public static string ToLatex(this Vector<double> x)
        {
            return x.ToColumnMatrix().ToLatex();
        }

        public static string ToLatex(this string[] x)
        {
            StringBuilder result = new StringBuilder();
            if (x.Length == 1)
                result.Append($"{x[0]} ");
            else
            {
                result.AppendLine(@"\left(\begin{array}{c}");
                for (int i = 0; i < x.Length; i++)
                {
                    result.Append($"{x[i]} ");
                    result.AppendLine(@"\\");
                }
                result.AppendLine(@"\end{array}\right)");
            }
            return result.ToString();
        }

        public static string ToLatex(this string[][] x)
        {
            StringBuilder result = new StringBuilder();
            if (x.Length == 1)
                result.Append($"{x[0][0]} ");
            else
            {
                result.AppendLine(@"\left(\begin{array}{" + new String('c', x.Length) + @"}");
                for (int i = 0; i < x.Length; i++)
                {
                    for (int j = 0; j < x[i].Length; j++)
                    {
                        result.Append($"{x[i][j]} ");
                        if (j < x[i].Length - 1) result.Append("& ");
                    }
                    result.AppendLine(@"\\");
                }
                result.AppendLine(@"\end{array}\right)");
            }
            return result.ToString();
        }

    }
}
