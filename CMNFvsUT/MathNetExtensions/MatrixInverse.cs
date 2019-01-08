using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathNetExtensions
{
    public static partial class Exts
    {
        /// <summary>
        /// Inverts a matrix using Jacobi method
        /// </summary>
        /// <param name="x">Square matrix to invert</param>
        /// <param name="tol1">upper bound of off - diagonal elements (default 1e-32)</param>
        /// <param name="tol2">lower bound of diagonal non-zero elements (default 1e-32)</param>
        /// <returns></returns>
        public static Matrix<double> Inverse(this Matrix<double> x, double tol1 = 1e-32, double tol2 = 1e-32, int maxiter = 10000)
        {
            int iter = 0;
            int N = x.RowCount;
            if (x.RowCount != x.ColumnCount)
            {
                throw new ArgumentException("Matrix should be square");
            }

            Matrix<double> R = Matrix<double>.Build.DenseIdentity(N);   // Result
            Matrix<double> X = Matrix<double>.Build.DenseIdentity(N);   // Eigenvectors matrix 
            Matrix<double> V = Matrix<double>.Build.DenseIdentity(N);   // Elementary rotation matrix
            Matrix<double> L = x;                                       // Eigenvalue matrix 
            Matrix<double> U = L.StrictlyUpperTriangle().PointwiseAbs();        // upper triangular matrix, containing absolute values of off-diagonal elements

            (double max, (int i, int j)) = U.FindMax();                 // maximal absolute value of off - diagonal elements and its position
            while ((max / L.PointwiseAbs().Enumerate().Max() > tol1) && iter < maxiter) // Main rotation loop
            {
                double phi = Math.Atan2(1.0, (L[j, j] - L[i, i]) / L[i, j] / 2) / 2; // rotation angle
                double c = Math.Cos(phi); 
                double s = Math.Sin(phi);
                V[i, i] = c;
                V[i, j] = s;
                V[j, i] = -s;
                V[j, j] = c; 
                X = X * V;          
                L = V.Transpose() * L * V;   
                V[i, i] = 1.0; 
                V[i, j] = 0.0;
                V[j, i] = 0.0;
                V[j, j] = 1.0; 
                U = L.StrictlyUpperTriangle().PointwiseAbs();
                (max, (i, j)) = U.FindMax();
                iter++;
            }
            if (iter == maxiter)
            {
                Console.WriteLine($"maxiter reached in Matrix.Inverse() with tolerance {max / L.PointwiseAbs().Enumerate().Max()}");
            }
            L = Exts.Diag(L.Diagonal().ToArray());            // annigilation of off-diagonal elements
            for (int k = 0; k < N; k++)
            {
                if (L[k, k] > tol2)
                {
                    L[k, k] = 1.0 / L[k, k];
                }
                else
                {
                    L[k, k] = 0.0;
                }
            }
            R = X * L * X.Transpose();
            return R;
        }
        public static (double, (int, int)) FindMax(this Matrix<double> x)
        {
            (int, int) loc = (int.MinValue, int.MinValue);
            double max = double.MinValue;
            for (int i = 0; i < x.RowCount; i++)
            {
                for (int j = 0; j < x.ColumnCount; j++)
                {
                    if (max < x[i, j])
                    {
                        max = x[i, j];
                        loc = (i, j);
                    }
                }
            }
            return (max, loc);
        }

        public static (double, (int, int)) FindMin(this Matrix<double> x)
        {
            (int, int) loc = (int.MinValue, int.MinValue);
            double min = double.MaxValue;
            for (int i = 0; i < x.RowCount; i++)
            {
                for (int j = 0; j < x.ColumnCount; j++)
                {
                    if (min > x[i, j])
                    {
                        min = x[i, j];
                        loc = (i, j);
                    }
                }
            }
            return (min, loc);
        }

    }
}
