using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using MathNetExtensions;
using MathNet.Numerics.Distributions;
using System.Threading.Tasks;

namespace EKF
{
    public class ExtendedKalmanFilter
    {
        private Func<int, Vector<double>, Vector<double>> Phi1; // Phi1(t, X)
        private Func<int, Vector<double>, Matrix<double>> Phi2;
        private Func<int, Vector<double>, Vector<double>> Psi1;
        private Func<int, Vector<double>, Matrix<double>> Psi2;

        private Func<int, Vector<double>, Matrix<double>> dPhi; // linearized dynamics
        private Func<int, Vector<double>, Matrix<double>> dPsi; // linearized observations

        public Vector<double> MW;
        public Matrix<double> DW;
        public Vector<double> MNu;
        public Matrix<double> DNu;

        // if a specific funcion is provided we use it for the prediction and its covariance calculation. 
        // This is for the continuous dynamics case, where predicted covariance estimate kTilde is calculated as a solution to the Riccati equation
        Func<int, Vector<double>, Matrix<double>, (Vector<double>, Matrix<double>)> Predict = null;
        public ExtendedKalmanFilter(Func<int, Vector<double>, Vector<double>> Phi1,
                                    Func<int, Vector<double>, Matrix<double>> Phi2,
                                    Func<int, Vector<double>, Vector<double>> Psi1,
                                    Func<int, Vector<double>, Matrix<double>> Psi2,
                                    Func<int, Vector<double>, Matrix<double>> dPhi,
                                    Func<int, Vector<double>, Matrix<double>> dPsi,
                                    Vector<double> MW,
                                    Matrix<double> DW,
                                    Vector<double> MNu,
                                    Matrix<double> DNu,
                                    Func<int, Vector<double>, Matrix<double>, (Vector<double>, Matrix<double>)> Predict = null
                                    )
        {
            this.Phi1 = Phi1;
            this.Phi2 = Phi2;
            this.Psi1 = Psi1;
            this.Psi2 = Psi2;
            if (dPhi == null) // if the derivative is not provided, the function is considered linear and its derivative is a constant matrix which is calculated columnwise:  dPhi_k = Phi1(e_k), e_k - unit vector
            {
                this.dPhi = (i, x) =>
                {
                    Vector<double>[] dPhi_vectors = new Vector<double>[MW.Count];
                    for (int k = 0; k < MW.Count; k++)
                    {
                        Vector<double> e_k = Exts.ZeroOfShape(MW);
                        e_k[k] = 1.0;
                        dPhi_vectors[k] = Phi1(i, e_k);
                    }
                    return Matrix<double>.Build.DenseOfColumnVectors(dPhi_vectors);
                };
            }
            else
                this.dPhi = dPhi;
            if (dPsi == null) // if the derivative is not provided, the function is considered linear and its derivative is a constant matrix which is calculated columnwise:  dPsi_k = Psi1(e_k), e_k - unit vector
            {
                this.dPsi = (i, x) =>
                {
                    Vector<double>[] dPsi_vectors = new Vector<double>[MW.Count];
                    for (int k = 0; k < MW.Count; k++)
                    {
                        Vector<double> e_k = Exts.ZeroOfShape(MW);
                        e_k[k] = 1.0;
                        dPsi_vectors[k] = Psi1(i, e_k);
                    }
                    return Matrix<double>.Build.DenseOfColumnVectors(dPsi_vectors);
                };
            }
            else
                this.dPsi = dPsi;
            this.MW = MW;
            this.DW = DW;
            this.MNu = MNu;
            this.DNu = DNu;
            this.Predict = Predict;
        }

        public (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat_, Matrix<double> kHat_)
        {

            // Predict
            Vector<double> xTilde;
            Matrix<double> kTilde;

            if (Predict == null)
            {
                Matrix<double> F = dPhi(t, xHat_);
                Matrix<double> Q = Phi2(t, xHat_) * DW * Phi2(t, xHat_).Transpose();
                xTilde = Phi1(t, xHat_) + Phi2(t, xHat_) * MW;
                kTilde = F * kHat_ * F.Transpose() + Q; 
            }
            else
                (xTilde, kTilde) = Predict(t, xHat_, kHat_); 
                // if special function was provided. This is for the continuous dynamics case, where kTilde is calculated as a solution to the Riccati equation

            // Update
            Matrix<double> H = dPsi(t, xTilde);
            Matrix<double> R = Psi2(t, xTilde) * DNu * Psi2(t, xTilde).Transpose();
            Matrix<double> I = Matrix<double>.Build.DenseIdentity(xTilde.Count);
            Matrix<double> K = kTilde * H.Transpose() * (H * kTilde * H.Transpose() + R).Inverse(1e-10, 1e-10);
            Vector<double> xHat__ = xTilde + K * (y - Psi1(t, xTilde) - Psi2(t, xTilde) * MNu);
            Matrix<double> kHat = (I - K * H) * kTilde;
            return (xHat__, kHat);
        }
    }
}
