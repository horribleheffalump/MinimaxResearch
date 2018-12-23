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

        public ExtendedKalmanFilter(Func<int, Vector<double>, Vector<double>> Phi1,
                                    Func<int, Vector<double>, Matrix<double>> Phi2,
                                    Func<int, Vector<double>, Vector<double>> Psi1,
                                    Func<int, Vector<double>, Matrix<double>> Psi2,
                                    Func<int, Vector<double>, Matrix<double>> dPhi,
                                    Func<int, Vector<double>, Matrix<double>> dPsi,
                                    Vector<double> MW,
                                    Matrix<double> DW,
                                    Vector<double> MNu,
                                    Matrix<double> DNu
                                    )
        {
            this.Phi1 = Phi1;
            this.Phi2 = Phi2;
            this.Psi1 = Psi1;
            this.Psi2 = Psi2;
            this.dPhi = dPhi;
            this.dPsi = dPsi;
            this.MW = MW;
            this.DW = DW;
            this.MNu = MNu;
            this.DNu = DNu;
        }

        public (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat_, Matrix<double> kHat_)
        {
            Matrix<double> F = dPhi(t, xHat_);
            Matrix<double> H = dPsi(t, xHat_);
            Matrix<double> Q = Phi2(t, xHat_) * DW * Phi2(t, xHat_).Transpose();
            Matrix<double> R = Psi2(t, xHat_) * DNu * Psi2(t, xHat_).Transpose();
            Matrix<double> I = Matrix<double>.Build.DenseIdentity(xHat_.Count);

            // Predict
            Vector<double> xTilde = Phi1(t, xHat_) + Phi2(t, xHat_) * MW;
            Matrix<double> KTilde = F * kHat_ * F.Transpose() + Q;
            // Update
            Matrix<double> K = KTilde * H.Transpose() * (H * KTilde * H.Transpose() + R).PseudoInverse();
            Vector<double> xHat__ = xTilde + K * (y - Psi1(t, xHat_) - Psi2(t, xHat_) * MNu);
            Matrix<double> kHat = (I - K * H) * KTilde;
            return (xHat__, kHat);
        }
    }
}
