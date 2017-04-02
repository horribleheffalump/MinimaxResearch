using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;

namespace UKF
{
    public class UKFilter
    {
        private int L; //state dim
        private double Alpha; //spread param
        private double Beta;  //distribution param
        private double Kappa; //secondary scale param 3 - L;
        private double Lambda; //scale param

        private Vector<double> Wm;
        private Vector<double> Wc;

        public Func<Vector<double>, Vector<double>> Phi;
        public Func<Vector<double>, Vector<double>> Psi;
        public Matrix<double> Rw;
        public Matrix<double> Rnu;

        public UKFilter(int l, Func<Vector<double>, Vector<double>> phi, Func<Vector<double>, Vector<double>> psi, Matrix<double> Rw, Matrix<double> Rnu, double alpha = 1e-3, double beta = 2.0, double kappa = double.NaN)
        {
            L = l;
            Phi = phi;
            Psi = psi;
            this.Rw = Rw;
            this.Rnu = Rnu;

            Alpha = alpha;
            Beta = beta;
            if (double.IsNaN(kappa))
            {
                kappa = 3.0 - L;
                //kappa = 0;
            }
            Kappa = kappa;
            Lambda = Math.Pow(alpha, 2.0) * (L + Kappa) - L;

            Wm = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wm[0] = Lambda / (Lambda + L);

            Wc = Vector<double>.Build.Dense(2 * L + 1, 0.5 / (Lambda + L));
            Wc[0] = Lambda / (Lambda + L) + 1.0 - Math.Pow(Alpha, 2.0) + Beta;

        }

        public Matrix<double> GenerateSigmaPoints(Vector<double> x, Matrix<double> P)
        {
            Matrix<double> Sqrt = (Math.Sqrt(Lambda + L)) * P.Cholesky().Factor.Transpose();

            Matrix<double> Xi = x.ToColumnMatrix();
            for (int i = 0; i < L; i++)
            {
                Xi = Xi.Append((x + Sqrt.Column(i)).ToColumnMatrix());
            }
            for (int i = 0; i < L; i++)
            {
                Xi = Xi.Append((x - Sqrt.Column(i)).ToColumnMatrix());
            }

            return Xi;
        }

        public void UT(Func<Vector<double>, Vector<double>> f, Matrix<double> Xi, Matrix<double> R, out Matrix<double> Upsilon, out Vector<double> y, out Matrix<double> Py)
        {
            Upsilon = f(Xi.Column(0)).ToColumnMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {
                Upsilon = Upsilon.Append(f(Xi.Column(i)).ToColumnMatrix());
            }

            y = Wm[0] * Upsilon.Column(0);
            for (int i = 1; i < 2 * L + 1; i++)
            {
                y = y + Wm[i] * Upsilon.Column(i);
            }

            Py = Wc[0] * (Upsilon.Column(0) - y).ToColumnMatrix() * (Upsilon.Column(0) - y).ToRowMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {
                Py = Py + Wc[i] * (Upsilon.Column(i) - y).ToColumnMatrix() * (Upsilon.Column(i) - y).ToRowMatrix();
            }
            Py = Py + R;


        }

        public void Step(Vector<double> y, Vector<double> xHat_, Matrix<double> P_, out Vector<double> xHat, out Matrix<double> PHat)
        {

            Matrix<double> Xi_ = GenerateSigmaPoints(xHat_, P_);

            Matrix<double> XiStar;
            Vector<double> Xtilde;
            Matrix<double> Ptilde;

            UT(Phi, Xi_, Rw, out XiStar, out Xtilde, out Ptilde);

            //Matrix<double> Xiplus = GenerateSigmaPoints(XiStar.Column(0), Rw);
            Matrix<double> Xi = GenerateSigmaPoints(Xtilde, Ptilde);

            Matrix<double> Upsilon;
            Vector<double> Ytilde;
            Matrix<double> PYtilde;

            UT(Psi, Xi, Rnu, out Upsilon, out Ytilde, out PYtilde);

            Matrix<double> PXY = Wc[0] * (Xi.Column(0) - Xtilde).ToColumnMatrix() * (Upsilon.Column(0) - Ytilde).ToRowMatrix();
            for (int i = 1; i < 2 * L + 1; i++)
            {
                PXY = PXY + Wc[i] * (Xi.Column(i) - Xtilde).ToColumnMatrix() * (Upsilon.Column(i) - Ytilde).ToRowMatrix();
            }

            Matrix<double> K = PXY * PYtilde.Inverse();

            xHat = Xtilde + K * (y - Ytilde);
            PHat = Ptilde - K * PYtilde * K.Transpose();
        }


    }
}
