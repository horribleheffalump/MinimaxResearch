﻿using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using System.Text;
using System.Threading.Tasks;
using EKF;

namespace TestEnvironments.Filters
{
    internal class EKFWrapper : BasicFilter
    {
        private ExtendedKalmanFilter EKF;
        public int T;
        public Func<int, Vector<double>, Vector<double>> Phi1;
        public Func<int, Vector<double>, Matrix<double>> Phi2;
        public Func<int, Vector<double>, Vector<double>> Psi1;
        public Func<int, Vector<double>, Matrix<double>> Psi2;
        public Func<int, Vector<double>, Matrix<double>> dPhi;
        public Func<int, Vector<double>, Matrix<double>> dPsi;
        public Vector<double> MW;
        public Matrix<double> DW;
        public Vector<double> MNu;
        public Matrix<double> DNu;

        public string outputFolder;
        public Func<int, Vector<double>, Matrix<double>, (Vector<double>, Matrix<double>)> Predict;

        public override void Initialize()
        {
            FilterName = "EKF";
            EKF = new ExtendedKalmanFilter(Phi1, Phi2, Psi1, Psi2, dPhi, dPsi, MW, DW, MNu, DNu, Predict);
        }

        public override void InitializeAndTrain()
        {
            Initialize();
        }

        public override (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat)
        {
            return EKF.Step(t, y, xHat, kHat);
        }
    }
}
