﻿using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UKF;

namespace TestEnvironments
{
    public class UKFWrapper : BasicFilter
    {
        private UKFilter UKF;
        public int T;
        public Vector<double> X0Hat;
        public DiscreteVectorModel[] Models;
        public Func<int, Vector<double>, Vector<double>> Phi1;
        public Func<int, Vector<double>, Matrix<double>> Phi2;
        public Func<int, Vector<double>, Vector<double>> Psi1;
        public Func<int, Vector<double>, Matrix<double>> Psi2;
        public Vector<double> MW;
        public Matrix<double> DW;
        public Vector<double> MNu;
        public Matrix<double> DNu;
        public Func<Matrix<double>, double> Crit;
        public Matrix<double> DX0Hat;
        public string outputFolder;
        public bool doOptimizeWithRandomShoot;
        public bool doCalculateUKFStepwise;

        public override void Initialize()
        {
            if (doOptimizeWithRandomShoot)
                UKF = new UKFilter(UTDefinitionType.ImplicitAlphaBetaKappa, OptimizationMethod.RandomShoot);
            else
                UKF = new UKFilter(UTDefinitionType.ImplicitAlphaBetaKappa, OptimizationMethod.NelderMeed);

            if (doCalculateUKFStepwise)
            {
                Console.WriteLine($"UKF estimate parameters stepwise");
                UKF.EstimateParametersStepwise(Phi1, Phi2, Psi1, Psi2, MW, DW, MNu, DNu, Crit, T, Models, X0Hat, DX0Hat, outputFolder);
            }
            else
            {
                Console.WriteLine($"UKF estimate parameters");
                UKF.EstimateParameters(Phi1, Phi2, Psi1, Psi2, MW, DW, MNu, DNu, Crit, T, Models, X0Hat, DX0Hat, outputFolder);
            }
        }

        public override (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat)
        {
            return UKF.Step(Phi1, Phi2, Psi1, Psi2, MW, DW, MNu, DNu, t, y, xHat, kHat);
        }
    }
}