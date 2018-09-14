using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using PythonInteract;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CMNFTest
{
    class TestPolar: TestEnvironmentStatic
    {
        //Vector<double> mX = Vector(30, 40); Matrix<double> KX = Diag(30 * 30, 30 * 30);
        //Vector<double> mNu = Vector(0, 0); Matrix<double> KNu = Diag(Math.Pow(5 * Math.PI / 180.0, 2.0), 30 * 30);
        //Normal[] NormalX = new Normal[2] { new Normal(mX[0], Math.Sqrt(KX[0, 0])), new Normal(mX[1], Math.Sqrt(KX[1, 1])) };
        //Normal[] NormalNu = new Normal[2] { new Normal(mNu[0], Math.Sqrt(KNu[0, 0])), new Normal(mNu[1], Math.Sqrt(KNu[1, 1])) }; ;

        ////Console.WriteLine(mX.ToLine());

        //TestEnvironmentStatic testPolar = new TestEnvironmentStatic
        //{
        //    Phi = x => Extensions.cart2pol(x),
        //    InvPhi = y => Extensions.pol2cart(y),
        //    W = () => Vector(NormalX[0].Sample(), NormalX[1].Sample()),
        //    Nu = () => Vector(NormalNu[0].Sample(), NormalNu[1].Sample()),
        //    MX = mX,
        //    KX = KX,
        //    KNu = KNu,
        //    utOptimizationType = UTOptimizationType.ImplicitAlphaBetaKappa
        //};

    }
}
