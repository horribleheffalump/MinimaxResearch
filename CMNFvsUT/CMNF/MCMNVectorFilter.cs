using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using MathNetExtensions;

namespace CMNF
{
    public class ModifiedCMNFilter
    {
        //Func<int, Vector<double>, Vector<double>> Xi;
        //Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> Zeta;

        //public ModifiedCMNFilter(Func<int, Vector<double>, Vector<double>> xi, Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> zeta)
        //{
        //    Xi = xi;
        //    Zeta = zeta;
        //}

        //public Vector<double> Step(int t, Vector<double> y, Vector<double> xHat_)
        //{



        //    Vector<double> xTilde = FHat[t] * Xi(t, xHat_) + fHat[t];
        //    Vector<double> xHat = xTilde + HHat[t] * Zeta(t, xTilde, y, KTilde[t]) + hHat[t];
        //    return xHat;
        //}


        //public void EstimateParameters(DiscreteVectorModel[] models, Vector<double> xhat0, int T)
        //{
        //    int n = models.Count();

        //    Vector<double>[] xHat = Enumerable.Repeat(xhat0, n).ToArray();

        //    for (int t = 0; t < T; t++)
        //    {
        //        Console.WriteLine($"CMNF estimate parameters: t={t}");
        //        Vector<double>[] x = new Vector<double>[n];
        //        Vector<double>[] y = new Vector<double>[n];
        //        Vector<double>[] xiHat = new Vector<double>[n];
        //        for (int i = 0; i < n; i++)
        //        {
        //            //models[i].Step();
        //            //x[i] = models[i].State;
        //            //y[i] = models[i].Obs;
        //            x[i] = models[i].Trajectory[t][0];
        //            y[i] = models[i].Trajectory[t][1];
        //            xiHat[i] = Xi(t, xHat[i]);
        //        }

        //        Matrix<double> F = Exts.Cov(x, xiHat) * (Exts.Cov(xiHat, xiHat).PseudoInverse());
        //        Vector<double> f = x.Average() - F * xiHat.Average();
        //        Matrix<double> kTilde = Exts.Cov(x, x) - Exts.Cov(x, xiHat) * F.Transpose();

        //        Vector<double>[] xTilde = new Vector<double>[n];
        //        Vector<double>[] zetaTilde = new Vector<double>[n];
        //        for (int i = 0; i < n; i++)
        //        {
        //            xTilde[i] = F * xiHat[i] + f;
        //            zetaTilde[i] = Zeta(t, xTilde[i], y[i], kTilde);
        //        }

        //        Matrix<double> H = Exts.Cov(x.Subtract(xTilde), zetaTilde) * (Exts.Cov(zetaTilde, zetaTilde).PseudoInverse());
        //        Vector<double> h = -H * zetaTilde.Average();

        //        Matrix<double> kHat = kTilde - Exts.Cov(x.Subtract(xTilde), zetaTilde) * H.Transpose();
        //        for (int i = 0; i < n; i++)
        //        {
        //            //xHat[i] = F* xiHat[i] +f + H * zetaTilde[i] + h;
        //            xHat[i] = xTilde[i] + H * zetaTilde[i] + h;
        //        }
        //        FHat.Add(t, F);
        //        fHat.Add(t, f);
        //        HHat.Add(t, H);
        //        hHat.Add(t, h);


        //        KTilde.Add(t, kTilde);
        //        KHat.Add(t, kHat);
        //        //KHat.Add(t, Exts.Cov(x, x) - Exts.Cov(x, xiHat) * F - Exts.Cov(x.Subtract(xTilde), zetaTilde) * H);
        //        // KHat.Add(t, cov(x, x) - cov(x, xiHat) * F - cov(x - xTilde, zetaTilde) * H);

        //    }

        //}
    }


}
