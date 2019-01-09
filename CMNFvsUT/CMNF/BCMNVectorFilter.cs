using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using MathNetExtensions;
using System.Xml.Serialization;
using System.Xml;
using System.Xml.Schema;

namespace CMNF
{
    public class BCMNFilter
    {
        Func<int, Vector<double>, Vector<double>> Alpha;
        Func<int, Vector<double>, Vector<double>, Vector<double>> Gamma;

        public Dictionary<int, Matrix<double>> FHat;
        public Dictionary<int, Vector<double>> fHat;
        public Dictionary<int, Matrix<double>> HHat;
        public Dictionary<int, Vector<double>> hHat;
        public Dictionary<int, Matrix<double>> GainHat;
        public Dictionary<int, Matrix<double>> KTilde;
        public Dictionary<int, Matrix<double>> KHat;

        public BCMNFilter(Func<int, Vector<double>, Vector<double>> Alpha, Func<int, Vector<double>, Vector<double>, Vector<double>> Gamma)
        {
            this.Alpha = Alpha;
            this.Gamma = Gamma;

            FHat = new Dictionary<int, Matrix<double>>();
            fHat = new Dictionary<int, Vector<double>>();
            HHat = new Dictionary<int, Matrix<double>>();
            hHat = new Dictionary<int, Vector<double>>();
            GainHat = new Dictionary<int, Matrix<double>>();
            KTilde = new Dictionary<int, Matrix<double>>();
            KHat = new Dictionary<int, Matrix<double>>();
        }

        public void EstimateParameters(DiscreteVectorModel[] models, Vector<double> xhat0, int T)
        {

            int n_total = models.Count();

            Vector<double>[] xHat = Enumerable.Repeat(xhat0, n_total).ToArray();
            Console.WriteLine($"BCMNF estimate parameters start");
            DateTime start = DateTime.Now;
            for (int t = 1; t < T; t++) // start from 1 because for 0 we do not have observations
            {
                DateTime startiteration = DateTime.Now;
                Vector<double>[] x = new Vector<double>[n_total];
                Vector<double>[] y = new Vector<double>[n_total];
                Vector<double>[] alpha = new Vector<double>[n_total];
                Vector<double>[] gamma = new Vector<double>[n_total];
                for (int i = 0; i < n_total; i++)
                {
                    x[i] = models[i].Trajectory[t][0];
                    y[i] = models[i].Trajectory[t][1];
                    alpha[i] = Alpha(t, xHat[i]);
                    gamma[i] = Gamma(t, xHat[i], y[i]);
                }

                Vector<double> mX = x.Average();
                Vector<double> mAlpha = alpha.Average();
                Vector<double> mGamma = gamma.Average();

                Matrix<double> covXX = Exts.Cov(x, x);
                Matrix<double> covAlphaAlpha = Exts.Cov(alpha, alpha);
                Matrix<double> covGammaGamma = Exts.Cov(gamma, gamma);

                Matrix<double> covXAlpha = Exts.Cov(x, alpha);
                Matrix<double> covXGamma = Exts.Cov(x, gamma);
                Matrix<double> covGammaAlpha = Exts.Cov(gamma, alpha);

                Matrix<double> invCovAlphaAlpha = covAlphaAlpha.Inverse(1e-32, 1e-32);
                //Matrix<double> invCovAlphaAlpha = covAlphaAlpha.PseudoInverse();


                Matrix<double> F = covXAlpha * invCovAlphaAlpha;
                Vector<double> f = mX - F * mAlpha;
                Matrix<double> kTildeXX = covXX - F * covXAlpha.Transpose();

                Matrix<double> H = covGammaAlpha * invCovAlphaAlpha;
                Vector<double> h = mGamma - H * mAlpha;

                Matrix<double> kTildeXGamma = covXGamma - F * covGammaAlpha.Transpose();
                Matrix<double> kTildeGammaGamma = covGammaGamma - H * covGammaAlpha.Transpose();

                Matrix<double> invKTildeGammaGamma = kTildeGammaGamma.Inverse(1e-32, 1e-32);
                //Matrix<double> invKTildeGammaGamma = kTildeGammaGamma.PseudoInverse();

                Matrix<double> Gain = kTildeXGamma * invKTildeGammaGamma;

                Matrix<double> kHat = kTildeXX - Gain * kTildeXGamma.Transpose();

                for (int i = 0; i < n_total; i++)
                {
                    xHat[i] = F * alpha[i] + f + Gain * (gamma[i] - H * alpha[i] - h);
                }
                FHat.Add(t, F);
                fHat.Add(t, f);
                HHat.Add(t, H);
                hHat.Add(t, h);
                GainHat.Add(t, Gain);


                KTilde.Add(t, kTildeXX);
                KHat.Add(t, kHat);
                Console.WriteLine($"BCMNF estimate parameters for t={t}, done in {(DateTime.Now - startiteration).ToString(@"hh\:mm\:ss\.fff")}");
                x = null;
                y = null;
                alpha = null;
                gamma = null;
            }
            DateTime finish = DateTime.Now;
            Console.WriteLine($"BCMNF estimate parameters finished in {(finish - start).ToString(@"hh\:mm\:ss\.fff")}");


        }

        public (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat_)
        {
            Vector<double> alpha = Alpha(t, xHat_);
            Vector<double> xTilde = FHat[t] * alpha + fHat[t];
            Vector<double> xHat = xTilde + GainHat[t] *(Gamma(t, xHat_, y) - HHat[t] * alpha - hHat[t]);
            return (xHat, KHat[t]);
        }
    }


}
