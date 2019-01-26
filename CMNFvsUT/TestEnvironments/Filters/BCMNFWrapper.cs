using MathNet.Numerics.LinearAlgebra;
using NonlinearSystem;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CMNF;
using System.Xml.Serialization;
using System.IO;
using MathNetExtensions;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;

namespace TestEnvironments.Filters
{
    internal class BCMNFWrapper: BasicFilter
    {
        private BCMNFilter BCMNF;
        public int T;
        public Vector<double> X0Hat;
        public DiscreteVectorModel[] Models;
        public Func<int, Vector<double>, Vector<double>> Alpha;
        public Func<int, Vector<double>, Vector<double>, Vector<double>> Gamma;

        public override void Initialize()
        {
            FilterName = "BCMNF";
            BCMNF = new BCMNFilter(Alpha, Gamma);
        }

        public override void InitializeAndTrain()
        {
            Initialize();
            BCMNF.EstimateParameters(Models, X0Hat, T);
        }

        public override (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat)
        {
            return BCMNF.Step(t, y, xHat);
        }

        public BCMNVectorFilterParams GetParams()
        {
            BCMNVectorFilterParams p = new BCMNVectorFilterParams();
            p.FHat = BCMNF.FHat.Select(x => x.Value).ToArray();
            p.fHat = BCMNF.fHat.Select(x => x.Value.ToColumnMatrix()).ToArray();
            p.HHat = BCMNF.HHat.Select(x => x.Value).ToArray();
            p.hHat = BCMNF.hHat.Select(x => x.Value.ToColumnMatrix()).ToArray();
            p.GainHat = BCMNF.GainHat.Select(x => x.Value).ToArray();
            p.KTilde = BCMNF.KTilde.Select(x => x.Value).ToArray();
            p.KHat = BCMNF.KHat.Select(x => x.Value).ToArray();
            return p;
        }

        public void SetParams(BCMNVectorFilterParams p)
        {
            for (int t = 0; t < p.FHat.Length; t++)
            {
                // t+1 since we start filtring from t = 1
                BCMNF.FHat.Add(t + 1, p.FHat[t]);
                BCMNF.fHat.Add(t + 1, p.fHat[t].Column(0));
                BCMNF.HHat.Add(t + 1, p.HHat[t]);
                BCMNF.hHat.Add(t + 1, p.hHat[t].Column(0));
                BCMNF.GainHat.Add(t + 1, p.GainHat[t]);
                BCMNF.KTilde.Add(t + 1, p.KTilde[t]);
                BCMNF.KHat.Add(t + 1, p.KHat[t]);
            }

        }

        public override void SaveParams()
        {
            //XmlSerializer formatter = new XmlSerializer(typeof(CMNVectorFilterParams));
            IFormatter formatter = new BinaryFormatter();
            using (Stream stream = new FileStream(FileName, FileMode.Create))
            {
                formatter.Serialize(stream, GetParams());
                stream.Close();
            }
        }
        public override void LoadParams()
        {
            BCMNVectorFilterParams p;
            IFormatter formatter = new BinaryFormatter();
            using (Stream stream = new FileStream(FileName, FileMode.Open))
            {
                p = (BCMNVectorFilterParams)formatter.Deserialize(stream);
                stream.Close();
            }
            SetParams(p);
        }

    }

    [Serializable]
    public class BCMNVectorFilterParams
    {
        public Matrix<double>[] FHat;
        public Matrix<double>[] fHat;
        public Matrix<double>[] HHat;
        public Matrix<double>[] hHat;
        public Matrix<double>[] GainHat;
        public Matrix<double>[] KTilde;
        public Matrix<double>[] KHat;

        public BCMNVectorFilterParams()
        {
        }
    }
}
