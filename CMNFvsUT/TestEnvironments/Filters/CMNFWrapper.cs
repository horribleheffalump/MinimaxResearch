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

namespace TestEnvironments
{
    public class CMNFWrapper: BasicFilter
    {
        private CMNFilter CMNF;
        public int T;
        public Vector<double> X0Hat;
        public DiscreteVectorModel[] Models;
        public Func<int, Vector<double>, Vector<double>> Xi;
        public Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> Zeta;

        public override void Initialize()
        {
            CMNF = new CMNFilter(Xi, Zeta);
        }

        public override void InitializeAndTrain()
        {
            Initialize();
            CMNF.EstimateParameters(Models, X0Hat, T);
        }

        public override (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat)
        {
            return CMNF.Step(t, y, xHat);
        }

        public CMNVectorFilterParams GetParams()
        {
            CMNVectorFilterParams p = new CMNVectorFilterParams();
            p.FHat = CMNF.FHat.Select(x => x.Value).ToArray();
            p.fHat = CMNF.fHat.Select(x => x.Value.ToColumnMatrix()).ToArray();
            p.HHat = CMNF.HHat.Select(x => x.Value).ToArray();
            p.hHat = CMNF.hHat.Select(x => x.Value.ToColumnMatrix()).ToArray();
            p.KTilde = CMNF.KTilde.Select(x => x.Value).ToArray();
            p.KHat = CMNF.KHat.Select(x => x.Value).ToArray();
            return p;
        }

        public void SetParams(CMNVectorFilterParams p)
        {
            for (int t = 0; t < p.FHat.Length; t++)
            {
                CMNF.FHat.Add(t, p.FHat[t]);
                CMNF.fHat.Add(t, p.fHat[t].Column(0));
                CMNF.HHat.Add(t, p.HHat[t]);
                CMNF.hHat.Add(t, p.hHat[t].Column(0));
                CMNF.KTilde.Add(t, p.KTilde[t]);
                CMNF.KHat.Add(t, p.KHat[t]);
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
            CMNVectorFilterParams p;
            IFormatter formatter = new BinaryFormatter();
            using (Stream stream = new FileStream(FileName, FileMode.Open))
            {
                p = (CMNVectorFilterParams)formatter.Deserialize(stream);
                stream.Close();
            }
            SetParams(p);
        }

    }

    [Serializable]
    public class CMNVectorFilterParams
    {
        public Matrix<double>[] FHat;
        public Matrix<double>[] fHat;
        public Matrix<double>[] HHat;
        public Matrix<double>[] hHat;
        public Matrix<double>[] KTilde;
        public Matrix<double>[] KHat;

        public CMNVectorFilterParams()
        {
        }
    }
}
