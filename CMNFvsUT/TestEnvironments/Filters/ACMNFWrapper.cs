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
    public class ACMNFWrapper: BasicFilter
    {
        private AdaptiveCMNFilter ACMNF;
        public int T;
        public Vector<double> X0Hat;
        public DiscreteVectorModel[] Models;
        public Func<int, Vector<double>, Vector<double>> Xi;
        public Func<int, Vector<double>, Vector<double>, Matrix<double>, Vector<double>> Zeta;

        public Func<int, Vector<double>, Vector<double>> Phi1;
        public Func<int, Vector<double>, Matrix<double>> Phi2;
        public Func<int, Vector<double>, Vector<double>> Psi1;
        public Func<int, Vector<double>, Matrix<double>> Psi2;
        public Func<int, Vector<double>> W;
        public Func<int, Vector<double>> Nu;

        public override void Initialize()
        {
            ACMNF = new AdaptiveCMNFilter(Xi, Zeta, Phi1, Phi2, Psi1, Psi2, W, Nu);
        }

        public override void InitializeAndTrain()
        {
            Initialize();
            ACMNF.EstimateParameters(Models, X0Hat, T);
        }

        public override (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat)
        {
            return ACMNF.Step(t, y, xHat, kHat);
        }

        public CMNVectorFilterParams GetParams()
        {
            CMNVectorFilterParams p = new CMNVectorFilterParams();
            p.FHat = ACMNF.FHat.Select(x => x.Value).ToArray();
            p.fHat = ACMNF.fHat.Select(x => x.Value.ToColumnMatrix()).ToArray();
            p.HHat = ACMNF.HHat.Select(x => x.Value).ToArray();
            p.hHat = ACMNF.hHat.Select(x => x.Value.ToColumnMatrix()).ToArray();
            p.KTilde = ACMNF.KTilde.Select(x => x.Value).ToArray();
            p.KHat = ACMNF.KHat.Select(x => x.Value).ToArray();
            return p;
        }

        public void SetParams(CMNVectorFilterParams p)
        {
            for (int t = 0; t < p.FHat.Length; t++)
            {
                // t+1 since we start filtring from t = 1
                ACMNF.FHat.Add(t + 1, p.FHat[t]);
                ACMNF.fHat.Add(t + 1, p.fHat[t].Column(0));
                ACMNF.HHat.Add(t + 1, p.HHat[t]);
                ACMNF.hHat.Add(t + 1, p.hHat[t].Column(0));
                ACMNF.KTilde.Add(t + 1, p.KTilde[t]);
                ACMNF.KHat.Add(t + 1, p.KHat[t]);
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
}
