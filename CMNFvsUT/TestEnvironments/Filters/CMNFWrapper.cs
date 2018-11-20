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

        public CMFFilterParams GetParams()
        {
            CMFFilterParams p = new CMFFilterParams();
            p.FHat = CMNF.FHat.Select(x => x.Value).ToArray();
            p.fHat = CMNF.fHat.Select(x => x.Value.ToColumnMatrix()).ToArray();
            p.HHat = CMNF.HHat.Select(x => x.Value).ToArray();
            p.hHat = CMNF.hHat.Select(x => x.Value.ToColumnMatrix()).ToArray();
            p.KTilde = CMNF.KTilde.Select(x => x.Value).ToArray();
            p.KHat = CMNF.KHat.Select(x => x.Value).ToArray();
            return p;
        }
        public override void SaveParams()
        {
            //XmlSerializer formatter = new XmlSerializer(typeof(CMFFilterParams));
            IFormatter formatter = new BinaryFormatter();
            using (Stream stream = new FileStream(FileName, FileMode.Create))
            {
                formatter.Serialize(stream, GetParams());
                stream.Close();
            }
        }
        public override void LoadParams()
        {
            CMFFilterParams p;
            IFormatter formatter = new BinaryFormatter();
            using (Stream stream = new FileStream(FileName, FileMode.Open))
            {
                p = (CMFFilterParams)formatter.Deserialize(stream);
                stream.Close();
            }
        }

    }

    [Serializable]
    public class CMFFilterParams: ISerializable
    {
        public Matrix<double>[] FHat;
        public Matrix<double>[] fHat;
        public Matrix<double>[] HHat;
        public Matrix<double>[] hHat;
        public Matrix<double>[] KTilde;
        public Matrix<double>[] KHat;

        public CMFFilterParams()
        {
        }

        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("FHat", FHat);
            info.AddValue("fHat", FHat);
            info.AddValue("HHat", FHat);
            info.AddValue("hHat", FHat);
            info.AddValue("KTilde", FHat);
            info.AddValue("KHat", FHat);
        }

        public CMFFilterParams(SerializationInfo info, StreamingContext context)
        {
            FHat = (Matrix<double>[])info.GetValue("FHat", typeof(Matrix<double>[]));
            fHat = (Matrix<double>[])info.GetValue("FHat", typeof(Matrix<double>[]));
            HHat = (Matrix<double>[])info.GetValue("FHat", typeof(Matrix<double>[]));
            hHat = (Matrix<double>[])info.GetValue("FHat", typeof(Matrix<double>[]));
            KTilde = (Matrix<double>[])info.GetValue("FHat", typeof(Matrix<double>[]));
            KHat = (Matrix<double>[])info.GetValue("FHat", typeof(Matrix<double>[]));
        }
    }
}
