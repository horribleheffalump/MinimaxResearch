using MathNet.Numerics.LinearAlgebra;
using MathNetExtensions;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using System.Text;
using System.Threading.Tasks;

namespace TestEnvironments
{
    [Serializable]
    public class ProcessInfo
    {
        public int Count;
        public Matrix<double>[] mx;
        public Matrix<double>[] Dx;
        public FilterQualityInfo[] FilterQualityInfos;
        public ProcessInfo(int T, string[] Filters, Vector<double> X0Hat, Matrix<double> DX0Hat)
        {
            Count = 0;
            mx = Exts.ZerosArrayOfShape(X0Hat.ToColumnMatrix(), T);
            Dx = Exts.ZerosArrayOfShape(DX0Hat, T);
            FilterQualityInfos = new FilterQualityInfo[Filters.Length];
            for (int j = 0; j < Filters.Length; j++)
            {
                FilterQualityInfos[j] = new FilterQualityInfo(Filters[j], T, X0Hat, DX0Hat);
            }


        }

        public static ProcessInfo LoadFromFile(string fileName)
        {
            ProcessInfo p;
            IFormatter formatter = new BinaryFormatter();
            using (Stream stream = new FileStream(fileName, FileMode.Open))
            {
                p = (ProcessInfo)formatter.Deserialize(stream);
                stream.Close();
            }
            return p;
        }

        public ProcessInfo(ProcessInfo[] infos)
        {
            Count = (int)infos.Select(i => (double)i.Count).Sum();
            double[] TrCounts = infos.Select(i => (double)i.Count).ToArray();
            mx = Exts.Average(infos.Select(i => i.mx).ToArray(), TrCounts);
            Dx = Exts.Average(infos.Select(i => i.Dx).ToArray(), TrCounts);
            FilterQualityInfos = new FilterQualityInfo[infos[0].FilterQualityInfos.Count()];
            for (int j = 0; j < FilterQualityInfos.Count(); j++)
            {
                string filter = infos[0].FilterQualityInfos[j].FilterName;
                FilterQualityInfos[j] = new FilterQualityInfo(filter, infos[0].FilterQualityInfos[j].mError.Length, infos[0].FilterQualityInfos[j].mError[0], infos[0].FilterQualityInfos[j].DError[0]);
                double[] FTrCounts = infos.Select(i => (double)i.FilterQualityInfos[j].Count).ToArray();
                FilterQualityInfos[j].mError = Exts.Average(infos.Select(i => i.FilterQualityInfos[j].mError).ToArray(), FTrCounts);
                FilterQualityInfos[j].DError = Exts.Average(infos.Select(i => i.FilterQualityInfos[j].DError).ToArray(), FTrCounts);
                FilterQualityInfos[j].mxHat = Exts.Average(infos.Select(i => i.FilterQualityInfos[j].mxHat).ToArray(), FTrCounts);
                FilterQualityInfos[j].mKHat = Exts.Average(infos.Select(i => i.FilterQualityInfos[j].mKHat).ToArray(), FTrCounts);
            }
        }

        public void SaveToFile(string fileName)
        {
            IFormatter formatter = new BinaryFormatter();
            using (Stream stream = new FileStream(fileName, FileMode.Create))
            {
                formatter.Serialize(stream, this);
                stream.Close();
            }
        }

        public void SaveToText(string fileName)
        {
            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            for (int k = 0; k < mx[0].RowCount; k++)
            {
                using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString())))
                {
                    outputfile.Write(string.Format(provider, "{0} {1} {2}",
                        "t", "mx", "Dx"
                        ));
                    for (int j = 0; j < FilterQualityInfos.Count(); j++)
                    {
                        outputfile.Write(string.Format(provider, " {0} {1} {2} {3}",
                        FilterQualityInfos[j].FilterName + "_mxHat", FilterQualityInfos[j].FilterName + "_mError", FilterQualityInfos[j].FilterName + "_DError", FilterQualityInfos[j].FilterName + "_mKHat"
                        ));
                    }
                    outputfile.WriteLine();
                    outputfile.Close();
                }
            }
            for (int t = 0; t < mx.Length; t++)
            {
                for (int k = 0; k < mx[0].RowCount; k++)
                {
                    using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(fileName.Replace("{0}", k.ToString()), true))
                    {
                        outputfile.Write(string.Format(provider, "{0} {1} {2}",
                            t, mx[t][k, 0], Dx[t][k, k]
                            ));
                        for (int j = 0; j < FilterQualityInfos.Count(); j++)
                        {
                            outputfile.Write(string.Format(provider, " {0} {1} {2} {3}",
                            FilterQualityInfos[j].mxHat[t][k, 0], FilterQualityInfos[j].mError[t][k, 0], FilterQualityInfos[j].DError[t][k, k], FilterQualityInfos[j].mKHat[t][k, k]
                            ));
                        }
                        outputfile.WriteLine();
                        outputfile.Close();
                    }
                }
            }
        }
    }
}
