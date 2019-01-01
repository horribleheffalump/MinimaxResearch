using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestEnvironments
{
    [Serializable]
    class SingleTrajectoryInfo
    {
        //public bool Marked = false;
        public int T;
        public Matrix<double>[] x;
        public Matrix<double>[] y;
        public Dictionary<string, SingleFilterInfo> Filters;

        public SingleTrajectoryInfo(int T, string[] Filters)
        {
            this.T = T;
            x = new Matrix<double>[T];
            y = new Matrix<double>[T];
            this.Filters = new Dictionary<string, SingleFilterInfo>();
            foreach (var f in Filters)
            {
                this.Filters.Add(f, new SingleFilterInfo(T, f));
            }
        }

        public void SaveToText(string FileNameStateTemplate, string FileNameObsTemplate)
        {
            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            string[] filters = Filters.Keys.ToArray();

            int dimX = x[1].RowCount;
            int dimY = y[1].RowCount;

            for (int k = 0; k < dimX; k++)
            {
                using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(FileNameStateTemplate.Replace("{0}", k.ToString())))
                {
                    outputfile.Write(string.Format(provider, "{0} {1}",
                        "t", "x"
                        ));
                    for (int i = 0; i < filters.Length; i++)
                    {
                        outputfile.Write(string.Format(provider, " {0}",
                            filters[i] + "_Error"
                            ));
                    }
                    outputfile.WriteLine();

                    for (int t = 1; t < T; t++) // start from 1 because for 0 we do not have observations
                    {
                        outputfile.Write(string.Format(provider, "{0} {1}",
                            t, x[t][k, 0]
                            ));
                        for (int i = 0; i < filters.Length; i++)
                        {
                            outputfile.Write(string.Format(provider, " {0}",
                            Filters[filters[i]].Err[t][k, 0]
                            ));
                        }
                        outputfile.WriteLine();
                    }
                    outputfile.Close();
                }
            }
            for (int k = 0; k < dimY; k++)
            {
                using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(FileNameObsTemplate.Replace("{0}", k.ToString())))
                {
                    outputfile.Write(string.Format(provider, "{0} {1}",
                         "t", "y"
                    ));
                    outputfile.WriteLine();

                    for (int t = 1; t < T; t++) // start from 1 because for 0 we do not have observations
                    {
                        outputfile.Write(string.Format(provider, "{0} {1}",
                            t, y[t][k, 0]
                            ));
                        outputfile.WriteLine();
                    }
                    outputfile.Close();
                }
            }

        }
    }
}
