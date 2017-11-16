using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace MathNetExtensions
{
    public static class RandomOptimizer
    {
        const int PackSize = 10;
        public static (double min, Vector<double> argmin) Minimize(Func<Vector<double>, double> Objective, Vector<double> LowerBound, Vector<double> UpperBound, int PointsUniform = 100, int PointsNormal = 100, string OutputFileName = null)
        {
            int n = LowerBound.Count;
            if (PointsUniform <= 0 || PointsNormal < 0)
                throw new ArgumentException("Number of points on the first step must be positive. Number of points on the second must be positive or equal to zero.");
            if (LowerBound.Count != UpperBound.Count)
                throw new ArgumentException("Lower and upper bound dimensions must agree");
            for (int i = 0; i < n; i++)
                if (LowerBound[i] >= UpperBound[i])
                    throw new ArgumentException("Upper bound must be greater then the lower");


            //string filename = string.IsNullOrWhiteSpace(OutputFolder)? null : Path.Combine(OutputFolder, "Random_optimization.txt");

            IContinuousDistribution[] distr = new IContinuousDistribution[n];
            for (int i = 0; i < n; i++)
            {
                distr[i] = new ContinuousUniform(LowerBound[i], UpperBound[i]);
            }

            AsyncCalculatorPlanner acp = new AsyncCalculatorPlanner(PointsUniform, PackSize, () => CalculateSample(Objective, distr));

            List<(double val, Vector<double> x)> results1 = acp.DoCalculate();

            (double val, Vector<double> x) min1 = results1.Where(i => !double.IsNaN(i.val)).OrderBy(i => i.val).First();
            (double val, Vector<double> x) min2;
            if (PointsNormal == 0)
                min2 = min1;
            else
            {
                for (int i = 0; i < n; i++)
                {
                    distr[i] = new Normal(min1.x[i], (UpperBound[i] - LowerBound[i]) / (PointsUniform / n * 3)); // to ajust the standart deviation with average distanse between the uniformly generated points
                }

                acp = new AsyncCalculatorPlanner(PointsNormal, PackSize, () => CalculateSample(Objective, distr));
                List<(double val, Vector<double> x)> results2 = acp.DoCalculate();

                min2 = results2.Where(i => !double.IsNaN(i.val)).OrderBy(i => i.val).First();

                if (!string.IsNullOrWhiteSpace(OutputFileName))
                    using (System.IO.StreamWriter outputfile = new System.IO.StreamWriter(OutputFileName))
                    {
                        NumberFormatInfo provider;
                        provider = new NumberFormatInfo
                        {
                            NumberDecimalSeparator = "."
                        };

                        var results = results1.Concat(results2).Where(i => !double.IsNaN(i.val) && i.val < double.MaxValue);
                        foreach (var e in results.OrderBy(i => i.val))
                        {
                            outputfile.WriteLine(string.Format(provider, "{0}", e.val) + ", " + String.Join(",", e.x.Select(s => string.Format(provider, "{0}", s))));
                        }

                    }
            }
            return min2;
        }

        public static (double, Vector<double>) CalculateSample(Func<Vector<double>, double> Objective, IContinuousDistribution[] distribution)
        {
            int n = distribution.Count();

            Vector<double> x = Exts.Vector(distribution.Select(d => d.Sample()).ToArray());

            double crit = Objective(x);

            return (crit, x);
        }

    }



    public class AsyncCalculator
    {
        private ManualResetEvent doneEvent;
        private (double, Vector<double>) result;
        private Func<(double, Vector<double>)> calculate;


        public (double, Vector<double>) Result { get { return result; } }

        // Constructor.
        public AsyncCalculator(int n, ManualResetEvent doneEvent, Func<(double, Vector<double>)> calculate)
        {
            this.doneEvent = doneEvent;
            this.calculate = calculate;
        }

        // Wrapper method for use with thread pool.
        public void ThreadPoolCallback(Object threadContext)
        {
            int threadIndex = (int)threadContext;
            //Console.WriteLine("thread {0} started...", threadIndex);
            result = calculate();
            //Console.WriteLine("thread {0} result calculated...", threadIndex);
            doneEvent.Set();
        }
        
    }

    public class AsyncCalculatorPlanner
    {
        private int samplesCount;
        private int packCount;
        private Func<(double, Vector<double>)> calculate;

        public AsyncCalculatorPlanner(int samplesCount, int packCount, Func<(double, Vector<double>)> calculate)
        {
            this.samplesCount = samplesCount;
            this.packCount = packCount;
            this.calculate = calculate;
        }

        public List<(double, Vector<double>)> DoCalculate()
        {
            List<(double, Vector<double>)> result = new List<(double, Vector<double>)>();
            for (int pack = 0; pack <= samplesCount / packCount; pack++)
            {
                ManualResetEvent[] doneEvents = new ManualResetEvent[Math.Min(packCount, samplesCount - pack * packCount)];
                AsyncCalculator[] calcArray = new AsyncCalculator[Math.Min(packCount, samplesCount - pack * packCount)];

                // Configure and start threads using ThreadPool.
                //Console.WriteLine("launching {0} tasks...", packCount);
                for (int i = 0; i < Math.Min(packCount, samplesCount - pack * packCount); i++)
                {
                    doneEvents[i] = new ManualResetEvent(false);
                    AsyncCalculator calc = new AsyncCalculator(pack * packCount + i, doneEvents[i], calculate);
                    calcArray[i] = calc;
                    ThreadPool.QueueUserWorkItem(calc.ThreadPoolCallback, i);
                }

                // Wait for all threads in pool to calculate.
                if (doneEvents.Length > 0)
                    WaitHandle.WaitAll(doneEvents);
                //Console.WriteLine("All calculations are complete.");

                // Display the results.
                for (int i = 0; i < Math.Min(packCount, samplesCount - pack * packCount); i++)
                {
                    AsyncCalculator calc = calcArray[i];
                    result.Add(calc.Result);
                    //Console.WriteLine("({0}) = {1}", j.N, j.JOfN);
                }
            }
            return result;
        }
    }
}
