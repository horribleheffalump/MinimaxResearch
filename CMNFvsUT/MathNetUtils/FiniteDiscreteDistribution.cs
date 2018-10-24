using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Random;

namespace MathNetExtensions
{
    public static class FiniteDiscreteDistribution
    {
        private static double _tolerance = 1E-6;
        private static Random random = new SystemRandomSource(RandomSeed.Robust());

        public static int Sample(Vector<double> measure)
        {
            if (Math.Abs(measure.Sum() - 1.0) > _tolerance)
            {
                throw new ArgumentException("Sum of probabilities should be equal to 1");
            }
            int result = int.MinValue;

            double max = double.MinValue;
            double secondmax = double.MinValue;

            for (int i = 0; i < measure.Count; i++)
            {
                if (measure[i] >= secondmax)
                {
                    if (measure[i] >= max)
                    {
                        secondmax = max;
                        max = measure[i];
                    }
                    else
                    {
                        secondmax = measure[i];
                    }
                }
            }
            //double max = measure.OrderByDescending(x => x).ElementAt(0);
            //double secondmax = measure.OrderByDescending(x => x).ElementAt(1);

            if (secondmax / max < _tolerance) // if one value of measure is much bigger then the others, then the correxponding number is declared sample without drawing 
            {
                result = measure.MaximumIndex();
            }
            else // do the drawing
            {
                Vector<double> intervals = Vector<double>.Build.Dense(measure.ToArray());
                //measure.CopyTo(intervals);
                for (int i = 1; i < intervals.Count; i++)
                {
                    intervals[i] += intervals[i - 1];
                }

                double sample = random.NextDouble();

                for (int i = 0; i < intervals.Count; i++)
                {
                    if (measure[i] < 0)
                        throw new ArgumentException("Probabilities should be positive");
                    if (sample < intervals[i])
                    {
                        result = i;
                        break;
                    }
                }
            }
            return result;
        }
    }
}
