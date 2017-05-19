using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;

namespace UKF
{
    public class AsyncCalculator
    {
        private ManualResetEvent doneEvent;
        private double[] result;
        private Func<double[]> calculate;


        public double[] Result { get { return result; } }

        // Constructor.
        public AsyncCalculator(int n, ManualResetEvent _doneEvent, Func<double[]> _calculate)
        {
            doneEvent = _doneEvent;
            calculate = _calculate;
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

        // Recursive method that calculates the Nth Fibonacci number.
    }

    public class AsyncCalculatorPlanner
    {
        private int samplesCount;
        private int packCount;
        private Func<double[]> calculate;

        public AsyncCalculatorPlanner(int _samplesCount, int _packCount, Func<double[]> _calculate)
        {
            samplesCount = _samplesCount;
            packCount = _packCount;
            calculate = _calculate;
        }

        public List<double[]> DoCalculate()
        {
            List<double[]> result = new List<double[]>();
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
