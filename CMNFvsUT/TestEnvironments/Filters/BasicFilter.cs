using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestEnvironments.Filters
{
    /// <summary>
    /// Basic class for filters' standardization
    /// Siutable for filters with recursive structure, 
    /// where the estimate on the current step depend only 
    /// on the current observations, estimate on the previous step,
    /// and estimated error covariance
    /// </summary>
    internal abstract class BasicFilter
    {
        public string FilterName; 
        public string FileName; // file name to save/load filter parameters
        
        /// <summary>
        /// Abstract method for filter initialization
        /// </summary>
        public abstract void Initialize();

        /// <summary>
        /// Abstract method for filter initialization and parameter calculation
        /// </summary>
        public abstract void InitializeAndTrain();

        /// <summary>
        /// Abstract method for estimate calculation
        /// </summary>
        /// <param name="t">Current time instant</param>
        /// <param name="y">Current observation vector</param>
        /// <param name="xHat">Estimate on the previous step</param>
        /// <param name="kHat">Error covariance matrix</param>
        /// <returns>Current state estimate</returns>
        public abstract (Vector<double>, Matrix<double>) Step(int t, Vector<double> y, Vector<double> xHat, Matrix<double> kHat);

        /// <summary>
        /// Virtual method to save parameters as text
        /// </summary>
        public virtual void SaveParamsText() { }

        /// <summary>
        /// Virtual method to save filter parameters as binary
        /// for future use
        /// </summary>
        public virtual void SaveParams() { }

        /// <summary>
        /// Virtual method to load filter parameters 
        /// </summary>
        public virtual void LoadParams() { }
    }
}
