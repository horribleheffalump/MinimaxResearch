using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;

namespace PythonInteract
{
    /// <summary>
    /// Simple helper class for running Python scripts
    /// </summary>
    public static class Python
    {
        public static string RunScript(string scriptPath, string[] args, string pythonPath = @"C:\Program Files (x86)\Microsoft Visual Studio\Shared\Anaconda3_64\python.exe")
        {
            string output = string.Empty;
            try
            {
                ProcessStartInfo pythonProcessStartInfo = new ProcessStartInfo(pythonPath)
                {
                    UseShellExecute = false,
                    RedirectStandardOutput = true,
                    Arguments = scriptPath + " " + string.Join(" ", args)
                };

                Process pythonProcess = new Process
                {
                    StartInfo = pythonProcessStartInfo
                };
                pythonProcess.Start();

                StreamReader outputStreamReader = pythonProcess.StandardOutput;
                output = outputStreamReader.ReadLine();

                pythonProcess.WaitForExit();
                pythonProcess.Close();
            }
            catch (Exception e)
            {
                output = e.Message;
            }
            return output;
        }
    }
}
