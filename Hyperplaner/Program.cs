using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;

namespace Hyperplaner
{
    class Program
    {
        static string wikilinksdatasetpath = @"D:\v-zw\wikilinks\train.txt";
        static string embeddingresultspath = @"D:\v-zw\results\embedding";
        static Freebase KB = null;
        static double[] learningrate = { 0.001, 0.01, 0.05, 0.1 };
        static int[] dimension = { 50, 100 };
        static double[] C = {0.0625, 0.125, 0.25, 0.5, 1.0};
        static int[] batchsize = { 30, 60, 90, 120, 150 };
        static List<Tuple<double, int, double, int>> hyperpara = new List<Tuple<double, int, double, int>>();
        static int maxConcurrent = 24;
        static int curIdx = 0;
        static int conThread = 0;
        static object curIdxLock = new object();
        static object conThreadLock = new object();

        class Trainer
        {
            public void TrainModel()
            {
                int idx = 0;
                lock (curIdxLock)
                {
                    idx = curIdx;
                    curIdx++;
                }
                //Hyperplaner model = new Hyperplaner(idx, KB);
                TransE model = new TransE(idx, KB);
                model.Optimization(hyperpara[idx]);
                model.Output(embeddingresultspath + "_entity_" + idx.ToString(), embeddingresultspath + "_relation_" + idx.ToString());
                Console.WriteLine("Completed training a model with " + idx.ToString() + "-th setting.");
                lock (conThreadLock)
                {
                    conThread--;
                }
            }
        }

        static void Main(string[] args)
        {
            #region
            KB = new Freebase(wikilinksdatasetpath);
            KB.BuildNegative();
            foreach (double alpha in learningrate)
            {
                foreach (int k in dimension)
                {
                    foreach (double c in C)
                    {
                        foreach (int b in batchsize)
                        {
                            hyperpara.Add(new Tuple<double, int, double, int>(alpha, k, c, b));
                        }
                    }
                }
            }
            #endregion

            #region
            curIdx = conThread = 0;
            while (true)
            {
                lock (curIdxLock)
                {
                    if (curIdx == hyperpara.Count)
                    {
                        break;
                    }
                }
                lock (conThreadLock)
                {
                    if (conThread < maxConcurrent)
                    {
                        Trainer ob = new Trainer();
                        Thread processor = new Thread(new ThreadStart(ob.TrainModel));
                        processor.Start();
                        while (!processor.IsAlive) ;
                        conThread++;
                        Thread.Sleep(1);
                    }
                }
            }
            Console.WriteLine("Has Provoked All.");
            #endregion

            Console.ReadLine();
        }
    }
}
