using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

using Accord.Math;

namespace Hyperplaner
{
    class Hyperplaner
    {
        public int mid = -1;
        public Freebase KB = null;
        //public PersonBase KB = null;
        public List<double[]> eeb = new List<double[]>();
        public List<double[]> web = new List<double[]>();
        public List<double> deb = new List<double>();
        public int[] dist = null;
        public List<double> energylist = new List<double>();

        public Hyperplaner(int modelid, Freebase mikeminz)
        {
            double[] test = Matrix.Vector(10, 0.2);
            test.Subtract(Matrix.Vector(10, 0.3), true);
            foreach (double ele in test)
            {
                Console.Write(ele.ToString() + "\t");
            }
            Console.WriteLine();

            mid = modelid;
            KB = mikeminz;
        }

        /* public Hyperplaner(int modelid, PersonBase people)
        {
            mid = modelid;
            KB = people;
        } */

        public void Init(int Dim)
        {
            Random rd = new Random();
            double z = Math.Sqrt((double)Dim);
            for (int e = 0; e < KB.E.Count; e++)
            {
                double[] eb = new double[Dim];
                for (int j = 0; j < eb.Length; j++)
                {
                    eb[j] = (-1.0 + 2.0 * rd.NextDouble()) / z;
                }
                eeb.Add(eb);
            }
            for (int r = 0; r < KB.R.Count; r++)
            {
                double[] eb = new double[Dim];
                for (int j = 0; j < eb.Length; j++)
                {
                    eb[j] = (-1.0 + 2.0 * rd.NextDouble()) / z;
                }
                web.Add(eb);
                deb.Add(-2.0 + 4.0 * rd.NextDouble());
            }
            dist = new int[KB.posTriples.Count];
            for (int i = 0; i < dist.Length; i++) dist[i] = i;
        }

        public void Shuffle(int nbSwap)
        {
            Random rd = new Random();
            for (int i = 0; i < nbSwap; i++)
            {
                int x = rd.Next(dist.Length), y = rd.Next(dist.Length);
                int tmp = dist[x];
                dist[x] = dist[y];
                dist[y] = tmp;
            }
        }

        public double f(double[] w, double[] e, double[] v, double d)
        {
            return w.InnerProduct(e.Subtract(v)) + d;
        }

        public double CalcEnergy(int B, double wt)
        {
            double marginenergy = 0;
            for (int i = 0; i < KB.posTriples.Count; i++)
            {
                double res = f(web[KB.posTriples[i].Item2], eeb[KB.posTriples[i].Item1], eeb[KB.posTriples[i].Item3], deb[KB.posTriples[i].Item2]);
                if (res < 1.0)
                {
                    marginenergy += 1.0 - res;
                }
                res = f(web[KB.negTriples[i].Item2], eeb[KB.negTriples[i].Item1], eeb[KB.negTriples[i].Item3], deb[KB.negTriples[i].Item2]);
                if (res > -1.0)
                {
                    marginenergy += 1.0 + res;
                }
            }
            double enormenergy = 0;
            for (int e = 0; e < KB.E.Count; e++)
            {
                enormenergy += eeb[e].InnerProduct(eeb[e]);
            }
            enormenergy /= B;
            enormenergy /= KB.E.Count;
            enormenergy *= wt;
            enormenergy *= 2;
            double rnormenergy = 0;
            for (int r = 0; r < KB.R.Count; r++)
            {
                rnormenergy += web[r].InnerProduct(web[r]);
                rnormenergy += deb[r] * deb[r];
            }
            rnormenergy /= B;
            rnormenergy /= KB.R.Count;
            rnormenergy *= wt;
            return marginenergy + enormenergy + rnormenergy;
        }

        public void Round(double alpha, int Dim, double wt, int B)
        {
            double[] scalv = Matrix.Vector(Dim, alpha);
            double[] scalc = Matrix.Vector(Dim, 2.0 * wt);
            for (int i = 0; i < dist.Length; i += B)
            {
                Dictionary<int, double[]> egrad = new Dictionary<int, double[]>();
                Dictionary<int, double[]> wgrad = new Dictionary<int, double[]>();
                Dictionary<int, double> dgrad = new Dictionary<int, double>();
                #region
                for (int j = 0; j < B && i + j < dist.Length; j++)
                {
                    Tuple<int, int, int> post = KB.posTriples[i + j];
                    if (egrad.ContainsKey(post.Item1) == false)
                    {
                        egrad.Add(post.Item1, Matrix.Vector(Dim, 0.0));
                    }
                    if (wgrad.ContainsKey(post.Item2) == false)
                    {
                        wgrad.Add(post.Item2, Matrix.Vector(Dim, 0.0));
                        dgrad.Add(post.Item2, 0.0);
                    }
                    if (egrad.ContainsKey(post.Item3) == false)
                    {
                        egrad.Add(post.Item3, Matrix.Vector(Dim, 0.0));
                    }
                    double res = f(web[post.Item2], eeb[post.Item1], eeb[post.Item3], deb[post.Item2]);
                    if (res < 1.0)
                    {
                        egrad[post.Item1] = egrad[post.Item1].Subtract(web[post.Item2]);
                        wgrad[post.Item2] = wgrad[post.Item2].Subtract(eeb[post.Item1].Subtract(eeb[post.Item3]));
                        egrad[post.Item3] = egrad[post.Item3].Add(web[post.Item2]);
                        dgrad[post.Item2] -= 1;
                    }
                    Tuple<int, int, int> negt = KB.ConstructNegEg(post);
                    if (egrad.ContainsKey(negt.Item1) == false)
                    {
                        egrad.Add(negt.Item1, Matrix.Vector(Dim, 0.0));
                    }
                    if (wgrad.ContainsKey(negt.Item2) == false)
                    {
                        wgrad.Add(negt.Item2, Matrix.Vector(Dim, 0.0));
                        dgrad.Add(negt.Item2, 0.0);
                    }
                    if (egrad.ContainsKey(negt.Item3) == false)
                    {
                        egrad.Add(negt.Item3, Matrix.Vector(Dim, 0.0));
                    }
                    res = f(web[negt.Item2], eeb[negt.Item1], eeb[negt.Item3], deb[negt.Item2]);
                    if (res > -1.0)
                    {
                        egrad[negt.Item1] = egrad[negt.Item1].Add(web[negt.Item2]);
                        wgrad[negt.Item2] = wgrad[negt.Item2].Add(eeb[negt.Item1].Subtract(eeb[negt.Item3]));
                        egrad[negt.Item3] = egrad[negt.Item3].Subtract(web[negt.Item2]);
                        dgrad[negt.Item2] += 1;
                    }
                }
                #endregion

                List<int> elist = egrad.Keys.ToList();
                List<int> rlist = wgrad.Keys.ToList();
                foreach (int e in elist)
                {
                    egrad[e] = egrad[e].Add(eeb[e].ElementwiseMultiply(scalc));
                }
                foreach (int r in rlist)
                {
                    wgrad[r] = wgrad[r].Add(web[r].ElementwiseMultiply(scalc));
                    dgrad[r] = dgrad[r] += 2.0 * wt * deb[r];
                }
                foreach (int e in elist)
                {
                    eeb[e].Subtract(egrad[e].ElementwiseMultiply(scalv), true);
                }
                foreach (int r in rlist)
                {
                    web[r].Subtract(wgrad[r].ElementwiseMultiply(scalv), true);
                    deb[r] -= alpha * dgrad[r];
                }
            }
        }

        public bool NoDecrease()
        {
            double recent = 0;
            for (int i = energylist.Count - 50, j = 0; j < 50; j++)
                recent += energylist[i + j];
            recent /= 50;
            double longlongago = 0;
            for (int i = energylist.Count - 100, j = 0; j < 50; j++)
                longlongago += energylist[i + j];
            longlongago /= 50;
            if (recent > longlongago - 1000) return true;
            return false;
        }

        //learning rate, dimension, C, batch_size 
        public void Optimization(Tuple<double, int, double, int> hyperpara)
        {
            Init(hyperpara.Item2);
            Shuffle(2 * dist.Length);
            Console.WriteLine("Model " + mid.ToString() + " begins its iterations.");

            //at most 1000 epoches 
            for (int t = 0; t < 1000; t++)
            {
                double energy = CalcEnergy(hyperpara.Item4, hyperpara.Item3);
                Console.WriteLine(t.ToString() + "-th round of model " + mid.ToString() + " has an energy of " + energy.ToString());
                energylist.Add(energy);
                Shuffle(KB.posTriples.Count);
                Round(hyperpara.Item1, hyperpara.Item2, hyperpara.Item3, hyperpara.Item4);
                if (t >= 500 && NoDecrease())
                {
                    break;
                }
            }
        }

        public void Output(string efilepath, string rfilepath)
        {
            TextWriter ops = new StreamWriter(efilepath);
            foreach (string ename in KB.E.Keys)
            {
                ops.Write(ename + "\t");
                StringBuilder vec = new StringBuilder();
                foreach (double cp in eeb[KB.E[ename]])
                {
                    vec.Append(cp.ToString());
                    vec.Append(",");
                }
                ops.WriteLine(vec.ToString());
            }
            ops.Close();
            ops = new StreamWriter(rfilepath);
            foreach (string rname in KB.R.Keys)
            {
                ops.Write(rname + "\t" + deb[KB.R[rname]].ToString() + "\t");
                StringBuilder vec = new StringBuilder();
                foreach (double cp in web[KB.R[rname]])
                {
                    vec.Append(cp.ToString());
                    vec.Append(",");
                }
                ops.WriteLine(vec.ToString());
            }
            ops.Close();
        }
    }
}
