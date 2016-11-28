using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

using Accord.Math;



namespace Hyperplaner 
{
    class TransE
    {
        public int mid = -1;
        Freebase KB = null;
        public List<double[]> eeb = new List<double[]>();
        public List<double[]> reb = new List<double[]>();
        public int[] dist = null;
        public List<double> energylist = new List<double>();

        /* public List<Tuple<int, int, int>> Triples = null;
        public HashSet<Tuple<int, int, int>> LookupTable = null;
        static public List<Tuple<int, int, int>> heldoutTriples = null;
        static public List<Tuple<int, int, int>> valsetnegs = new List<Tuple<int, int, int>>();
        static public List<Tuple<int, int, int>> heldoutvalsetnegs = new List<Tuple<int, int, int>>();
        static public int firstheldoute = -1;
        static public int firstheldoutr = -1; */



        public TransE(int modelid, Freebase people)
        {
            mid = modelid;
            KB = people;
        }

        /* private double BoxMullerTransform(double u1, double u2)
        {
            return Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
        } */

        private double UniformRange(double l, double r, double pos)
        {
            return l + (r - l) * pos;
        }

        private void Init(int Dim)
        {
            Random rd = new Random();
            double z = Math.Sqrt((double)Dim);
            for (int i = 0; i < KB.E.Count; i++)
            {
                double[] eb = new double[Dim];
                for (int j = 0; j < eb.Length; j++)
                {
                    eb[j] = UniformRange(-1.0, 1.0, rd.NextDouble()) / z;
                }
                eeb.Add(eb);
            }
            for (int i = 0; i < KB.R.Count; i++)
            {
                double[] eb = new double[Dim];
                for (int j = 0; j < eb.Length; j++)
                {
                    eb[j] = UniformRange(-1.0, 1.0, rd.NextDouble()) / z;
                }
                reb.Add(eb);
            }
            dist = new int[KB.posTriples.Count];
            for (int i = 0; i < dist.Length; i++)
            {
                dist[i] = i;
            }
        }

        private void Shuffle(int nbSwap)
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

        /* private List<Tuple<int, int, int>> GenerateNegEg(int l, int r, List<Tuple<int,int,int>> poslist)
        {
            List<Tuple<int, int, int>> result = new List<Tuple<int, int, int>>();

            for (int i = l; i < r; i++)
            {
                result.Add(ConstructNegativeEg(poslist[i]));
            }

            return result;
        } */

        private void Round(double alpha, double C, int B)
        {
            double[] scalv = Matrix.Vector(eeb[0].Length, 2 * alpha);
            double[] scalc = Matrix.Vector(eeb[0].Length, C);
            for (int i = 0; i < dist.Length; i += B)
            {
                Dictionary<int, double[]> egrad = new Dictionary<int, double[]>();
                Dictionary<int, double[]> rgrad = new Dictionary<int, double[]>();

                for (int j = 0; j < B && i + j < dist.Length; j++)
                {
                    Tuple<int, int, int> en = KB.posTriples[dist[i + j]];
                    Tuple<int, int, int> negen = KB.ConstructNegEg(en);
                    if (egrad.ContainsKey(en.Item1) == false)
                    {
                        egrad.Add(en.Item1, Matrix.Vector(eeb[0].Length, 0.0));
                    }
                    if (egrad.ContainsKey(en.Item3) == false)
                    {
                        egrad.Add(en.Item3, Matrix.Vector(eeb[0].Length, 0.0));
                    }
                    if (egrad.ContainsKey(negen.Item1) == false)
                    {
                        egrad.Add(negen.Item1, Matrix.Vector(eeb[0].Length, 0.0));
                    }
                    if (egrad.ContainsKey(negen.Item3) == false)
                    {
                        egrad.Add(negen.Item3, Matrix.Vector(eeb[0].Length, 0.0));
                    }
                    if (rgrad.ContainsKey(en.Item2) == false)
                    {
                        rgrad.Add(en.Item2, Matrix.Vector(reb[0].Length, 0.0));
                    }
                    if (rgrad.ContainsKey(negen.Item2) == false)
                    {
                        rgrad.Add(negen.Item2, Matrix.Vector(reb[0].Length, 0.0));
                    }
                    double[] Skb = eeb[en.Item1].Add(reb[en.Item2].Subtract(eeb[en.Item3]));
                    double skbNorm = Skb.InnerProduct(Skb);
                    double[] Skbp = eeb[negen.Item1].Add(reb[negen.Item2].Subtract(eeb[negen.Item3]));
                    double skbpNorm = Skbp.InnerProduct(Skbp);
                    if (skbNorm + 1 - skbpNorm > 0)
                    {
                        egrad[en.Item1] = egrad[en.Item1].Add(Skb);
                        rgrad[en.Item2] = rgrad[en.Item2].Add(Skb);
                        egrad[en.Item3] = egrad[en.Item3].Subtract(Skb);
                        egrad[negen.Item1] = egrad[negen.Item1].Subtract(Skbp);
                        rgrad[negen.Item2] = rgrad[negen.Item2].Subtract(Skbp);
                        egrad[negen.Item3] = egrad[negen.Item3].Add(Skbp);
                    }
                }
                //check norm constraints 
                List<int> eidlist = egrad.Keys.ToList();
                List<int> ridlist = rgrad.Keys.ToList();
                foreach (int eid in eidlist)
                {
                    double norm = eeb[eid].InnerProduct(eeb[eid]);
                    if (norm > 1)
                    {
                        egrad[eid] = egrad[eid].Add(eeb[eid].ElementwiseMultiply(scalc));
                    }
                }
                foreach (int rid in ridlist)
                {
                    double norm = reb[rid].InnerProduct(reb[rid]);
                    if (norm > 1)
                    {
                        rgrad[rid] = rgrad[rid].Add(reb[rid].ElementwiseMultiply(scalc));
                    }
                }
                //update
                foreach (int eid in eidlist)
                {
                    eeb[eid].Subtract(egrad[eid].ElementwiseMultiply(scalv), true);
                }
                foreach (int rid in ridlist)
                {
                    reb[rid].Subtract(rgrad[rid].ElementwiseMultiply(scalv), true);
                }
            }
        }

        /* private double CheckAcc()
        {
            Random rd = new Random();
            int cor = 0;

            for (int i = 0; i < 1000000; i++)
            {
                int idx = rd.Next(Triples.Count);
                double[] Skb = eeb[Triples[idx].Item1].Add(reb[Triples[idx].Item2].Subtract(eeb[Triples[idx].Item3]));
                double Norm = Skb.InnerProduct(Skb);
                bool flag = true;

                for (int r = 0; r < R.Count; r++)
                {
                    if (r != Triples[idx].Item2)
                    {
                        Tuple<int, int, int> tmp = new Tuple<int, int, int>(Triples[idx].Item1,r,Triples[idx].Item3);
                        if (LookupTable.Contains(tmp)) continue;

                        double[] Skbp = eeb[Triples[idx].Item1].Add(reb[r].Subtract(eeb[Triples[idx].Item3]));
                        double NormP = Skbp.InnerProduct(Skbp);
                        if (NormP < Norm)
                        {
                            flag = false;
                            break;
                        }
                    }
                }
                if (flag) cor++;
            }

            return (double)cor / (double)1000000;
        } */

        /* private double CheckENorm()
        {
            Random rd = new Random();

            int cor = 0;
            for (int i = 0; i < 100000; i++)
            {
                int idx = rd.Next(E.Count);
                double norm = eeb[idx].InnerProduct(eeb[idx]);
                if (norm <= 1) cor++;
            }

            return (double)cor / (double)100000;
        } */

        /* private double CheckRNorm()
        {
            int cor = 0;
            for (int i = 0; i < R.Count; i++)
            {
                double norm = reb[i].InnerProduct(reb[i]);
                if (norm <= 1) cor++;
            }
            return (double)cor / (double)R.Count;
        } */

        /* private void ConstructValidationSet(List<Tuple<int,int,int>> poslist, List<Tuple<int,int,int>> neglist)
        {
            for (int i = 0; i < poslist.Count; i++)
            {
                neglist.Add(ConstructNegativeEg(poslist[i])); 
            }
        } */

        private double CalcEnergy(int B, double C)
        {
            double energy = 0, entitynormenergy = 0, relationnormenergy = 0;
            for (int i = 0; i < KB.posTriples.Count; i++)
            {
                double[] Skb = eeb[KB.posTriples[i].Item1].Add(reb[KB.posTriples[i].Item2].Subtract(eeb[KB.posTriples[i].Item3]));
                double skbNorm = Skb.InnerProduct(Skb);
                double[] Skbp = eeb[KB.negTriples[i].Item1].Add(reb[KB.negTriples[i].Item2].Subtract(eeb[KB.negTriples[i].Item3]));
                double skbpNorm = Skbp.InnerProduct(Skbp);
                if (skbNorm + 1 - skbpNorm > 0)
                {
                    energy += skbNorm + 1 - skbpNorm;
                }
            }
            for (int e = 0; e < KB.E.Count; e++)
            {
                double norm = eeb[e].InnerProduct(eeb[e]);
                if (norm > 1)
                {
                    entitynormenergy += (norm - 1);
                }
            }
            entitynormenergy /= B;
entitynormenergy /= KB.E.Count;
            entitynormenergy *= C;
entitynormenergy *= 4.0;
            
            for (int r = 0; r < KB.R.Count; r++)
            {
                double norm = reb[r].InnerProduct(reb[r]);
                if (norm > 1)
                {
                    energy += (norm - 1);
                }
            }
relationnormenergy /= KB.R.Count;
            relationnormenergy /= B;
            relationnormenergy *= C;
            relationnormenergy *= 2.0;

            return energy + entitynormenergy + relationnormenergy;
        }

        /* public void IncInitialize(List<Tuple<int, int, int>> heldoutKB, Dictionary<string,int> ES, Dictionary<string,int> RS, int ebb, int rbb)
        {
            heldoutTriples = heldoutKB;
            E = ES;
            R = RS;
            firstheldoute = ebb;
            firstheldoutr = rbb;
            foreach (var en in heldoutTriples)
            {
                LookupTable.Add(en);
            }
        } */

        /* public void AdjustTestEg(int Dim, int nbSwap, int T, int B, double alpha)
        {
            TextWriter ops = new StreamWriter(@"D:\v-zw\EnergyRecord_test.dat");
            ConstructValidationSet(heldoutTriples, heldoutvalsetnegs);
            Init(Dim, E.Count - firstheldoute, R.Count - firstheldoutr);
            Shuffle(3 * nbSwap, heldoutTriples, heldoutvalsetnegs);
            Console.WriteLine("Begin Iteration:");

            for (int t = 0; t < T; t++)
            {
                double energy = CalcEnergy(B, heldoutTriples, heldoutvalsetnegs);
                Console.WriteLine("Current energy is " + energy.ToString());
                ops.WriteLine(t.ToString() + "\t" + energy.ToString());

                Shuffle(nbSwap, heldoutTriples, heldoutvalsetnegs);
                Round(B, alpha, heldoutTriples);
                Console.WriteLine("Complete Round "+t.ToString());
            }

            ops.Close();
        } */

        private bool NoDecrease()
        {
            double recent = 0;
            for (int i = energylist.Count - 50, j = 0; j < 50; j++)
            {
                recent += energylist[i + j];
            }
            recent /= 50;
            double longlongago = 0;
            for (int i = energylist.Count - 100, j = 0; j < 50; j++)
            {
                longlongago += energylist[i + j];
            }
            longlongago /= 50;
            if (recent > longlongago - 1000) return true;
            return false;
        }

        public void Optimization(Tuple<double, int, double, int> hyperpara)
        {
            Init(hyperpara.Item2);
            Shuffle(4 * KB.posTriples.Count);
            Console.WriteLine("Model " + mid.ToString() + " begins its iterations:");
            for (int t = 0; t < 1000; t++)
            {
                double energy = CalcEnergy(hyperpara.Item4, hyperpara.Item3);
                Console.WriteLine(t.ToString() + "-th round of model " + mid.ToString() + " has energy: " + energy.ToString());
                energylist.Add(energy);
                Shuffle(KB.posTriples.Count);
                //learning rate, C, batch size
                Round(hyperpara.Item1, hyperpara.Item3, hyperpara.Item4);
                if (t >= 500 && NoDecrease())
                {
                    break;
                }
            }
        }

        public void Output(string path0, string path1)
        {
            TextWriter ops = new StreamWriter(path0);
            foreach (string ename in KB.E.Keys)
            {
                ops.Write(ename + "\t");
                for (int i = 0; i < eeb[KB.E[ename]].Length; i++)
                {
                    ops.Write(eeb[KB.E[ename]][i].ToString() + ",");
                }
                ops.WriteLine();
            }
            ops.Close();
            ops = new StreamWriter(path1);
            foreach (string rname in KB.R.Keys)
            {
                ops.Write(rname + "\t");
                for (int i = 0; i < reb[KB.R[rname]].Length; i++)
                {
                    ops.Write(reb[KB.R[rname]][i].ToString() + ",");
                }
                ops.WriteLine();
            }
            ops.Close();
        }
    }
}
