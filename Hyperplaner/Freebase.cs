using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Hyperplaner
{
    class Freebase
    {
        public Dictionary<string, int> E = new Dictionary<string, int>();
        public Dictionary<string, int> R = new Dictionary<string, int>();
        public HashSet<Tuple<int, int, int>> I = new HashSet<Tuple<int, int, int>>();
        public List<Tuple<int, int, int>> posTriples = new List<Tuple<int, int, int>>();
        public List<Tuple<int, int, int>> negTriples = new List<Tuple<int, int, int>>();

        public Freebase(string inputpath)
        {
            TextReader ips = new StreamReader(inputpath);
            while (ips.Peek() != -1)
            {
                string[] eles = ips.ReadLine().Split(new char[] { '\t' }, StringSplitOptions.RemoveEmptyEntries);
                if (eles.Length == 3)
                {
                    if (E.ContainsKey(eles[0]) == false) E.Add(eles[0], E.Count);
                    if (E.ContainsKey(eles[2]) == false) E.Add(eles[2], E.Count);
                    if (R.ContainsKey(eles[1]) == false) R.Add(eles[1], R.Count);
                    Tuple<int, int, int> curt = new Tuple<int, int, int>(E[eles[0]], R[eles[1]], E[eles[2]]);
                    I.Add(curt);
                    posTriples.Add(curt);
                }
            }
            ips.Close();
        }

        public Tuple<int, int, int> ConstructNegEg(Tuple<int, int, int> post)
        {
            Random rd = new Random();
            int ch = rd.Next(E.Count + E.Count + R.Count - 5);
            if (ch < R.Count - 1)
            {
                while (true)
                {
                    Tuple<int, int, int> negt = new Tuple<int, int, int>(post.Item1, rd.Next(R.Count), post.Item3);
                    if (I.Contains(negt) == false) return negt;
                }
            }
            else
            {
                while (true)
                {
                    int prime = rd.Next(E.Count);
                    if (prime != post.Item1 && prime != post.Item3)
                    {
                        Tuple<int, int, int> negt = ch % 2 == 1 ? new Tuple<int, int, int>(post.Item1, post.Item2, prime) : new Tuple<int, int, int>(prime, post.Item2, post.Item3);
                        if (I.Contains(negt) == false) return negt;
                    }
                }
            }
        }

        public void BuildNegative()
        {
            foreach (var en in posTriples)
            {
                negTriples.Add(ConstructNegEg(en));
            }
        }

        public void Split(string inputpath)
        {
            TextReader ips = new StreamReader(inputpath);
            TextWriter ops0 = new StreamWriter(inputpath + "_train");
            TextWriter ops1 = new StreamWriter(inputpath + "_validation");
            TextWriter ops2 = new StreamWriter(inputpath + "_test");
            Random rd = new Random();
            string curR = "";
            List<Tuple<string,string>> buf = new List<Tuple<string,string>>();
            int ctR = 0;

            while (ips.Peek() != -1)
            {
                string[] eles = ips.ReadLine().Split(new char[]{'\t'},StringSplitOptions.RemoveEmptyEntries);
                if ( eles.Length == 1 )
                {
                    if (buf.Count >= 100)
                    {
                        ctR++;
                        foreach (var en in buf)
                        {
                            int ch = rd.Next(10);
                            if (ch <= 5)
                            {
                                ops0.WriteLine(en.Item1 + "\t" + curR + "\t" + en.Item2);
                            }
                            else if (ch == 6)
                            {
                                ops1.WriteLine(en.Item1 + "\t" + curR + "\t" + en.Item2);
                            }
                            else
                            {
                                ops2.WriteLine(en.Item1 + "\t" + curR + "\t" + en.Item2);
                            }
                        }
                    }
                    curR = eles[0];
                    buf.Clear();
                }
                else if ( eles.Length == 2)
                {
                    buf.Add(new Tuple<string, string>(eles[0], eles[1]));
                }
            }
            if (buf.Count >= 100)
                    {
                        ctR++;
                        foreach (var en in buf)
                        {
                            int ch = rd.Next(10);
                            if (ch <= 5)
                            {
                                ops0.WriteLine(en.Item1 + "\t" + curR + "\t" + en.Item2);
                            }
                            else if (ch == 6)
                            {
                                ops1.WriteLine(en.Item1 + "\t" + curR + "\t" + en.Item2);
                            }
                            else
                            {
                                ops2.WriteLine(en.Item1 + "\t" + curR + "\t" + en.Item2);
                            }
                        }
                    }
            ips.Close();
            ops0.Close();
            ops1.Close();
            ops2.Close();
            Console.WriteLine("There are " + ctR.ToString() + " relations in this filtered dataset.");
        }
    }
}
