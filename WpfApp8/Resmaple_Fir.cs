using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace WpfApp8
{
    public class Resmaple_Fir
    {
        // function[y, h] = uniformResample(x, p, q, N, bta )


        /// <summary>
        /// 
        /// </summary>
        /// <param name="xdata">待重采样信号</param>
        /// <param name="targetP">目标采样率</param>
        /// <param name="orginQ">待重采样信号的采样率</param>
        /// <param name="n">滤波器长度与n成正比，采用chebyshevIIR型低通滤波器的阶数，缺省值为10</param>
        /// <param name="bta">设置低通滤波器使使用Kaiser窗的参数，缺省值为5</param>
        public double[] uniformResample(double[] xdata, Int64 targetP, Int64 orginQ, int n = 10, int bta = 5)
        {
            Int64 greatestDivisor = greatest_common_divisor(targetP, orginQ);
            targetP /= greatestDivisor;
            orginQ /= greatestDivisor;

            Int64 maxPQ = Math.Max(targetP, orginQ);
            //fc = 1 / 2 / pqmax;
            double fc = 1.0 / 2 / maxPQ;
            //L = 2 * N * pqmax + 1;
            int L = (int)(2 * n * maxPQ + 1);
            //h = firls(L - 1, [0 2 * fc 2 * fc 1], [1 1 0 0]).* kaiser(L, bta)' ;
            //h = p * h / sum(h);

            double Lhalf = (L - 1) / 2;

            int Lx = xdata.Length;

            //滤波器的设计
            double[] fir = firls(L - 1, new double[] { 0, 2 * fc, 2 * fc, 1 }, new double[] { 1, 1, 0, 0 });

            double[] kwind = kaiser(L, bta);

            for (int i = 0; i < fir.Length; i++)
            {
                fir[i] *= kwind[i];
            }

            double firSum = fir.Sum();
            fir = fir.Select(a => a * targetP / firSum).ToArray();

            int firlength = fir.Length;
            //% Need to delay output so that downsampling by q hits center tap of filter.
            int nz = (int)Math.Floor(orginQ - Lhalf % orginQ);

            //调整滤波因子
             if (nz > 0)
            {
                Lhalf = Lhalf + nz;
                firlength += nz;
            }
            int delay = (int)Math.Floor(Math.Ceiling((Lhalf) / orginQ));
            int nz1 = 0;
            while (Math.Ceiling((double)(((Lx - 1) * targetP + firlength + nz1) / orginQ)) - delay < Math.Ceiling((double)(Lx * targetP / orginQ)))
            {
                nz1 += nz1;
            }
            double[] tempFir = new double[fir.Length];
            Array.Copy(fir, tempFir, fir.Length);
            fir = new double[fir.Length + nz + nz1];
            Array.Copy(tempFir, 0, fir, nz, tempFir.Length);

            //滤波后的数据
            double[] ydata = upfirdn(xdata, fir, (int)targetP, (int)orginQ);


            //output length
            int Ly = (int)Math.Ceiling((double)(Lx * targetP / orginQ));
            double[] yResult = new double[Ly];

            Array.Copy(ydata, delay, yResult, 0, Ly);
            return yResult;
            //if (nz > 0 && nz1 > 0)
            //{
            //    fir = new double[fir.Length + nz + nz1];
            //    Array.Copy(tempFir, 0, fir, nz, tempFir.Length);
            //}
            //else if (nz > 0 && nz1 == 0)
            //{
            //    fir = new double[fir.Length + nz+nz1];
            //    Array.Copy(tempFir, 0, fir, nz, tempFir.Length);
            //}
            //else if (nz <= 0 && nz1 > 0)
            //{
            //    fir = new double[fir.Length + nz + nz1];
            //    Array.Copy(tempFir, 0, fir, nz, tempFir.Length);
            //}


        }

        /// <summary>
        /// 求两个数最大公约数
        /// </summary>
        /// <param name="num1"></param>
        /// <param name="num2"></param>
        /// <returns></returns>
        private Int64 greatest_common_divisor(Int64 num1, Int64 num2)
        {
            Int64 a = num1;
            Int64 b = num2;
            Int64 tem;
            while (b != 0)
            {
                tem = a % b;
                a = b;
                b = tem;
            }
            return a;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n">Order of the filter</param>
        /// <param name="f">Pairs of frequency points</param>
        /// <param name="a">Amplitude values of the function at each frequency point, specified as a vector of the same length as f</param>
        public double[] firls(int n, double[] f, double[] m)
        {
            #region 判断条件是否满足

            //f 的值是否在0到1之间
            // 判断 f 个数是否为偶数
            // 判断 f ，m 长度是否为相等
            #endregion

            #region 确定参数 W 和ftype 并初始化为全1

            int[,] W = new int[2, 1] { { 1 }, { 1 } };

            int ftype = 0;
            int differ = 0;

            //% filter length
            n = n + 1;
            //make these guys columns
            double[] fcol = f.Select((a) => a / 2).ToArray();
            double[] mcol = m.ToArray();
            double[] wcol = { 1, 1 };

            //计算 f 的差分
            double[] df = new double[fcol.Length - 1];
            for (int i = 1; i < fcol.Length; i++)
            {
                df[i - 1] = fcol[i] - fcol[i - 1];
            }

            //判断差分是否小于0
            if (df.Min() < 0)
            {
                return null;
            }
            #endregion

            bool fullband = true;
            bool constant_weights = true;

            #region 算法
            int L = (n - 1) / 2;

            bool Nodd = n % 2 == 1 ? true : false;

            double b0 = 0;
            double[] b = new double[L];
            double[] k = new double[b.Length];
            for (int i = 0; i < k.Length; i++)
            {
                k[i] = i + 1;
            }

            var kk = k.Select((a) => a * a).ToArray();
            for (int i = 0; i < fcol.Length; i += 2)
            {
                double tempm = (mcol[i + 1] - mcol[i]) / (fcol[i + 1] - fcol[i]);
                double b1 = mcol[i] - tempm * fcol[i];
                if (Nodd)
                {
                    b0 = b0 + (b1 * (fcol[i + 1] - fcol[i]) + tempm / 2 * (fcol[i + 1] * fcol[i + 1] - fcol[i] * fcol[i])) *
                        Math.Abs(wcol[(i + 1) / 2] * wcol[(i + 1) / 2]);
                }

                //第一步处理系数；
                double pipi = tempm / (4 * Math.PI * Math.PI);
                var ktemp = k.Select((a) => pipi * (Math.Cos(2 * Math.PI * a * fcol[i + 1]) - Math.Cos(2 * Math.PI * a * fcol[i]))).ToArray();

                for (int bi = 0; bi < b.Length; bi++)
                {
                    b[bi] = b[bi] + (ktemp[bi] / kk[bi]) * Math.Abs(wcol[(i + 1) / 2] * wcol[(i + 1) / 2]);
                }


                //第二步处理
                var ktemp2 = k.Select((a) => (fcol[i + 1] * (tempm * fcol[i + 1] + b1) * (sinc(2 * a * fcol[i + 1])) -
                fcol[i] * (tempm * fcol[i] + b1) * (sinc(2 * a * fcol[i]))
                )).ToArray();

                for (int bi = 0; bi < b.Length; bi++)
                {
                    b[bi] = b[bi] + (ktemp2[bi]) * Math.Abs(wcol[(i + 1) / 2] * wcol[(i + 1) / 2]);
                }
            }


            #endregion

            #region 因子处理
            if (Nodd)
            {

                double[] aArray = new double[b.Length + 1];
                Array.Copy(b, 0, aArray, 1, b.Length);
                aArray[0] = b0;
                aArray = aArray.Select((a) => wcol[1] * wcol[1] * 4 * a / 2).ToArray();
                //aArray[0] = aArray[0] / 2;

                double[] h = new double[aArray.Length * 2 - 1];
                //复制 aArray；
                h[aArray.Length-1] = aArray[0];
                Array.Copy(aArray, 1, h, aArray.Length, aArray.Length - 1);
                //复制 aArray并反转；
                Array.Reverse(aArray);
                Array.Copy(aArray, 0, h, 0, aArray.Length - 1);

                return h;
            }


            #endregion
            return null;
        }


        private double sinc(double temp)
        {
            if (temp == 0)
            {
                return 1;
            }
            else
            {
                return Math.Sin(temp * Math.PI) / (temp * Math.PI);
            }
        }
        /// <summary>
        /// returns an L-point Kaiser window with shape factor beta
        /// </summary>
        /// <param name="L">Window length</param>
        /// <param name="bta">Shape factor</param>
        /// <returns></returns>
        public double[] kaiser(double L, double bta)
        {

            //nw = round(nn);
            //bes = abs(besseli(0, bta));
            //odd = rem(nw, 2);
            //xind = (nw - 1) ^ 2;
            //n = fix((nw + 1) / 2);
            int nw = (int)L;
            double bes = Math.Abs(Bessel.BesselI(0, bta));
            int odd = nw % 2;
            int xind = (nw - 1) * (nw - 1);
            int n = (nw + 1) / 2;

            double[] xi = new double[n];
            for (int i = 0; i < xi.Length; i++)
            {
                xi[i] = i + 0.5 * (1 - odd);
                xi[i] = 4 * xi[i] * xi[i];
                //besseli(0, bta * sqrt(1 - xi / xind)) / bes
                xi[i] = Bessel.BesselI(0, bta * Math.Sqrt(1 - xi[i] / xind)) / bes;
            }

            double[] resultW = new double[2 * n - odd];
            Array.Copy(xi, 0, resultW, n - 1, n);
            Array.Reverse(xi);
            Array.Copy(xi, 0, resultW, 0, n - odd);

            return resultW;

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="x">输入信号</param>
        /// <param name="h">滤波器</param>
        /// <param name="p">上采样</param>
        /// <param name="q">下采样</param>
        /// <returns></returns>
        public static double[] upfirdn(double[] x, double[] h, int p, int q)
        {
            int lx = x.Length;
            int lh = h.Length;

            int igv = p * q;
            int iv = q;
            int yStep = p;
            int ig = p;

            int ylength = h.Length + (x.Length - 1) * p / q;

            if (q > p)
            {
                ylength = (x.Length + h.Length) / (q /p);
            }
            //计算并定义结果的数据 暂定义为输出结果
            double[] y = new double[ylength];
   
            for (int i = 0; i < p; i++)
            {
                int yIndex = i;
                //权重的计算位置
                int pgIndex = (i * q) % p;
                int Lg = h.Length - pgIndex;
                Lg = Lg % p > 0 ? Lg / p + 1 : Lg / p;
                int rpq_offset = (i * q) / p;
                int pvEndIndex = rpq_offset;
                /*
               * PSEUDO-CODE for CONVOLUTION with GENERAL INCREMENTS:
               *
                *   w[n] = v[n] * g[n]
                *
                * Given:
                *   pointers:   pg, pv, and pw
                *   or arrays:  g[ ], v[ ], and w[ ]
                *   increments: ig, iv, and iw
                *   end points: h+Lh, x+Lx
                */
                /*
             * Region #1 (running onto the data):
             */
                int pgEndIndex = pgIndex + p * rpq_offset;
                while (pvEndIndex < x.Length && pgEndIndex < h.Length)
                {
                    double acc = 0.0;
                    int pvTempIndex = pvEndIndex;
                    int pgTempIndex = pgIndex;
                    while (pgTempIndex <= pgEndIndex)
                    {
                        acc += h[pgTempIndex] * x[pvTempIndex];
                        pvTempIndex--;
                        pgTempIndex += ig;

                    }
                    y[yIndex] += acc;
                    yIndex += yStep;
                    pvEndIndex += iv;
                    pgEndIndex += igv;
                }
                /*
                * Do we need to drain rest of signal?
                */
                if (pvEndIndex < x.Length)
                {
                    while (pgEndIndex >h.Length - 1)
                    {
                        pgEndIndex -= igv;
                    }
                    while (pvEndIndex < x.Length)
                    {
                        double acc = 0.0;
                        int pvTempIndex = pvEndIndex;
                        int pgTempIndex = pgIndex;
                        while (pgTempIndex <= pgEndIndex)
                        {
                            acc += h[pgTempIndex] * x[pvTempIndex];
                            pvTempIndex--;
                            pgTempIndex += ig;
                        }
                        y[yIndex] += acc;
                        yIndex += yStep;
                        pvEndIndex += iv;
                    }
                }
                else if (pgEndIndex < h.Length)
                {
                    double acc = 0.0;
                    int pvTempIndex = 0;
                    int pgTempIndex = x.Length - 1;
                    while (pvTempIndex < x.Length)
                    {
                        acc += h[pgTempIndex] * x[pvTempIndex];
                        pvTempIndex++;
                        pgTempIndex -= ig;
                    }
                    y[yIndex] += acc;
                    yIndex += yStep;
                    pgEndIndex += igv;
                    pvEndIndex += iv;

                }

                while (pgEndIndex > h.Length - 1)
                {
                    pgEndIndex -= igv;
                }
                int pvIndex = x.Length - Lg + 1;
                while (pvIndex < x.Length)
                {
                    double acc = 0.0;
                    int pvTempIndex = pvIndex;
                    int pgTempIndex = pgEndIndex;
                    while (pvTempIndex < x.Length)
                    {
                        acc += h[pgTempIndex] * x[pvTempIndex];
                        pvTempIndex++;
                        pgTempIndex -= ig;
                    }
                    y[yIndex] += acc;
                    yIndex += yStep;
                    pvIndex += iv;
                }

            }
            return y;
        }
    }
}