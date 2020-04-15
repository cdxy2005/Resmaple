/* =============================================
 * Copyright @ 2017 北京中创锐科技术有限公司 
 * CLR 版本：4.0.30319.42000 
 * 当前版本：1.0.0.1 
 * 名    称：Bessel  
 * 功    能：Bessel  
 * 作    者：CHEN XF Administrator
 * 添加时间：2020/3/31 11:01:28
 * =============================================*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using Complex = System.Numerics.Complex;

namespace WpfApp8
{
    public static partial class Bessel
    {

        /// <summary>
        /// Returns the modified Bessel function of the first kind.
        /// <para>BesselI(n, z) is a solution to the modified Bessel differential equation.</para>
        /// </summary>
        /// <param name="n">The order of the modified Bessel function.</param>
        /// <param name="z">The value to compute the modified Bessel function of.</param>
        /// <returns>The modified Bessel function of the first kind.</returns>
        public static double BesselI(double n, double z)
        {
            return BesselI(n, new Complex(z, 0)).Real;
        }

        /// <summary>
        /// Returns the modified Bessel function of the first kind.
        /// <para>BesselI(n, z) is a solution to the modified Bessel differential equation.</para>
        /// </summary>
        /// <param name="n">The order of the modified Bessel function.</param>
        /// <param name="z">The value to compute the modified Bessel function of.</param>
        /// <returns>The modified Bessel function of the first kind.</returns>
        public static Complex BesselI(double n, Complex z)
        {
            return Amos.Cbesi(n, z);
        }
        private static class Amos
        {
            #region BesselI

            // Returns I(v, z)
            public static Complex Cbesi(double v, Complex z)
            {
                if (double.IsNaN(v) || double.IsNaN(z.Real) || double.IsNaN(z.Imaginary))
                {
                    return new Complex(double.NaN, double.NaN);
                }

                int sign = 1;
                if (v < 0)
                {
                    v = -v;
                    sign = -1;
                }

                int n = 1;
                int kode = 1;
                int nz = 0;
                int ierr = 0;

                double[] cyir = new double[n];
                double[] cyii = new double[n];
                for (int i = 0; i < n; i++)
                {
                    cyir[i] = double.NaN;
                    cyii[i] = double.NaN;
                }

                AmosHelper.zbesi(z.Real, z.Imaginary, v, kode, n, cyir, cyii, ref nz, ref ierr);
                Complex cyi = new Complex(cyir[0], cyii[0]);

                if (ierr == 2)
                {
                    //overflow
                    if (z.Imaginary == 0 && (z.Real >= 0 || v == Math.Floor(v)))
                    {
                        if (z.Real < 0 && v / 2 != Math.Floor(v / 2))
                            cyi = new Complex(double.NegativeInfinity, 0);
                        else
                            cyi = new Complex(double.PositiveInfinity, 0);
                    }
                    else
                    {
                        cyi = ScaledCbesi(v * sign, z);
                        cyi = new Complex(cyi.Real * double.PositiveInfinity, cyi.Imaginary * double.PositiveInfinity);
                    }
                }

                if (sign == -1)
                {
                    if (!ReflectI(cyi, v))
                    {
                        double[] cykr = new double[n];
                        double[] cyki = new double[n];
                        AmosHelper.zbesk(z.Real, z.Imaginary, v, kode, n, cykr, cyki, ref nz, ref ierr);
                        Complex cyk = new Complex(cykr[0], cyki[0]);

                        cyi = RotateI(cyi, cyk, v);
                    }
                }

                return cyi;
            }

            // Return Exp(-Abs(x)) * I(v, z) where x = z.Real
            public static Complex ScaledCbesi(double v, Complex z)
            {
                if (double.IsNaN(v) || double.IsNaN(z.Real) || double.IsNaN(z.Imaginary))
                {
                    return new Complex(double.NaN, double.NaN);
                }

                int sign = 1;
                if (v < 0)
                {
                    v = -v;
                    sign = -1;
                }

                int n = 1;
                int kode = 2;
                int nz = 0;
                int ierr = 0;

                double[] cyir = new double[n];
                double[] cyii = new double[n];
                for (int i = 0; i < n; i++)
                {
                    cyir[i] = double.NaN;
                    cyii[i] = double.NaN;
                }

                AmosHelper.zbesi(z.Real, z.Imaginary, v, kode, n, cyir, cyii, ref nz, ref ierr);
                Complex cyi = new Complex(cyir[0], cyii[0]);

                if (sign == -1)
                {
                    if (!ReflectI(cyi, v))
                    {
                        double[] cykr = new double[n];
                        double[] cyki = new double[n];
                        AmosHelper.zbesk(z.Real, z.Imaginary, v, kode, n, cykr, cyki, ref nz, ref ierr);
                        Complex cyk = new Complex(cykr[0], cyki[0]);

                        //adjust scaling to match zbesi
                        cyk = Rotate(cyk, -z.Imaginary / Math.PI);
                        if (z.Real > 0)
                        {
                            cyk = new Complex(cyk.Real * Math.Exp(-2 * z.Real), cyk.Imaginary * Math.Exp(-2 * z.Real));
                        }
                        //v -> -v
                        cyi = RotateI(cyi, cyk, v);
                    }
                }

                return cyi;
            }

            #region utilities

            private static double SinPi(double x)
            {
                if (Math.Floor(x) == x && Math.Abs(x) < 1.0e14)
                {
                    //Return 0 when at exact zero, as long as the floating point number is
                    //small enough to distinguish integer points from other points.

                    return 0;
                }
                return Math.Sin(Math.PI * x);
            }

            private static double CosPi(double x)
            {
                if (Math.Floor(x + 0.5) == x + 0.5 && Math.Abs(x) < 1.0E14)
                {
                    //Return 0 when at exact zero, as long as the floating point number is
                    //small enough to distinguish integer points from other points.

                    return 0;
                }
                return Math.Cos(Math.PI * x);
            }

            private static Complex Rotate(Complex z, double v)
            {
                double c = CosPi(v);
                double s = SinPi(v);
                return new Complex(z.Real * c - z.Imaginary * s, z.Real * s + z.Imaginary * c);
            }

            private static Complex RotateJY(Complex j, Complex y, double v)
            {
                double c = CosPi(v);
                double s = SinPi(v);
                return new Complex(j.Real * c - y.Real * s, j.Imaginary * c - y.Imaginary * s);
            }

            private static bool ReflectJY(ref Complex jy, double v)
            {
                //NB: Y_v may be huge near negative integers -- so handle exact
                //     integers carefully

                if (v != Math.Floor(v))
                {
                    return false;
                }

                int i = (int)(v - 16384.0 * Math.Floor(v / 16384.0));
                if (i % 2 == 1)
                {
                    jy = new Complex(-jy.Real, -jy.Imaginary);
                }

                return true;
            }

            private static bool ReflectI(Complex ik, double v)
            {
                if (v != Math.Floor(v))
                {
                    return false;
                }

                return true; //I is symmetric for integer v
            }

            private static Complex RotateI(Complex i, Complex k, double v)
            {
                double s = Math.Sin(v * Math.PI) * (2.0 / Math.PI);
                return new Complex(i.Real + s * k.Real, i.Imaginary + s * k.Imaginary);
            }

            #endregion

            #endregion
        }
    }

    
    }
