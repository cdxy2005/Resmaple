using Autofac;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.IO;
using MathNet.Numerics.IntegralTransforms;
using MathNet.Numerics;

namespace WpfApp8
{
    /// <summary>
    /// MainWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MainWindow : System.Windows.Window
    {
        public MainWindow()
        {
            InitializeComponent();
            //aaa();
           // double xxx = Bessel.BesselI(0, 5);
        }

        private ConcurrentBag<tempaa> tempList = new ConcurrentBag<tempaa>();
        private void aaa()
        {
            tempList = new ConcurrentBag<tempaa>();
            int deltaAngleCount = 300;
            int deltaRCount = 60000;
            double deltaAngle = 0.1 / 180 * Math.PI;
            double azimuth= 30 / 180 * Math.PI;
            double deltaR = 150;
            int add = 1;
            double locationX = 1000;
            double locationY = 1000;
            double x;
            double y;

            Stopwatch st = new Stopwatch();
            st.Start();
            Parallel.For(0, deltaAngleCount, (i) => {
                for (int j = 1; j <= deltaRCount; j++)
                {
                    y = Math.Sin((azimuth + i * deltaAngle)) * (j * deltaR) * add + locationX;
                    x = Math.Cos((azimuth + i * deltaAngle)) * (j * deltaR) * add + locationY;
                    //double y = Math.Sin((azimuth + i * deltaAngle) / 180 * Math.PI) * (j * deltaR) * add + locationX;
                    //double x = Math.Cos((azimuth + i * deltaAngle) / 180 * Math.PI) * (j * deltaR) * add + locationY;
                    //if (Math.Abs(x) > earthR || Math.Abs(x) > earthR)
                    //{
                    //    break;
                    //}
                    //AddRCS(x, y);
                    tempList.Add(new tempaa() { Myx = x, MyY = y });
                }
            });

            st.Stop();
          this.sss.Text=  st.Elapsed.TotalMilliseconds.ToString();
        }

        private List<tempaa> tempList2 = new List<tempaa>();
        private void bbb()
        {
            tempList2.Clear();
            int deltaAngleCount = 300;
            int deltaRCount = 60000;
            double deltaAngle = 0.1 / 180 * Math.PI;
            double azimuth = 30 / 180 * Math.PI;
            double deltaR = 150;
            int add = 1;
            double locationX = 1000;
            double locationY = 1000;
            double x;
            double y;

            Stopwatch st = new Stopwatch();
            st.Start();
            for (int i = 0; i < deltaAngleCount; i++)
            {
                for (int j = 1; j <= deltaRCount; j++)
                {
                    y = Math.Sin((azimuth + i * deltaAngle)) * (j * deltaR) * add + locationX;
                    x = Math.Cos((azimuth + i * deltaAngle)) * (j * deltaR) * add + locationY;
                    //double y = Math.Sin((azimuth + i * deltaAngle) / 180 * Math.PI) * (j * deltaR) * add + locationX;
                    //double x = Math.Cos((azimuth + i * deltaAngle) / 180 * Math.PI) * (j * deltaR) * add + locationY;
                    //if (Math.Abs(x) > earthR || Math.Abs(x) > earthR)
                    //{
                    //    break;
                    //}
                    //AddRCS(x, y);
                    tempList2.Add(new tempaa() { Myx = x, MyY = y });
                }
            }
            st.Stop();
            this.sss.Text = st.Elapsed.TotalMilliseconds.ToString();
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            aaa();
        }

        private void Button_Click_1(object sender, RoutedEventArgs e)
        {
            bbb();
        }

        private void ccc()
        {
            // Create your builder.
            var builder = new ContainerBuilder();

            // Usually you're only interested in exposing the type
            // via its interface:
          //  builder.RegisterType<SomeType>().As<IService>();

            // However, if you want BOTH services (not as common)
            // you can say so:
           // builder.RegisterType<SomeType>().AsSelf().As<IService>();
        }

        private void btnResmaple_Click(object sender, RoutedEventArgs e)
        {
            //double[] singalData= LFMSingnal();

            List<double> idatList = new List<double>();
            List<double> qdataList = new List<double>();

            readIQFromCSV(idatList, qdataList);

            double[] singalIData = idatList.ToArray();
            double[] singalQData = qdataList.ToArray();

            Resmaple_Fir resmaple = new Resmaple_Fir();
            double[] idata= resmaple.uniformResample(singalIData, (Int64)(4 *1e9),(Int64)(200*1e6));
            double[] qdata = resmaple.uniformResample(singalQData, (Int64)(4 * 1e9), (Int64)(200 * 1e6));

            downFreq(idata,qdata);

            string csvFile =  "tempIQdata.csv";


            writeToCSV(csvFile,idata,qdata);

        

            freqDomi("UPresompe.csv", idata, qdata, 200 * 1e6);
             
        }


        /// <summary>
        /// 下变频数据
        /// </summary>
        /// <param name="singalIData"></param>
        /// <param name="singalQData"></param>
        private void downFreq(double[] singalIData, double[] singalQData)
        {
            Resmaple_Fir resmaple = new Resmaple_Fir();
            double[] idata = resmaple.uniformResample(singalIData, (Int64)(160 * 1e6), (Int64)(4 * 1e9));

            //double[] singalDataq = LFMSingnal(false);
            //  double[] singalDataqDouble = singalDataq.Select(a => (double)a).ToArray();
            double[] qdata = resmaple.uniformResample(singalQData, (Int64)(160 * 1e6), (Int64)(4 * 1e9));

            freqDomi("Downresompe.csv", idata, qdata, 160 * 1e6);
        }

        /// <summary>
        /// 获取线性调频信号中的 IQ数据
        /// </summary>
        /// <param name="isIdata"></param>
        /// <returns></returns>
        private double[] LFMSingnal(bool isIdata=true)
        {
            double T = 10 * 1e-6; // 0.00001;
            double B = 1 * 1e6;  // 1000000;
            double k = B / T;    // 调频斜率
            double t = 5 * 1e-6;// 0.000005;% 空信号时长5us
            double fs = 4 * 1e9;// 4000000000;

            // 生成200Msa线性调频信号 5us空 + 10us线性调频
            double fs1 = 200 * 1e6 ;
            int n1 =(int)Math.Round(T * fs1); //采样点个数
            //double[] t1 = new double[n1];
            double tspan = T / n1;
            MathNet.Numerics.Complex32[] y1 = new MathNet.Numerics.Complex32[n1];
            for (int i = 0; i < n1; i++)
            {
                 
                MathNet.Numerics.Complex32 tempComple = new MathNet.Numerics.Complex32(0, (float)(Math.PI * k * (i * tspan)* (i * tspan)));

                // double aa= tempComple* k;
                y1[i] = MathNet.Numerics.Complex32.Exp(tempComple);
            }



            double[] idata = y1.Select(a => (double)a.Real).ToArray();
            double[] qdata = y1.Select(a => (double)a.Imaginary).ToArray();

           // freqDomi("LFM_freq.csv",  idata, qdata, 200 * 1e6);
            double[] result = new double[1000 + n1];
            if (isIdata)
            {
                Array.Copy(idata, 0, result, 1000, n1);
                //result = idata;
            }
            else
            {
                 Array.Copy(qdata, 0, result, 1000, n1);
                //result = qdata;
            }


            return result;
            //y1 = exp(1j * pi * k * t1.^ 2);% LFM信号
            // MathNet.Numerics.Complex32.Exp()

        }

        private void readIQFromCSV(List<double> idataList, List<double> qdataList)
        {
            string fileName = @"多音信号.csv";
            string csvFile = AppDomain.CurrentDomain.BaseDirectory + fileName;
            string strLin = "";
           // List<double> idataList = new List<double>();
           // List<double> qdataList = new List<double>();
            double idata = double.NaN;
            double qdata = double.NaN;
            using (FileStream fs = new FileStream(csvFile, System.IO.FileMode.Open, FileAccess.ReadWrite))
            {
                using (StreamReader sr = new StreamReader(fs))
                {
                    while ((strLin = sr.ReadLine()) != null)
                    {
                        string[] iqStr = strLin.Split(new char[] { ',' });
                        if (double.TryParse(iqStr[0], out idata) && double.TryParse(iqStr[1], out qdata))
                        {
                            idataList.Add(idata);
                            qdataList.Add(qdata);
                        }
                    }

                }
            }


        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="signal">采样点数</param>
        /// <param name="freq">采样频率</param>
        private void freqDomi(string fileName, double[] idata, double[] qdata,double freq)
        {

            Fourier.Forward(idata,qdata);

            double freqStep = freq / idata.Length;
            double[] freqs = new double[idata.Length];
            for (int i = 0; i < freqs.Length; i++)
            {
                freqs[i]= i *freqStep;
            }

            Complex32[] complexData = new Complex32[idata.Length];
            for (int i = 0; i < idata.Length; i++)
            {
                complexData[i] = new Complex32((float)idata[i], (float)qdata[i]);
            }

            double[] amt = complexData.Select(a => (double)a.Magnitude).ToArray();
            double[] phases = complexData.Select(a => (double)a.Phase).ToArray();
            writeToCSV(fileName,  amt, phases);
        }

        /// <summary>
        /// 写入到CSV文件
        /// </summary>
        private void writeToCSV(string fileName,double[]idata,double[]qdata)
        {
            string csvFile = AppDomain.CurrentDomain.BaseDirectory + fileName;
            using (FileStream fs = new FileStream(csvFile, System.IO.FileMode.Create, FileAccess.ReadWrite))
            {
                using (StreamWriter sw = new StreamWriter(fs))
                {
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < idata.Length && i < qdata.Length; i++)
                    {
                        sb.Append(idata[i]);
                        sb.Append(",");
                        sb.Append(qdata[i]);
                        sb.Append("\n");

                    }
                    sw.Write(sb.ToString());
                }
            }
        }

        private void btnResmapleLFm_Click(object sender, RoutedEventArgs e)
        {
           
            double[] singalIData = LFMSingnal(); ;
            double[] singalQData = LFMSingnal(false);

            Resmaple_Fir resmaple = new Resmaple_Fir();
            double[] idata = resmaple.uniformResample(singalIData, (Int64)(4 * 1e9), (Int64)(200 * 1e6));
            double[] qdata = resmaple.uniformResample(singalQData, (Int64)(4 * 1e9), (Int64)(200 * 1e6));

            //downFreq(idata,qdata);

            string csvFile = "tempIQdata.csv";

            writeToCSV(csvFile, idata, qdata);

            freqDomi("UPresompe.csv", idata, qdata, 200 * 1e6);
        }
    }

    public class tempaa
    {
        private double _myX;

        public double Myx
        {
            get { return _myX; }
            set { _myX = value; }
        }

        private double _myY;

        public double MyY
        {
            get { return _myY; }
            set { _myY = value; }
        }


    }
}
