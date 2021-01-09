using System;
using System.Collections.Generic;
using System.Text;

namespace ElectronDensity
{
    public class DensityPoint : IComparable
    {
        public EdVector XYZ = new EdVector(0, 0, 0);
        public EdVector CRS = new EdVector(0, 0, 0);

        public Dictionary<int, double> values = new Dictionary<int, double>();
        public double v;
        public double distance = 0;
        public bool valid = true;

        public DensityPoint(double i, double j, double k, double vv, string pointtype)
        {
            if (pointtype == "XYZ")
            {
                XYZ.x = i;
                XYZ.y = j;
                XYZ.z = k;
            }
            else
            {
                CRS.x = i;
                CRS.y = j;
                CRS.z = k;
            }
            v = vv;
        }

        public DensityPoint(DensityPoint copyPoint)
        {
            XYZ.x = copyPoint.XYZ.x;
            XYZ.y = copyPoint.XYZ.y;
            XYZ.z = copyPoint.XYZ.z;
            CRS.x = copyPoint.CRS.x;
            CRS.y = copyPoint.CRS.y;
            CRS.z = copyPoint.CRS.z;
            v = copyPoint.v;
            values = copyPoint.values;
            valid = copyPoint.valid;

        }

        public void setDerivativeasMain(int diff)
        {
            v = values[diff];
        }

        public void setXYZ(double xx, double yy, double zz)
        {
            XYZ.x = xx;
            XYZ.y = yy;
            XYZ.z = zz;
        }

        public void setCRS(double cc, double rr, double ss)
        {
            CRS.x = cc;
            CRS.y = rr;
            CRS.z = ss;
        }

        public EdVector getXYZ()
        {
            /*DenseVector<Double> p = Vector.Create<Double>(3);
            p[0] = x;
            p[1] = y;
            p[2] = z;*/
            return XYZ;
        }

        public void movePoint(double distance, EdVector direction)
        {
            //TODO for now just add the amount to every point
            XYZ.x += distance * direction.x;
            XYZ.y += distance * direction.y;
            XYZ.z += distance * direction.z;
        }

        public string getPointStringXYZ()
        {
            return "(" + Convert.ToString(Math.Round(X, 2)) + "," + Convert.ToString(Math.Round(Y, 2)) + "," + Convert.ToString(Math.Round(Z, 2)) + ")";
        }
        public string getPointStringCRS()
        {
            return "(" + Convert.ToString(Math.Round(C, 0)) + "," + Convert.ToString(Math.Round(R, 0)) + "," + Convert.ToString(Math.Round(S, 0)) + ")";
        }

        public string getPointStringDetails()
        {
            return "Density=" + Convert.ToString(V) + " XYZ=:" + getPointStringXYZ() + " CRS=:" + getPointStringCRS();
        }

        //get and sets
        public double X
        {
            get
            {
                return XYZ.x;
            }
        }
        public double X2
        {
            get
            {
                return Math.Round(XYZ.x, 2);
            }
        }
        public double Y
        {
            get
            {
                return XYZ.y;
            }
        }
        public double Y2
        {
            get
            {
                return Math.Round(XYZ.y, 2);
            }
        }
        public double Z
        {
            get
            {
                return XYZ.z;
            }
        }
        public double Z2
        {
            get
            {
                return Math.Round(XYZ.z, 2);
            }
        }
        public double C
        {
            get
            {
                return CRS.x;
            }
        }
        public double C0
        {
            get
            {
                return Math.Round(CRS.x, 0);
            }
        }
        public double R
        {
            get
            {
                return CRS.y;
            }
        }
        public double R0
        {
            get
            {
                return Math.Round(CRS.y, 0);
            }
        }
        public double S
        {
            get
            {
                return CRS.z;
            }
        }
        public double S0
        {
            get
            {
                return Math.Round(CRS.z, 0);
            }
        }
        public double V
        {
            get
            {
                return Math.Round(v, 3);
            }
        }

        public double D
        {
            get
            {
                return Math.Round(distance, 3);
            }
        }

        public int CompareTo(object obj)
        {
            if (obj == null) return 1;
            if (v < (obj as DensityPoint).v)
                return 1;
            if (v > (obj as DensityPoint).v)
                return -1;
            else
                return 0;
        }


    }
}
