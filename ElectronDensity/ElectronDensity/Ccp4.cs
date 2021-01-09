using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace ElectronDensity
{
    public class Ccp4
    {
        /*         
        * https://www.ccpem.ac.uk/mrc_format/mrc2014.php 
        */
        string _pdbCode;
        string _edFilePath;
        string _edWebPath;

        const int MagicInvalid = -1000;//this is a number if invalid, it would be better to use the minimum of the structure -1 perhaps.        
        public double[,,] DensityMatrix;

        /*
         * These are the "words" from the ccp4 file
         */
        // 1-3

        public int w01_NX = 0; //Column dim, the fastest changing axis
        public int w02_NY = 0; //Row dim, the fastest changing axis
        public int w03_NZ = 0; //Sector dim, the slowest changing axis
        public int w05_NXSTART = 0;
        public int w06_NYSTART = 0;
        public int w07_NZSTART = 0;
        public int w08_MX = 0;
        public int w09_MY = 0;
        public int w10_MZ = 0;
        public double w11_CELLA_X = 0;
        public double w12_CELLA_Y = 0;
        public double w13_CELLA_Z = 0;
        public double w14_CELLB_X = 0;
        public double w15_CELLB_Y = 0;
        public double w16_CELLB_Z = 0;
        public int w17_MAPC = 0;
        public int w18_MAPR = 0;
        public int w19_MAPS = 0;
        public double w20_DMIN = 0;
        public double w21_DMAX = 0;
        public double w22_DMEAN = 0;

        //Calculated from the "Words"

        private EdMatrix _orthoMat;
        private EdMatrix _deOrthoMat;
        private EdVector _origin = new EdVector();
        private int[] _map2xyz;
        private int[] _map2crs;
        private double[] _cellDims;
        private int[] _axisSampling;
        private int[] _crsStart;
        private int[] _dimOrder;

        public string InfoDim = "";
        public string InfoCellLengths = "";
        public string InfoCellAngles = "";

        public int Degree = 0;
        public int Differential = 0;

        public Ccp4(string pdbCode, string edPath, bool isDiff)
        {            
            if (pdbCode != "")            
                setup(pdbCode, edPath, isDiff);            
        }
        private void setup(string pdbCode, string edPath, bool isDiff)
        {
            _pdbCode = pdbCode;
            if (isDiff)
            {
                _edFilePath = edPath + pdbCode + "_diff.ccp4";
                _edWebPath = @"https://www.ebi.ac.uk/pdbe/coordinates/files/" + _pdbCode + "_diff.ccp4";
            }
            else
            {
                _edFilePath = edPath + pdbCode + ".ccp4";
                _edWebPath = @"https://www.ebi.ac.uk/pdbe/coordinates/files/" + _pdbCode + ".ccp4";
            }

            if (!File.Exists(_edFilePath))
            {
                // Create a new WebClient instance.
                using (System.Net.WebClient myWebClient = new System.Net.WebClient())
                {
                    // Download the Web resource and save it into the current filesystem folder.
                    myWebClient.DownloadFile(_edWebPath, _edFilePath);
                }
            }

            byte[] fileInBinary = ReadBinaryFile(_edFilePath);
            w01_NX = bytesToInt(fileInBinary, 0); // 1
            w02_NY = bytesToInt(fileInBinary, 4); // 2
            w03_NZ = bytesToInt(fileInBinary, 8); // 3
            int MODE = bytesToInt(fileInBinary, 12); // 4
            w05_NXSTART = bytesToInt(fileInBinary, 16); // 5
            w06_NYSTART = bytesToInt(fileInBinary, 20); // 6
            w07_NZSTART = bytesToInt(fileInBinary, 24); // 7
            w08_MX = bytesToInt(fileInBinary, 28); // 8
            w09_MY = bytesToInt(fileInBinary, 32); // 9 
            w10_MZ = bytesToInt(fileInBinary, 36); // 10
            w11_CELLA_X = bytesToSingle(fileInBinary, 40); // 11
            w12_CELLA_Y = bytesToSingle(fileInBinary, 44); // 12
            w13_CELLA_Z = bytesToSingle(fileInBinary, 48); // 13
            w14_CELLB_X = Convert.ToSingle(Math.Round(Convert.ToDouble(bytesToSingle(fileInBinary, 52)), 2)); // 14
            w15_CELLB_Y = Convert.ToSingle(Math.Round(Convert.ToDouble(bytesToSingle(fileInBinary, 56)), 2)); // 15
            w16_CELLB_Z = Convert.ToSingle(Math.Round(Convert.ToDouble(bytesToSingle(fileInBinary, 60)), 2)); // 16
            w17_MAPC = bytesToInt(fileInBinary, 64) - 1; // 17
            w18_MAPR = bytesToInt(fileInBinary, 68) - 1; // 18
            w19_MAPS = bytesToInt(fileInBinary, 72) - 1; // 19
            w20_DMIN = bytesToSingle(fileInBinary, 76); // 20
            w21_DMAX = bytesToSingle(fileInBinary, 80); // 21
            w22_DMEAN = bytesToSingle(fileInBinary, 84); // 22
            int ISPG = bytesToInt(fileInBinary, 88); // 23
            int NYSYMBT = bytesToInt(fileInBinary, 92); // 24
            //EXTRA
            int EXTTYP = bytesToInt(fileInBinary, 104); // 27
            int NVERSION = bytesToInt(fileInBinary, 108); // 28
            //EXTRA
            Single ORIGIN_X = bytesToSingle(fileInBinary, 196); // 50
            Single ORIGIN_Y = bytesToSingle(fileInBinary, 200); // 51
            Single ORIGIN_Z = bytesToSingle(fileInBinary, 204); // 52
            string MAP = bytesToString(fileInBinary, 208, 4); // 53
            Single RMS = bytesToSingle(fileInBinary, 216); // 55
            int NLABL = bytesToInt(fileInBinary, 220); // 56
            string LABEL1 = bytesToString(fileInBinary, 224, 80); // 57
            string LABEL2 = bytesToString(fileInBinary, 304, 80); // 58
            string LABEL3 = bytesToString(fileInBinary, 384, 80); // 59
            string LABEL4 = bytesToString(fileInBinary, 464, 80); // 60
            string LABEL5 = bytesToString(fileInBinary, 544, 80); // 61
            string LABEL6 = bytesToString(fileInBinary, 624, 80); // 62
            string LABEL7 = bytesToString(fileInBinary, 704, 80); // 63
            string LABEL8 = bytesToString(fileInBinary, 784, 80); // 64
            string LABEL9 = bytesToString(fileInBinary, 864, 80); // 65
            string LABEL10 = bytesToString(fileInBinary, 944, 80); // 66

            //And now get the actual matrx data as a list
            List<Single> theMatrixData = bytesToSingles(fileInBinary, NYSYMBT);

            int count = 0;
            DensityMatrix = new double[w03_NZ, w02_NY, w01_NX];
            for (int i = 0; i < w03_NZ; ++i)
            {
                for (int j = 0; j < w02_NY; ++j)
                {
                    for (int k = 0; k < w01_NX; ++k)
                    {
                        DensityMatrix[i, j, k] = theMatrixData[count];
                        count += 1;
                    }
                }
            }
            //And now make some calculations
            calculateOrthoMat();
            calculateOrigin();
            _map2xyz = new int[3];
            _map2xyz[w17_MAPC] = 0;
            _map2xyz[w18_MAPR] = 1;
            _map2xyz[w19_MAPS] = 2;
            
            _map2crs = new int[3];
            _map2crs[0] = w17_MAPC;
            _map2crs[1] = w18_MAPR;
            _map2crs[2] = w19_MAPS;
            
            _cellDims = new double[3];
            _cellDims[0] = w11_CELLA_X;
            _cellDims[1] = w12_CELLA_Y;
            _cellDims[2] = w13_CELLA_Z;

            _axisSampling = new int[3];
            _axisSampling[0] = w08_MX;
            _axisSampling[1] = w09_MY;
            _axisSampling[2] = w10_MZ;

            _crsStart = new int[3];
            _crsStart[0] = w05_NXSTART;
            _crsStart[1] = w06_NYSTART;
            _crsStart[2] = w07_NZSTART;

            _dimOrder = new int[3];
            _dimOrder[0] = w01_NX;
            _dimOrder[1] = w02_NY;
            _dimOrder[2] = w03_NZ;

            InfoDim = Convert.ToString(w01_NX) + ":" + Convert.ToString(w02_NY) + ":" + Convert.ToString(w03_NZ);
            InfoCellLengths = Convert.ToString(Math.Round(w11_CELLA_X / w08_MX, 2)) + ":" + Convert.ToString(Math.Round(w12_CELLA_Y / w09_MY, 2)) + ":" + Convert.ToString(Math.Round(w13_CELLA_Z / w10_MZ, 2));
            InfoCellAngles = Convert.ToString(Math.Round(w14_CELLB_X, 0)) + ":" + Convert.ToString(Math.Round(w15_CELLB_Y, 0)) + ":" + Convert.ToString(Math.Round(w16_CELLB_Z, 0));
        }

        public static byte[] ReadBinaryFile(string filePath)
        {
            byte[] buffer;
            FileStream fileStream = new FileStream(filePath, FileMode.Open, FileAccess.Read);
            try
            {
                int length = (int)fileStream.Length;  // get file length
                buffer = new byte[length];            // create buffer
                int count;                            // actual number of bytes read
                int sum = 0;                          // total number of bytes read

                // read until Read method returns 0 (end of the stream has been reached)
                while ((count = fileStream.Read(buffer, sum, length - sum)) > 0)
                    sum += count;  // sum is a buffer offset for next reading
            }
            finally
            {
                fileStream.Close();
            }

            //if (BitConverter.IsLittleEndian)
            //    Array.Reverse(buffer);

            return buffer;
        }

        private string bytesToString(byte[] bytes, int start, int length)
        {
            byte[] result = new byte[length];
            Array.Copy(bytes, start, result, 0, length);
            using (var stream = new MemoryStream(result))
            {
                using (var streamReader = new StreamReader(stream))
                {
                    return streamReader.ReadToEnd();
                }
            }
        }
        private int bytesToInt(byte[] bytes, int start)
        {
            int i = BitConverter.ToInt32(bytes, start);
            return i;
        }

        private Single bytesToSingle(byte[] bytes, int start)
        {
            Single value = BitConverter.ToSingle(bytes, start);
            return value;
        }

        private List<Single> bytesToSingles(byte[] bytes, int start)
        {
            int len = w01_NX * w02_NY * w03_NZ;
            start = bytes.Length - (4 * len);

            List<Single> matvals = new List<Single>();
            for (int i = start; i < bytes.Length; i += 4)
            {
                Single value = BitConverter.ToSingle(bytes, i);
                matvals.Add(value);
            }
            return matvals;
        }

        public DensityPoint getDensityPointFromCRS(double c, double r, double s, bool calcValue)
        {
            DensityPoint dp = new DensityPoint(c, r, s, 0, "CRS");
            DensityPoint dp2 = getXYZFromCRS(c, r, s);
            dp.setXYZ(dp2.X, dp2.Y, dp2.Z);
            if (calcValue)
            {
                DensityPoint v = getCRSValue(dp.C, dp.R, dp.S);
                dp.v = v.v;
                dp.values = v.values;
            }
            return dp;
        }

        

        public double getCRSGridPoint(int c, int r, int s)
        {
            int getC = c < 0 ? 0 : c;
            int getR = r < 0 ? 0 : r;
            int getS = s < 0 ? 0 : s;
            getC = getC >= w01_NX ? w01_NX - 1 : getC;
            getR = getR >= w02_NY ? w02_NY - 1 : getR;
            getS = getS >= w03_NZ ? w03_NZ - 1 : getS;
            double val2 = DensityMatrix[getS, getR, getC];
            return val2;
        }
        public bool isValidPoint(DensityPoint p)
        {
            return isValidCRS((int)p.C, (int)p.R, (int)p.S);
        }
        public bool isValidCRS(int c, int r, int s)
        {
            bool valid = true;
            if (c < 0)
                valid = false;
            if (r < 0)
                valid = false;
            if (s < 0)
                valid = false;
            if (c >= w01_NX)
                valid = false;
            if (r >= w02_NY)
                valid = false;
            if (s >= w03_NZ)
                valid = false;
            return valid;
        }

        /*
         * The function below uses an interpolation method written by me :-) 
         */
        public DensityPoint getCRSValue(Double c, Double r, Double s)
        {
            DensityPoint centre = new DensityPoint(c, r, s, 0, "CRS");
            DensityPoint vals = getInterpolatedDensity(centre);
            return vals;
        }

        private DensityPoint getInterpolatedDensity(DensityPoint centre)
        {
            //1. Build the points from which we will interpolate
            //Going to create the interpolated points from the perspective of 3 axes and average

            List<DensityPoint> pointsA = new List<DensityPoint>();
            List<DensityPoint> pointsB = new List<DensityPoint>();
            List<DensityPoint> pointsC = new List<DensityPoint>();

            int offset = 1;
            int numPoints = 2 * offset;
            numPoints = 2;
            int clStart = Convert.ToInt16(Math.Floor(centre.C)) - offset + 1;
            int rlStart = Convert.ToInt16(Math.Floor(centre.R)) - offset + 1;
            int slStart = Convert.ToInt16(Math.Floor(centre.S)) - offset + 1;
            if (numPoints == 1)
            {
                clStart = Convert.ToInt16(Math.Round(centre.C));
                rlStart = Convert.ToInt16(Math.Round(centre.R));
                slStart = Convert.ToInt16(Math.Round(centre.S));
            }
            for (int i = 0; i < numPoints; ++i)
            {
                int cA = clStart + i;
                int rB = rlStart + i;
                int sC = slStart + i;
                for (int j = 0; j < numPoints; ++j)
                {
                    int rA = rlStart + j;
                    int sB = slStart + j;
                    int cC = clStart + j;
                    for (int k = 0; k < numPoints; ++k)
                    {
                        int sA = slStart + k;
                        int cB = clStart + k;
                        int rC = rlStart + k;
                        if (isValidCRS(cA, rA, sA))
                        {
                            Double val = getCRSGridPoint(cA, rA, sA);
                            DensityPoint dp = new DensityPoint(cA, rA, sA, val, "CRS");
                            pointsA.Add(dp);
                        }
                        else
                        {
                            DensityPoint dp = new DensityPoint(cA, rA, sA, MagicInvalid, "CRS");
                            dp.valid = false;
                            pointsA.Add(dp);
                        }                        
                    }
                }
            }
            MatrixInterpolator mi = new MatrixInterpolator();
            List<int> diffs = new List<int>();
            diffs.Add(0);
            DensityPoint A = mi.getSplinedDensity(pointsA, centre,diffs);

            List<DensityPoint> finals = new List<DensityPoint>();

            finals.Add(A);
            
            bool isValid = true;
            double vAv = 0;
            double xAv = 0;
            double yAv = 0;
            double zAv = 0;
            double cAv = 0;
            double rAv = 0;
            double sAv = 0;

            Dictionary<int, double> vals = new Dictionary<int, double>();

            foreach (DensityPoint dp in finals)
            {
                isValid = dp.valid ? isValid : false;
                vAv += dp.V;
                xAv += dp.X;
                yAv += dp.Y;
                zAv += dp.Z;
                cAv += dp.C;
                rAv += dp.R;
                sAv += dp.S;
                
                foreach (var v in A.values)
                {
                    double vv = dp.values[v.Key]; // does a negative gradient mean anything for firsrt derivative from this perspective?                    
                    if (vals.ContainsKey(v.Key))
                        vals[v.Key] += vv;
                    else
                        vals[v.Key] = vv;
                }                
            }
            vAv /= finals.Count;
            xAv /= finals.Count;
            yAv /= finals.Count;
            zAv /= finals.Count;
            cAv /= finals.Count;
            rAv /= finals.Count;
            sAv /= finals.Count;
            foreach (var v in A.values)
            {
                vals[v.Key] /= finals.Count;
            }

            if (!isValid)
            {
                DensityPoint interped = new DensityPoint(cAv, rAv, sAv, -1000, "CRS");
                interped.valid = false;
                interped.values = A.values;
                return interped;
            }
            else
            {
                A.v = vAv;
                A.setXYZ(xAv, yAv, zAv);
                A.setCRS(cAv, rAv, sAv);
                A.values = vals;
                return A;
            }
        }



        private Double getFraction(DensityPoint centre, DensityPoint p1, DensityPoint p2)
        {//Using the cosine rule, this interpolating within the crs matrix
            Double fraction = 0;
            Double a = Math.Sqrt(Math.Pow(centre.C - p1.C, 2) + Math.Pow(centre.R - p1.R, 2) + Math.Pow(centre.S - p1.S, 2));
            Double b = Math.Sqrt(Math.Pow(centre.C - p2.C, 2) + Math.Pow(centre.R - p2.R, 2) + Math.Pow(centre.S - p2.S, 2));
            Double c = Math.Sqrt(Math.Pow(p1.C - p2.C, 2) + Math.Pow(p1.R - p2.R, 2) + Math.Pow(p1.S - p2.S, 2));
            if (c != 0)
            {
                Double x = (Math.Pow(a, 2) + Math.Pow(c, 2) - Math.Pow(b, 2)) / (2 * c);
                fraction = x / c;
            }
            return fraction;
        }

        /*
         *The following functionaility is based on the code in the python library
         * https://pdb-eda.readthedocs.io/en/latest/
         *
         */

        private void calculateOrthoMat()
        {
            _orthoMat = new EdMatrix();
            double alpha = Math.PI / 180 * w14_CELLB_X;
            double beta = Math.PI / 180 * w15_CELLB_Y;
            double gamma = Math.PI / 180 * w16_CELLB_Z;
            double temp = Math.Sqrt(1 - Math.Pow(Math.Cos(alpha), 2) - Math.Pow(Math.Cos(beta), 2) - Math.Pow(Math.Cos(gamma), 2) + 2 * Math.Cos(alpha) * Math.Cos(beta) * Math.Cos(gamma));
            _orthoMat.matrix[0, 0] = w11_CELLA_X;
            _orthoMat.matrix[0, 1] = w12_CELLA_Y * Math.Cos(gamma);
            _orthoMat.matrix[0, 2] = w13_CELLA_Z * Math.Cos(beta);
            _orthoMat.matrix[1, 0] = 0;
            _orthoMat.matrix[1, 1] = w12_CELLA_Y * Math.Sin(gamma);
            _orthoMat.matrix[1, 2] = w13_CELLA_Z * (Math.Cos(alpha) - Math.Cos(beta) * Math.Cos(gamma)) / Math.Sin(gamma);
            _orthoMat.matrix[2, 0] = 0;
            _orthoMat.matrix[2, 1] = 0;
            _orthoMat.matrix[2, 2] = w13_CELLA_Z * temp / Math.Sin(gamma);

            _deOrthoMat = _orthoMat.getInverse();

        }
        private void calculateOrigin()
        {
            // TODO I am ignoring the possibility of passing in the origin for now and using the dot product calc for non orthoganality.
            // The origin is perhaps used for cryoEM only and requires orthoganility

            // CRSSTART is w05_NXSTART, w06_NYSTART, w07_NZSTART             
            // Cell dims w08_MX, w09_MY, w10_MZ;            
            // Map of indices from crs to xyz is w17_MAPC, w18_MAPR, w19_MAPS

            EdVector vCRS = new EdVector();
            for (int i = 0; i < 3; ++i)
            {
                int startVal = 0;
                if (w17_MAPC == i)
                    startVal = w05_NXSTART;
                else if (w18_MAPR == i)
                    startVal = w06_NYSTART;
                else
                    startVal = w07_NZSTART;
             
                vCRS.vector[i] = startVal;
            }            
            vCRS.vector[0] /= w08_MX;
            vCRS.vector[1] /= w09_MY;
            vCRS.vector[2] /= w10_MZ;

            _origin = _orthoMat.multiply(vCRS);         
        }
        public DensityPoint getXYZFromCRS(double c, double r, double s)
        {            
            EdVector vXYZ = new EdVector();            
            EdVector vCRSIn = new EdVector();            
            vCRSIn.vector[0] = c;
            vCRSIn.vector[1] = r;
            vCRSIn.vector[2] = s;
            //If the axes are all orthogonal            
            if (w14_CELLB_X == 90 && w15_CELLB_Y == 90 && w16_CELLB_Z == 90)
            {
                for (int i = 0; i < 3; ++i)
                {                    
                    double startVal = vCRSIn.vector[_map2xyz[i]];
                    startVal *= _cellDims[i] / _axisSampling[i];                    
                    startVal += _origin.vector[i];
                    vXYZ.vector[i] = startVal;             
                }
            }
            else // they are not orthogonal
            {
                EdVector vCRS = new EdVector();
                for (int i = 0; i < 3; ++i)
                {
                    Double startVal = 0;
                    if (w17_MAPC == i)
                        startVal = w05_NXSTART + c;
                    else if (w18_MAPR == i)
                        startVal = w06_NYSTART + r;
                    else
                        startVal = w07_NZSTART + s;
                    
                    vCRS.vector[i] = startVal;
                }                
                vCRS.vector[0] /= w08_MX;
                vCRS.vector[1] /= w09_MY;
                vCRS.vector[2] /= w10_MZ;
                vXYZ = _orthoMat.multiply(vCRS);                
            }
            DensityPoint dp = new DensityPoint(Convert.ToSingle(vXYZ.vector[0]), Convert.ToSingle(vXYZ.vector[1]), Convert.ToSingle(vXYZ.vector[2]), 0, "XYZ");            
            return dp;
        }

        public DensityPoint getCRSFromXYZ(double x, double y, double z)
        {            
            EdVector vXYZIn = new EdVector();         
            EdVector vCRS = new EdVector();         
            vXYZIn.vector[0] = x;
            vXYZIn.vector[1] = y;
            vXYZIn.vector[2] = z;            
            //If the axes are all orthogonal            
            if (w14_CELLB_X == 90 && w15_CELLB_Y == 90 && w16_CELLB_Z == 90)
            {
                for (int i = 0; i < 3; ++i)
                {                    
                    double startVal = vXYZIn.vector[i] - _origin.vector[i];
                    startVal /= _cellDims[i] / _axisSampling[i];                 
                    vCRS.vector[i] = startVal;
                }
            }
            else // they are not orthogonal
            {             
                EdVector vFraction = _deOrthoMat.multiply(vXYZIn);
                for (int i = 0; i < 3; ++i)
                {             
                    double val = vFraction.vector[i] * _axisSampling[i] - _crsStart[_map2xyz[i]];             
                    vCRS.vector[i] = val;
                }
            }            
            double c = Convert.ToSingle(Math.Round(vCRS.vector[_map2crs[0]], 4));
            double r = Convert.ToSingle(Math.Round(vCRS.vector[_map2crs[1]], 4));
            double s = Convert.ToSingle(Math.Round(vCRS.vector[_map2crs[2]], 4));

            DensityPoint dp = new DensityPoint(c, r, s, 0, "CRS");
            return dp;                
        }
    }
}
