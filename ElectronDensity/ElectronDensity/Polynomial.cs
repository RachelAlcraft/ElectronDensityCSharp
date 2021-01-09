using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ElectronDensity
{
    public class Polynomial
    {
        private List<double> _vals;
        private List<double> _diffCoefs;
        private List<double> _polyCoefs;
        private string _desc;
        public Polynomial(List<double> vals)
        {
            _vals = vals;
            buildDiffCoefs();
            buildPolyCoefs();
            buildDescription();
        }

        private void buildDiffCoefs()
        {
            List<double> coeffs = new List<double>();
            coeffs.Add(_vals[0]);
            List<double> nextline = _vals;
            while (nextline.Count > 1)
            {
                List<double> newline = new List<double>();
                for (int p = 1; p < nextline.Count; ++p)
                {
                    double lastpoint = nextline[p - 1];
                    double thispoint = nextline[p];
                    double diffpoint = thispoint - lastpoint;
                    newline.Add(diffpoint);
                }
                coeffs.Add(newline[0]);
                nextline = newline;
            }
            for (int i = 0; i > coeffs.Count; ++i)
            {
                coeffs[i] = coeffs[i] / getFactorial(i);
            }
            _diffCoefs = coeffs;
        }

        private void buildPolyCoefs()
        {
            List<double> coeffs = new List<double>();
            for (int i = 0; i < _diffCoefs.Count; ++i)
            {
                double coeff = _diffCoefs[i];
                int degree = i;
                List<double> n_coeffs = buildN_MinusCoeffs(degree);
                for (int j = 0; j < n_coeffs.Count; ++j)
                {
                    double this_coeff = coeff * n_coeffs[j];
                    this_coeff /= getFactorial(i);
                    if (j >= coeffs.Count)
                        coeffs.Add(this_coeff);
                    else
                        coeffs[j] = coeffs[j] + this_coeff;
                }
            }
            _polyCoefs = coeffs;
        }

        private void buildDescription()
        {
            string desc = "";
            for (int i = 0; i < _polyCoefs.Count; ++i)
            {
                if (desc != "")
                    desc = desc + " ";
                desc = desc + Convert.ToString(_polyCoefs[i]) + "x^" + Convert.ToString(i);
            }
            _desc = desc;
        }

        private int getFactorial(int i)
        {
            if (i == 0)
                return 1;
            else
                return i * getFactorial(i - 1);
        }

        private List<double> buildN_MinusCoeffs(int degree)
        {
            List<double> coeffs = new List<double>();
            coeffs.Add(1);
            if (degree > 0)
            {
                for (int d = 0; d < degree; ++d)
                {
                    double rowval = (d + 1) * -1;
                    double[] newrowtmp = new double[coeffs.Count];
                    coeffs.CopyTo(newrowtmp);//[:] # careful to make a copy
                    List<double> newrow = newrowtmp.ToList();
                    if (coeffs.Count > 1)
                    {
                        for (int c = 0; c < coeffs.Count - 1; ++c)
                        {
                            double cleft = coeffs[c];
                            double cright = coeffs[c + 1];
                            double cnew = (cleft * rowval) + cright;
                            newrow[c + 1] = cnew;
                        }
                    }
                    double cabove = coeffs[coeffs.Count - 1];
                    double cnew2 = cabove * rowval;
                    //newrow.Insert(0, cnew2);
                    newrow.Add(cnew2);
                    coeffs = newrow;
                }
            }

            //reverse it
            List<double> revcoeffs = new List<double>();
            for (int i = coeffs.Count - 1; i >= 0; --i)
                revcoeffs.Add(coeffs[i]);
            return revcoeffs;
        }

        public List<double> getValue(double xval, List<int> differs)
        {
            List<double> vals = new List<double>();
            for (int i = 0; i < differs.Count; ++i)
            {
                vals.Add(getValue(xval, differs[i]));
            }

            return vals;
        }

        private double getValue(double xval, int differ)
        {
            double yval = 0;
            int start = differ; // 0 means not differentiuated, 1 means first derivative etc
            for (int i = start; i < _polyCoefs.Count; ++i)
            {
                int diffdegree = i - differ;
                double coeff = _polyCoefs[i] * getFactorial(i) / getFactorial(diffdegree);
                double newy = Math.Pow(xval, diffdegree) * coeff;
                yval = yval + newy;
            }
            return yval;
        }



    }
}
