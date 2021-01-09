using System;
using System.Collections.Generic;
using System.Text;

namespace ElectronDensity
{
    public class MatrixInterpolator
    {
        public DensityPoint getSplinedDensity(List<DensityPoint> points, DensityPoint centre, List<int> diffs) // #recursive function
        {
            /*
            RECURSIVE            
            List must be 2 ^ x long
            */
            int numPoints = 1 + 1;

            if (points.Count == 1) // # end of the recursion, return
            {
                foreach (int v in diffs)
                {
                    if (v == 0)
                        points[0].values[v] = points[0].v;
                    else
                        points[0].values[v] = 0;
                }

                return points[0];
            }
            else if (points.Count <= numPoints) // # end of the recursion, return
            {
                List<double> vs = new List<double>();
                bool valid = true;
                foreach (var point in points)
                {
                    double v = point.V;
                    if (point.valid)
                        vs.Add(v);
                    else
                        valid = false;
                }
                int half = Convert.ToInt32(points.Count / 2);
                DensityPoint p1 = points[half - 1];
                DensityPoint p2 = points[half];
                double fr = getFractionCosineRule(centre, p1, p2);
                double x = p1.C + fr * (p2.C - p1.C);
                double y = p1.R + fr * (p2.R - p1.R);
                double z = p1.S + fr * (p2.S - p1.S);
                if (vs.Count == 0 || !valid)
                {
                    DensityPoint interped = new DensityPoint(x, y, z, -1000, "CRS");
                    foreach (int v in diffs)
                    {
                        interped.values[v] = -1000;
                        interped.valid = false;
                    }
                    return interped;
                }
                else
                {
                    Polynomial poly = new Polynomial(vs);
                    //the poly is a sequence so the value we want is a fraction along from the halfway markers
                    double valPoint = half + fr;
                    List<double> finalvs = poly.getValue(valPoint, diffs);
                    if (finalvs.Count > 0)
                    {
                        DensityPoint interped = new DensityPoint(x, y, z, finalvs[0], "CRS");
                        for (int i = 0; i < diffs.Count; ++i)
                        {
                            interped.values[diffs[i]] = finalvs[i];
                        }
                        return interped;
                    }
                    else
                    {
                        DensityPoint interped = new DensityPoint(x, y, z, -1000, "CRS");
                        foreach (int v in diffs)
                        {
                            interped.values[v] = -1000;
                            interped.valid = false;
                        }
                        return interped;
                    }
                }
            }
            else//#split recursion down further
            {
                List<DensityPoint> ps = new List<DensityPoint>();
                List<DensityPoint> tmpps = new List<DensityPoint>();
                for (int i = 0; i < points.Count; ++i)
                {
                    tmpps.Add(points[i]);
                    if (tmpps.Count == numPoints)
                    {
                        DensityPoint newA = getSplinedDensity(tmpps, centre, diffs);
                        ps.Add(newA);
                        tmpps.Clear();
                    }
                }
                return getSplinedDensity(ps, centre, diffs);
            }
        }

        private double getFractionCosineRule(DensityPoint centre, DensityPoint p1, DensityPoint p2)
        {
            //The angle beta is found from the cosine rule
            //cos beta  equates x/a to (a^2 + c^2 - b^2) / 2ac
            double a2 = Math.Pow((centre.C - p1.C), 2) + Math.Pow((centre.R - p1.R), 2) + Math.Pow((centre.S - p1.S), 2);
            double b2 = Math.Pow((centre.C - p2.C), 2) + Math.Pow((centre.R - p2.R), 2) + Math.Pow((centre.S - p2.S), 2);
            double c2 = Math.Pow((p1.C - p2.C), 2) + Math.Pow((p1.R - p2.R), 2) + Math.Pow((p1.S - p2.S), 2);
            double a = Math.Sqrt(a2);
            double b = Math.Sqrt(b2);
            double c = Math.Sqrt(c2);
            double fraction = 0;
            if (c != 0 && a != 0)
            {
                double cosBeta = (a2 + c2 - b2) / (2 * a * c);
                double length = cosBeta * a;
                fraction = length / c;
            }
            return fraction;
        }
    }
    }
