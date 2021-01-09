using System;

namespace ElectronDensity
{
    class Program
    {
        static void Main(string[] args)
        {            
            string ccp4Path = "F:\\Code\\ProteinDataFiles\\ccp4_data\\";
            Ccp4 ed = new Ccp4("1ejg", ccp4Path, false);            
            
            DensityPoint dp = ed.getCRSFromXYZ(0, 0, 0);
            Console.WriteLine(dp.C + "," +  dp.R + "," + dp.S);
            
            DensityPoint dp2 = ed.getCRSValue(dp.C, dp.R, dp.S);
            Console.WriteLine(dp2.V);

            DensityPoint dp3 = ed.getXYZFromCRS(dp.C, dp.R, dp.S);
            Console.WriteLine(dp3.X + "," + dp.Y + "," + dp.Z);

        }
    }
}
