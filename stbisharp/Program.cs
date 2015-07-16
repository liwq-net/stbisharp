using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

class Program
{
    unsafe static void Main(string[] args)
    {
        var buffer = File.ReadAllBytes("F:\\p.png");

        fixed (byte* p = buffer)
        {
            int x; int y; int comp;
            byte* result = liwq.stbisharp.stbi_load_from_memory(p, buffer.Length, &x, &y, &comp, 0);
        }
    }
}
