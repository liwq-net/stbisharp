using System;
using System.Diagnostics;

namespace liwq
{
    unsafe public class stbisharp
    {
        /// <summary>重载 read skip eof，并且重写函数</summary>
        public struct stbi_io_callbacks
        {
            public bool read { get; set; }
            public bool skip { get; set; }
            public bool eof { get; set; }

            public int _read(void* user, sbyte* data, int size)    // fill 'data' with 'size' bytes.  return number of bytes actually read
            {
                throw new NotImplementedException();
            }
            public void _skip(void* user, int n)                   // skip the next 'n' bytes, or 'unget' the last -n bytes if negative
            {
                throw new NotImplementedException();
            }
            public int _eof(void* user)                            // returns nonzero if we are at end of file/data
            {
                throw new NotImplementedException();
            }
        }

        ///////////////////////////////////////////////
        //
        // stbi__context struct and start_xxx functions
        // stbi__context structure is our basic context used by all images, so it
        // contains all the IO context, plus some basic image information

        public struct stbi__context
        {
            public uint img_x;
            public uint img_y;
            public int img_n;
            public int img_out_n;
            public stbi_io_callbacks io;
            public void* io_user_data;
            public int read_from_callbacks;
            public int buflen;
            public fixed byte buffer_start[128];
            public byte* img_buffer;
            public byte* img_buffer_end;
            public byte* img_buffer_original;
        }

        // initialize a memory-decode context
        static void stbi__start_mem(stbi__context* s, byte* buffer, int len)
        {
            s->io.read = false;
            s->read_from_callbacks = 0;
            s->img_buffer = s->img_buffer_original = (byte*)buffer;
            s->img_buffer_end = (byte*)buffer + len;
        }


        static int stbi__vertically_flip_on_load = 0;
        public static void stbi_set_flip_vertically_on_load(int flag_true_if_should_flip)
        {
            stbi__vertically_flip_on_load = flag_true_if_should_flip;
        }

        static byte* stbi__load_main(stbi__context* s, int* x, int* y, int* comp, int req_comp)
        {
            //#if !STBI_NO_JPEG
            //            if (stbi__jpeg_test(s)) return stbi__jpeg_load(s, x, y, comp, req_comp);
            //#endif
            //#if !STBI_NO_PNG
            if (stbi__png_test(s) != 0) return stbi__png_load(s, x, y, comp, req_comp);
            //#endif
            //#if !STBI_NO_BMP
            //            if (stbi__bmp_test(s)) return stbi__bmp_load(s, x, y, comp, req_comp);
            //#endif
            //#if !STBI_NO_GIF
            //            if (stbi__gif_test(s)) return stbi__gif_load(s, x, y, comp, req_comp);
            //#endif
            //#if !STBI_NO_PSD
            //            if (stbi__psd_test(s)) return stbi__psd_load(s, x, y, comp, req_comp);
            //#endif
            //#if !STBI_NO_PIC
            //            if (stbi__pic_test(s)) return stbi__pic_load(s, x, y, comp, req_comp);
            //#endif
            //#if !STBI_NO_PNM
            //            if (stbi__pnm_test(s)) return stbi__pnm_load(s, x, y, comp, req_comp);
            //#endif

            //#if !STBI_NO_HDR
            //            if (stbi__hdr_test(s))
            //            {
            //                float* hdr = stbi__hdr_load(s, x, y, comp, req_comp);
            //                return stbi__hdr_to_ldr(hdr, *x, *y, req_comp ? req_comp : *comp);
            //            }
            //#endif

            //#if !STBI_NO_TGA
            //            // test tga last because it's a crappy test!
            //            if (stbi__tga_test(s))
            //                return stbi__tga_load(s, x, y, comp, req_comp);
            //#endif

            throw new Exception("unknown image type:Image not of any known type, or corrupt");
        }

        static byte* stbi__load_flip(stbi__context* s, int* x, int* y, int* comp, int req_comp)
        {
            byte* result = stbi__load_main(s, x, y, comp, req_comp);
            if (stbi__vertically_flip_on_load != 0 && result != null)
            {
                int w = *x, h = *y;
                int depth = req_comp != 0 ? req_comp : *comp;
                int row, col, z;
                byte temp;
                // @OPTIMIZE: use a bigger temp buffer and memcpy multiple pixels at once
                for (row = 0; row < (h >> 1); row++)
                {
                    for (col = 0; col < w; col++)
                    {
                        for (z = 0; z < depth; z++)
                        {
                            temp = result[(row * w + col) * depth + z];
                            result[(row * w + col) * depth + z] = result[((h - row - 1) * w + col) * depth + z];
                            result[((h - row - 1) * w + col) * depth + z] = temp;
                        }
                    }
                }
            }
            return result;
        }


        public static byte* stbi_load_from_memory(byte* buffer, int len, int* x, int* y, int* comp, int req_comp)
        {
            stbi__context s;
            stbi__start_mem(&s, buffer, len);
            return stbi__load_flip(&s, x, y, comp, req_comp);
        }
        //////////////////////////////////////////////////////////////////////////////
        //
        // Common code used by all image loaders
        //

        public const int STBI__SCAN_load = 0;
        public const int STBI__SCAN_type = 1;
        public const int STBI__SCAN_header = 2;

        static void stbi__refill_buffer(stbi__context* s)
        {
            int n = s->io._read(s->io_user_data, (sbyte*)s->buffer_start, s->buflen);
            if (n == 0)
            {
                // at end of file, treat same as if from memory, but need to handle case
                // where s->img_buffer isn't pointing to safe memory, e.g. 0-byte file
                s->read_from_callbacks = 0;
                s->img_buffer = s->buffer_start;
                s->img_buffer_end = s->buffer_start + 1;
                *s->img_buffer = 0;
            }
            else
            {
                s->img_buffer = s->buffer_start;
                s->img_buffer_end = s->buffer_start + n;
            }
        }

        static byte stbi__get8(stbi__context* s)
        {
            if (s->img_buffer < s->img_buffer_end)
                return *s->img_buffer++;
            if (s->read_from_callbacks != 0)
            {
                stbi__refill_buffer(s);
                return *s->img_buffer++;
            }
            return 0;
        }

        static int stbi__at_eof(stbi__context* s)
        {
            if (s->io.read)
            {
                if ((s->io._eof)(s->io_user_data) == 0)
                    return 0;
                // if feof() is true, check if buffer = end
                // special case: we've only got the special 0 character at the end
                if (s->read_from_callbacks == 0) return 1;
            }

            return s->img_buffer >= s->img_buffer_end ? 1 : 0;
        }

        static void stbi__skip(stbi__context* s, int n)
        {
            if (n < 0)
            {
                s->img_buffer = s->img_buffer_end;
                return;
            }
            if (s->io.read)
            {
                int blen = (int)(s->img_buffer_end - s->img_buffer);
                if (blen < n)
                {
                    s->img_buffer = s->img_buffer_end;
                    (s->io._skip)(s->io_user_data, n - blen);
                    return;
                }
            }
            s->img_buffer += n;
        }

        static int stbi__getn(stbi__context* s, byte* buffer, int n)
        {
            if (s->io.read)
            {
                int blen = (int)(s->img_buffer_end - s->img_buffer);
                if (blen < n)
                {
                    int res, count;

                    CLib.CString.memcpy(buffer, s->img_buffer, (uint)blen);

                    count = s->io._read(s->io_user_data, (sbyte*)buffer + blen, n - blen);
                    res = (count == (n - blen)) ? 1 : 0;
                    s->img_buffer = s->img_buffer_end;
                    return res;
                }
            }

            if (s->img_buffer + n <= s->img_buffer_end)
            {
                CLib.CString.memcpy(buffer, s->img_buffer, (uint)n);
                s->img_buffer += n;
                return 1;
            }
            else
                return 0;
        }

        static int stbi__get16be(stbi__context* s)
        {
            int z = stbi__get8(s);
            return (z << 8) + stbi__get8(s);
        }

        static uint stbi__get32be(stbi__context* s)
        {
            uint z = (uint)stbi__get16be(s);
            return (z << 16) + (uint)stbi__get16be(s);
        }

        static int stbi__get16le(stbi__context* s)
        {
            int z = stbi__get8(s);
            return z + (stbi__get8(s) << 8);
        }

        static uint stbi__get32le(stbi__context* s)
        {
            uint z = (uint)stbi__get16le(s);
            return z + (uint)(stbi__get16le(s) << 16);
        }

        static void stbi__rewind(stbi__context* s)
        {
            // conceptually rewind SHOULD rewind to the beginning of the stream,
            // but we just rewind to the beginning of the initial buffer, because
            // we only use it after doing 'test', which only ever looks at at most 92 bytes
            s->img_buffer = s->img_buffer_original;
        }

        static byte stbi__compute_y(int r, int g, int b)
        {
            return (byte)(((r * 77) + (g * 150) + (29 * b)) >> 8);
        }

        static byte* stbi__convert_format(byte* data, int img_n, int req_comp, uint x, uint y)
        {
            int i, j;
            byte* good;

            if (req_comp == img_n) return data;
            Debug.Assert(req_comp >= 1 && req_comp <= 4);

            good = (byte*)CLib.CStdlib.malloc((int)(req_comp * x * y));
            if (good == null)
            {
                CLib.CStdlib.free(data);
                throw new Exception("outofmem:Out of memory");
            }

            for (j = 0; j < (int)y; ++j)
            {
                byte* src = data + j * x * img_n;
                byte* dest = good + j * x * req_comp;
                // convert source image with img_n components to one with req_comp components;
                // avoid switch per pixel, so use switch per scanline and massive macros
                switch (img_n * 8 + req_comp)
                {
                    case (1 * 8 + 2): for (i = (int)x - 1; i >= 0; --i, src += 1, dest += 2) dest[0] = src[0]; dest[1] = 255; break;
                    case (1 * 8 + 3): for (i = (int)x - 1; i >= 0; --i, src += 1, dest += 3) dest[0] = dest[1] = dest[2] = src[0]; break;
                    case (1 * 8 + 4): for (i = (int)x - 1; i >= 0; --i, src += 1, dest += 4) dest[0] = dest[1] = dest[2] = src[0]; dest[3] = 255; break;
                    case (2 * 8 + 1): for (i = (int)x - 1; i >= 0; --i, src += 2, dest += 1) dest[0] = src[0]; break;
                    case (2 * 8 + 3): for (i = (int)x - 1; i >= 0; --i, src += 2, dest += 3) dest[0] = dest[1] = dest[2] = src[0]; break;
                    case (2 * 8 + 4): for (i = (int)x - 1; i >= 0; --i, src += 2, dest += 4) dest[0] = dest[1] = dest[2] = src[0]; dest[3] = src[1]; break;
                    case (3 * 8 + 4): for (i = (int)x - 1; i >= 0; --i, src += 3, dest += 4) dest[0] = src[0]; dest[1] = src[1]; dest[2] = src[2]; dest[3] = 255; break;
                    case (3 * 8 + 1): for (i = (int)x - 1; i >= 0; --i, src += 3, dest += 1) dest[0] = stbi__compute_y(src[0], src[1], src[2]); break;
                    case (3 * 8 + 2): for (i = (int)x - 1; i >= 0; --i, src += 3, dest += 2) dest[0] = stbi__compute_y(src[0], src[1], src[2]); dest[1] = 255; break;
                    case (4 * 8 + 1): for (i = (int)x - 1; i >= 0; --i, src += 4, dest += 1) dest[0] = stbi__compute_y(src[0], src[1], src[2]); break;
                    case (4 * 8 + 2): for (i = (int)x - 1; i >= 0; --i, src += 4, dest += 2) dest[0] = stbi__compute_y(src[0], src[1], src[2]); dest[1] = src[3]; break;
                    case (4 * 8 + 3): for (i = (int)x - 1; i >= 0; --i, src += 4, dest += 3) dest[0] = src[0]; dest[1] = src[1]; dest[2] = src[2]; break;
                    default: Debug.Assert(false); break;
                }
            }
            CLib.CStdlib.free(data);
            return good;
        }

        /// <summary>初始化静态数组</summary>
        static stbisharp()
        {
            stbi__zdefault_length = (byte*)CLib.CStdlib.malloc(288);
            stbi__zdefault_distance = (byte*)CLib.CStdlib.malloc(32);
        }

        #region zlib

#if !STBI_NO_ZLIB

        // fast-way is faster to check than jpeg huffman, but slow way is slower
        const int STBI__ZFAST_BITS = 9; // accelerate all cases in default tables
        const int STBI__ZFAST_MASK = ((1 << STBI__ZFAST_BITS) - 1);

        // zlib-style huffman encoding
        // (jpegs packs from left, zlib from right, so can't share code)
        public struct stbi__zhuffman
        {
            public fixed ushort fast[1 << STBI__ZFAST_BITS];
            public fixed ushort firstcode[16];
            public fixed int maxcode[17];
            public fixed ushort firstsymbol[16];
            public fixed byte size[288];
            public fixed ushort value[288];
        }

        static int stbi__bitreverse16(int n)
        {
            n = ((n & 0xAAAA) >> 1) | ((n & 0x5555) << 1);
            n = ((n & 0xCCCC) >> 2) | ((n & 0x3333) << 2);
            n = ((n & 0xF0F0) >> 4) | ((n & 0x0F0F) << 4);
            n = ((n & 0xFF00) >> 8) | ((n & 0x00FF) << 8);
            return n;
        }

        static int stbi__bit_reverse(int v, int bits)
        {
            Debug.Assert(bits <= 16);
            // to bit reverse n bits, reverse 16 and shift
            // e.g. 11 bits, bit reverse and shift away 5
            return stbi__bitreverse16(v) >> (16 - bits);
        }

        static int stbi__zbuild_huffman(stbi__zhuffman* z, byte* sizelist, int num)
        {
            int i; int k = 0;
            int code;
            int* next_code = stackalloc int[16];
            int* sizes = stackalloc int[17];

            // DEFLATE spec for generating codes
            CLib.CString.memset(sizes, 0, 17 * 4);
            CLib.CString.memset(z->fast, 0, (1 << STBI__ZFAST_BITS) * 2);

            for (i = 0; i < num; ++i)
                ++sizes[sizelist[i]];
            sizes[0] = 0;
            for (i = 1; i < 16; ++i)
                if (sizes[i] > (1 << i))
                    throw new Exception("bad sizes:Corrupt PNG");
            code = 0;
            for (i = 1; i < 16; ++i)
            {
                next_code[i] = code;
                z->firstcode[i] = (ushort)code;
                z->firstsymbol[i] = (ushort)k;
                code = (code + sizes[i]);
                if (sizes[i] != 0)
                    if (code - 1 >= (1 << i))
                        throw new Exception("bad codelengths:Corrupt PNG");
                z->maxcode[i] = code << (16 - i); // preshift for inner loop
                code <<= 1;
                k += sizes[i];
            }
            z->maxcode[16] = 0x10000; // sentinel
            for (i = 0; i < num; ++i)
            {
                int s = sizelist[i];
                if (s != 0)
                {
                    int c = next_code[s] - z->firstcode[s] + z->firstsymbol[s];
                    ushort fastv = (ushort)((s << 9) | i);
                    z->size[c] = (byte)s;
                    z->value[c] = (ushort)i;
                    if (s <= STBI__ZFAST_BITS)
                    {
                        int kk = stbi__bit_reverse(next_code[s], s);
                        while (kk < (1 << STBI__ZFAST_BITS))
                        {
                            z->fast[kk] = fastv;
                            kk += (1 << s);
                        }
                    }
                    ++next_code[s];
                }
            }
            return 1;
        }

        // zlib-from-memory implementation for PNG reading
        // because PNG allows splitting the zlib stream arbitrarily,
        // and it's annoying structurally to have PNG call ZLIB call PNG,
        // we require PNG read all the IDATs and combine them into a single
        // memory buffer
        public struct stbi__zbuf
        {
            public byte* zbuffer;
            public byte* zbuffer_end;
            public int num_bits;
            public uint code_buffer;
            public sbyte* zout;
            public sbyte* zout_start;
            public sbyte* zout_end;
            public int z_expandable;
            public stbi__zhuffman z_length;
            public stbi__zhuffman z_distance;
        }

        static byte stbi__zget8(stbi__zbuf* z)
        {
            if (z->zbuffer >= z->zbuffer_end) return 0;
            return *z->zbuffer++;
        }

        static void stbi__fill_bits(stbi__zbuf* z)
        {
            do
            {
                Debug.Assert(z->code_buffer < (1U << z->num_bits));
                z->code_buffer |= (uint)((int)stbi__zget8(z) << (int)z->num_bits);
                z->num_bits += 8;
            } while (z->num_bits <= 24);
        }

        static uint stbi__zreceive(stbi__zbuf* z, int n)
        {
            uint k;
            if (z->num_bits < n) stbi__fill_bits(z);
            k = (uint)((int)z->code_buffer & ((1 << n) - 1));
            z->code_buffer >>= n;
            z->num_bits -= n;
            return k;
        }

        static int stbi__zhuffman_decode_slowpath(stbi__zbuf* a, stbi__zhuffman* z)
        {
            int b, s, k;
            // not resolved by fast table, so compute it the slow way
            // use jpeg approach, which requires MSbits at top
            k = stbi__bit_reverse((int)a->code_buffer, 16);
            for (s = STBI__ZFAST_BITS + 1; ; ++s)
                if (k < z->maxcode[s])
                    break;
            if (s == 16) return -1; // invalid code!
            // code size is s, so:
            b = (k >> (16 - s)) - z->firstcode[s] + z->firstsymbol[s];
            Debug.Assert(z->size[b] == s);
            a->code_buffer >>= s;
            a->num_bits -= s;
            return z->value[b];
        }

        static int stbi__zhuffman_decode(stbi__zbuf* a, stbi__zhuffman* z)
        {
            int b, s;
            if (a->num_bits < 16) stbi__fill_bits(a);
            b = z->fast[a->code_buffer & STBI__ZFAST_MASK];
            if (b != 0)
            {
                s = b >> 9;
                a->code_buffer >>= s;
                a->num_bits -= s;
                return b & 511;
            }
            return stbi__zhuffman_decode_slowpath(a, z);
        }

        static int stbi__zexpand(stbi__zbuf* z, sbyte* zout, int n)  // need to make room for n bytes
        {
            sbyte* q;
            int cur, limit;
            z->zout = zout;
            if (z->z_expandable == 0) throw new Exception("output buffer limit:Corrupt PNG");
            cur = (int)(z->zout - z->zout_start);
            limit = (int)(z->zout_end - z->zout_start);
            while (cur + n > limit)
                limit *= 2;
            q = (sbyte*)CLib.CStdlib.realloc(z->zout_start, limit);
            if (q == null) throw new Exception("outofmem:Out of memory");
            z->zout_start = q;
            z->zout = q + cur;
            z->zout_end = q + limit;
            return 1;
        }

        static int[] stbi__zlength_base = new int[31] { 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258, 0, 0 };
        static int[] stbi__zlength_extra = new int[31] { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0, 0, 0 };
        static int[] stbi__zdist_base = new int[32] { 1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577, 0, 0 };
        static int[] stbi__zdist_extra = new int[32] { 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 0, 0 };

        static int stbi__parse_huffman_block(stbi__zbuf* a)
        {
            sbyte* zout = a->zout;
            for (; ; )
            {
                int z = stbi__zhuffman_decode(a, &a->z_length);
                if (z < 256)
                {
                    if (z < 0) throw new Exception("bad huffman code:Corrupt PNG"); // error in huffman codes
                    if (zout >= a->zout_end)
                    {
                        if (stbi__zexpand(a, zout, 1) == 0) return 0;
                        zout = a->zout;
                    }
                    *zout++ = (sbyte)z;
                }
                else
                {
                    byte* p;
                    int len, dist;
                    if (z == 256)
                    {
                        a->zout = zout;
                        return 1;
                    }
                    z -= 257;
                    len = stbi__zlength_base[z];
                    if (stbi__zlength_extra[z] != 0) len += (int)stbi__zreceive(a, stbi__zlength_extra[z]);
                    z = stbi__zhuffman_decode(a, &a->z_distance);
                    if (z < 0) throw new Exception("bad huffman code:Corrupt PNG");
                    dist = stbi__zdist_base[z];
                    if (stbi__zdist_extra[z] != 0) dist += (int)stbi__zreceive(a, stbi__zdist_extra[z]);
                    if (zout - a->zout_start < dist) throw new Exception("bad dist:Corrupt PNG");
                    if (zout + len > a->zout_end)
                    {
                        if (stbi__zexpand(a, zout, len) == 0) return 0;
                        zout = a->zout;
                    }
                    p = (byte*)(zout - dist);
                    if (dist == 1)
                    { // run of one byte; common in images.
                        byte v = *p;
                        if (len != 0) { do *zout++ = (sbyte)v; while ((--len) != 0); }
                    }
                    else
                    {
                        if (len != 0) { do *zout++ = (sbyte)*p++; while ((--len) != 0); }
                    }
                }
            }
        }

        static int stbi__compute_huffman_codes(stbi__zbuf* a)
        {
            byte* length_dezigzag = stackalloc byte[19];
            length_dezigzag[0] = 16;
            length_dezigzag[1] = 17;
            length_dezigzag[2] = 18;
            length_dezigzag[3] = 0;
            length_dezigzag[4] = 8;
            length_dezigzag[5] = 7;
            length_dezigzag[6] = 9;
            length_dezigzag[7] = 6;
            length_dezigzag[8] = 10;
            length_dezigzag[9] = 5;
            length_dezigzag[10] = 11;
            length_dezigzag[11] = 4;
            length_dezigzag[12] = 12;
            length_dezigzag[13] = 3;
            length_dezigzag[14] = 13;
            length_dezigzag[15] = 2;
            length_dezigzag[16] = 14;
            length_dezigzag[17] = 1;
            length_dezigzag[18] = 15;

            stbi__zhuffman z_codelength;
            byte* lencodes = stackalloc byte[286 + 32 + 137];//padding for maximum single op
            byte* codelength_sizes = stackalloc byte[19];
            int i, n;

            int hlit = (int)stbi__zreceive(a, 5) + 257;
            int hdist = (int)stbi__zreceive(a, 5) + 1;
            int hclen = (int)stbi__zreceive(a, 4) + 4;

            CLib.CString.memset(codelength_sizes, 0, 19);
            for (i = 0; i < hclen; ++i)
            {
                int s = (int)stbi__zreceive(a, 3);
                codelength_sizes[length_dezigzag[i]] = (byte)s;
            }
            if (stbi__zbuild_huffman(&z_codelength, codelength_sizes, 19) == 0) return 0;

            n = 0;
            while (n < hlit + hdist)
            {
                int c = stbi__zhuffman_decode(a, &z_codelength);
                if (c < 0 || c >= 19) throw new Exception("bad codelengths:Corrupt PNG");
                if (c < 16)
                    lencodes[n++] = (byte)c;
                else if (c == 16)
                {
                    c = (int)stbi__zreceive(a, 2) + 3;
                    CLib.CString.memset(lencodes + n, lencodes[n - 1], c);
                    n += c;
                }
                else if (c == 17)
                {
                    c = (int)stbi__zreceive(a, 3) + 3;
                    CLib.CString.memset(lencodes + n, 0, c);
                    n += c;
                }
                else
                {
                    Debug.Assert(c == 18);
                    c = (int)stbi__zreceive(a, 7) + 11;
                    CLib.CString.memset(lencodes + n, 0, c);
                    n += c;
                }
            }
            if (n != hlit + hdist) throw new Exception("bad codelengths:Corrupt PNG");
            if (stbi__zbuild_huffman(&a->z_length, lencodes, hlit) == 0) return 0;
            if (stbi__zbuild_huffman(&a->z_distance, lencodes + hlit, hdist) == 0) return 0;
            return 1;
        }

        static int stbi__parse_uncomperssed_block(stbi__zbuf* a)
        {
            byte* header = stackalloc byte[4];
            int len, nlen, k;
            if ((a->num_bits & 7) != 0)
                stbi__zreceive(a, a->num_bits & 7); // discard
            // drain the bit-packed data into header
            k = 0;
            while (a->num_bits > 0)
            {
                header[k++] = (byte)(a->code_buffer & 255); // suppress MSVC run-time check
                a->code_buffer >>= 8;
                a->num_bits -= 8;
            }
            Debug.Assert(a->num_bits == 0);
            // now fill header the normal way
            while (k < 4)
                header[k++] = stbi__zget8(a);
            len = header[1] * 256 + header[0];
            nlen = header[3] * 256 + header[2];
            if (nlen != (len ^ 0xffff)) throw new Exception("zlib corrupt:Corrupt PNG");
            if (a->zbuffer + len > a->zbuffer_end) throw new Exception("read past buffer:Corrupt PNG");
            if (a->zout + len > a->zout_end)
                if (stbi__zexpand(a, a->zout, len) == 0) return 0;
            CLib.CString.memcpy(a->zout, a->zbuffer, (uint)len);
            a->zbuffer += len;
            a->zout += len;
            return 1;
        }

        static int stbi__parse_zlib_header(stbi__zbuf* a)
        {
            int cmf = stbi__zget8(a);
            int cm = cmf & 15;
            /* int cinfo = cmf >> 4; */
            int flg = stbi__zget8(a);
            if ((cmf * 256 + flg) % 31 != 0) throw new Exception("bad zlib header:Corrupt PNG"); // zlib spec
            if ((flg & 32) != 0) throw new Exception("no preset dict:Corrupt PNG"); // preset dictionary not allowed in png
            if (cm != 8) throw new Exception("bad compression:Corrupt PNG"); // DEFLATE required for png
            // window = 1 << (8 + cinfo)... but who cares, we fully buffer output
            return 1;
        }

        // @TODO: should statically initialize these for optimal thread safety
        static byte* stbi__zdefault_length;
        static byte* stbi__zdefault_distance;

        static void stbi__init_zdefaults()
        {
            int i;   // use <= to match clearly with spec
            for (i = 0; i <= 143; ++i) stbi__zdefault_length[i] = 8;
            for (; i <= 255; ++i) stbi__zdefault_length[i] = 9;
            for (; i <= 279; ++i) stbi__zdefault_length[i] = 7;
            for (; i <= 287; ++i) stbi__zdefault_length[i] = 8;
            for (i = 0; i <= 31; ++i) stbi__zdefault_distance[i] = 5;
        }

        static int stbi__parse_zlib(stbi__zbuf* a, int parse_header)
        {
            int final, type;
            if (parse_header != 0)
                if (stbi__parse_zlib_header(a) == 0) return 0;
            a->num_bits = 0;
            a->code_buffer = 0;
            do
            {
                final = (int)stbi__zreceive(a, 1);
                type = (int)stbi__zreceive(a, 2);
                if (type == 0)
                {
                    if (stbi__parse_uncomperssed_block(a) == 0) return 0;
                }
                else if (type == 3)
                {
                    return 0;
                }
                else
                {
                    if (type == 1)
                    {
                        // use fixed code lengths
                        if (stbi__zdefault_distance[31] == 0) stbi__init_zdefaults();
                        if (stbi__zbuild_huffman(&a->z_length, stbi__zdefault_length, 288) == 0) return 0;
                        if (stbi__zbuild_huffman(&a->z_distance, stbi__zdefault_distance, 32) == 0) return 0;
                    }
                    else
                    {
                        if (stbi__compute_huffman_codes(a) == 0) return 0;
                    }
                    if (stbi__parse_huffman_block(a) == 0) return 0;
                }
            } while (final == 0);
            return 1;
        }

        static int stbi__do_zlib(stbi__zbuf* a, sbyte* obuf, int olen, int exp, int parse_header)
        {
            a->zout_start = obuf;
            a->zout = obuf;
            a->zout_end = obuf + olen;
            a->z_expandable = exp;
            return stbi__parse_zlib(a, parse_header);
        }

        public static sbyte* stbi_zlib_decode_malloc_guesssize(sbyte* buffer, int len, int initial_size, int* outlen)
        {
            stbi__zbuf a;
            sbyte* p = (sbyte*)CLib.CStdlib.malloc(initial_size);
            if (p == null) return null;
            a.zbuffer = (byte*)buffer;
            a.zbuffer_end = (byte*)buffer + len;
            if (stbi__do_zlib(&a, p, initial_size, 1, 1) != 0)
            {
                if (outlen != null) *outlen = (int)(a.zout - a.zout_start);
                return a.zout_start;
            }
            else
            {
                CLib.CStdlib.free(a.zout_start);
                return null;
            }
        }

        public static sbyte* stbi_zlib_decode_malloc(sbyte* buffer, int len, int* outlen)
        {
            return stbi_zlib_decode_malloc_guesssize(buffer, len, 16384, outlen);
        }

        public static sbyte* stbi_zlib_decode_malloc_guesssize_headerflag(sbyte* buffer, int len, int initial_size, int* outlen, int parse_header)
        {
            stbi__zbuf a;
            sbyte* p = (sbyte*)CLib.CStdlib.malloc(initial_size);
            if (p == null) return null;
            a.zbuffer = (byte*)buffer;
            a.zbuffer_end = (byte*)buffer + len;
            if (stbi__do_zlib(&a, p, initial_size, 1, parse_header) != 0)
            {
                if (outlen != null) *outlen = (int)(a.zout - a.zout_start);
                return a.zout_start;
            }
            else
            {
                CLib.CStdlib.free(a.zout_start);
                return null;
            }
        }

        public static int stbi_zlib_decode_buffer(sbyte* obuffer, int olen, sbyte* ibuffer, int ilen)
        {
            stbi__zbuf a;
            a.zbuffer = (byte*)ibuffer;
            a.zbuffer_end = (byte*)ibuffer + ilen;
            if (stbi__do_zlib(&a, obuffer, olen, 0, 1) != 0)
                return (int)(a.zout - a.zout_start);
            else
                return -1;
        }

        public static sbyte* stbi_zlib_decode_noheader_malloc(sbyte* buffer, int len, int* outlen)
        {
            stbi__zbuf a;
            sbyte* p = (sbyte*)CLib.CStdlib.malloc(16384);
            if (p == null) return null;
            a.zbuffer = (byte*)buffer;
            a.zbuffer_end = (byte*)buffer + len;
            if (stbi__do_zlib(&a, p, 16384, 1, 0) != 0)
            {
                if (outlen != null) *outlen = (int)(a.zout - a.zout_start);
                return a.zout_start;
            }
            else
            {
                CLib.CStdlib.free(a.zout_start);
                return null;
            }
        }

        public static int stbi_zlib_decode_noheader_buffer(sbyte* obuffer, int olen, sbyte* ibuffer, int ilen)
        {
            stbi__zbuf a;
            a.zbuffer = (byte*)ibuffer;
            a.zbuffer_end = (byte*)ibuffer + ilen;
            if (stbi__do_zlib(&a, obuffer, olen, 0, 0) != 0)
                return (int)(a.zout - a.zout_start);
            else
                return -1;
        }

#endif //!STBI_NO_ZLIB

        #endregion zlib


        #region png

#if !STBI_NO_PNG

        public struct stbi__pngchunk
        {
            public uint length;
            public uint type;
        }

        static stbi__pngchunk stbi__get_chunk_header(stbi__context* s)
        {
            stbi__pngchunk c;
            c.length = stbi__get32be(s);
            c.type = stbi__get32be(s);
            return c;
        }

        static int stbi__check_png_header(stbi__context* s)
        {
            byte[] png_sig = new byte[8] { 137, 80, 78, 71, 13, 10, 26, 10 };
            int i;
            for (i = 0; i < 8; ++i)
                if (stbi__get8(s) != png_sig[i]) throw new Exception("bad png sig:Not a PNG");
            return 1;
        }

        public struct stbi__png
        {
            public stbi__context* s;
            public byte* idata;
            public byte* expanded;
            public byte* Out;
        }

        public const int STBI__F_none = 0;
        public const int STBI__F_sub = 1;
        public const int STBI__F_up = 2;
        public const int STBI__F_avg = 3;
        public const int STBI__F_paeth = 4;
        // synthetic filters used for first scanline to avoid needing a dummy row of 0s
        public const int STBI__F_avg_first = 5;
        public const int STBI__F_paeth_first = 6;

        static byte[] first_row_filter = new byte[5] 
        {
           STBI__F_none,
           STBI__F_sub,
           STBI__F_none,
           STBI__F_avg_first,
           STBI__F_paeth_first
        };

        static int stbi__paeth(int a, int b, int c)
        {
            int p = a + b - c;
            int pa = Math.Abs(p - a);
            int pb = Math.Abs(p - b);
            int pc = Math.Abs(p - c);
            if (pa <= pb && pa <= pc) return a;
            if (pb <= pc) return b;
            return c;
        }

        static byte[] stbi__depth_scale_table = new byte[9] { 0, 0xff, 0x55, 0, 0x11, 0, 0, 0, 0x01 };

        // create the png data from post-deflated data
        static int stbi__create_png_image_raw(stbi__png* a, byte* raw, uint raw_len, int out_n, uint x, uint y, int depth, int color)
        {
            stbi__context* s = a->s;
            uint i, j, stride = x * (uint)out_n;
            uint img_len, img_width_bytes;
            int k;
            int img_n = s->img_n; // copy it into a local for later

            Debug.Assert(out_n == s->img_n || out_n == s->img_n + 1);
            a->Out = (byte*)CLib.CStdlib.malloc((int)x * (int)y * out_n); // extra bytes to write off the end into
            if (a->Out == null) throw new Exception("outofmem:Out of memory");

            img_width_bytes = ((((uint)img_n * x * (uint)depth) + 7) >> 3);
            img_len = (img_width_bytes + 1) * y;
            if (s->img_x == x && s->img_y == y)
            {
                if (raw_len != img_len) throw new Exception("not enough pixels:Corrupt PNG");
            }
            else
            { // interlaced:
                if (raw_len < img_len) throw new Exception("not enough pixels:Corrupt PNG");
            }

            for (j = 0; j < y; ++j)
            {
                byte* cur = a->Out + stride * j;
                byte* prior = cur - stride;
                int filter = *raw++;
                int filter_bytes = img_n;
                int width = (int)x;
                if (filter > 4)
                    throw new Exception("invalid filter:Corrupt PNG");

                if (depth < 8)
                {
                    Debug.Assert(img_width_bytes <= x);
                    cur += x * out_n - img_width_bytes; // store output to the rightmost img_len bytes, so we can decode in place
                    filter_bytes = 1;
                    width = (int)img_width_bytes;
                }

                // if first row, use special filter that doesn't sample previous row
                if (j == 0) filter = first_row_filter[filter];

                // handle first byte explicitly
                for (k = 0; k < filter_bytes; ++k)
                {
                    switch (filter)
                    {
                        case STBI__F_none: cur[k] = raw[k]; break;
                        case STBI__F_sub: cur[k] = raw[k]; break;
                        case STBI__F_up: cur[k] = (byte)(raw[k] + prior[k]); break;
                        case STBI__F_avg: cur[k] = (byte)(raw[k] + (prior[k] >> 1)); break;
                        case STBI__F_paeth: cur[k] = (byte)(raw[k] + stbi__paeth(0, prior[k], 0)); break;
                        case STBI__F_avg_first: cur[k] = raw[k]; break;
                        case STBI__F_paeth_first: cur[k] = raw[k]; break;
                    }
                }

                if (depth == 8)
                {
                    if (img_n != out_n)
                        cur[img_n] = 255; // first pixel
                    raw += img_n;
                    cur += out_n;
                    prior += out_n;
                }
                else
                {
                    raw += 1;
                    cur += 1;
                    prior += 1;
                }

                // this is a little gross, so that we don't switch per-pixel or per-component
                if (depth < 8 || img_n == out_n)
                {
                    int nk = (width - 1) * img_n;
                    switch (filter)
                    {
                        // "none" filter turns into a memcpy here; make that explicit.
                        case STBI__F_none: CLib.CString.memcpy(cur, raw, (uint)nk); break;
                        case STBI__F_sub: for (k = 0; k < nk; ++k) cur[k] = (byte)(raw[k] + cur[k - filter_bytes]); break;
                        case STBI__F_up: for (k = 0; k < nk; ++k) cur[k] = (byte)(raw[k] + prior[k]); break;
                        case STBI__F_avg: for (k = 0; k < nk; ++k) cur[k] = (byte)(raw[k] + ((prior[k] + cur[k - filter_bytes]) >> 1)); break;
                        case STBI__F_paeth: for (k = 0; k < nk; ++k) cur[k] = (byte)(raw[k] + stbi__paeth(cur[k - filter_bytes], prior[k], prior[k - filter_bytes])); break;
                        case STBI__F_avg_first: for (k = 0; k < nk; ++k) cur[k] = (byte)(raw[k] + (cur[k - filter_bytes] >> 1)); break;
                        case STBI__F_paeth_first: for (k = 0; k < nk; ++k) cur[k] = (byte)(raw[k] + stbi__paeth(cur[k - filter_bytes], 0, 0)); break;
                    }
                    raw += nk;
                }
                else
                {
                    switch (filter)
                    {
                        case STBI__F_none: for (i = x - 1; i >= 1; --i, cur[img_n] = 255, raw += img_n, cur += out_n, prior += out_n) for (k = 0; k < img_n; ++k) cur[k] = raw[k]; break;
                        case STBI__F_sub: for (i = x - 1; i >= 1; --i, cur[img_n] = 255, raw += img_n, cur += out_n, prior += out_n) for (k = 0; k < img_n; ++k) cur[k] = (byte)(raw[k] + cur[k - out_n]); break;
                        case STBI__F_up: for (i = x - 1; i >= 1; --i, cur[img_n] = 255, raw += img_n, cur += out_n, prior += out_n) for (k = 0; k < img_n; ++k) cur[k] = (byte)(raw[k] + prior[k]); break;
                        case STBI__F_avg: for (i = x - 1; i >= 1; --i, cur[img_n] = 255, raw += img_n, cur += out_n, prior += out_n) for (k = 0; k < img_n; ++k) cur[k] = (byte)(raw[k] + ((prior[k] + cur[k - out_n]) >> 1)); break;
                        case STBI__F_paeth: for (i = x - 1; i >= 1; --i, cur[img_n] = 255, raw += img_n, cur += out_n, prior += out_n) for (k = 0; k < img_n; ++k) cur[k] = (byte)(raw[k] + stbi__paeth(cur[k - out_n], prior[k], prior[k - out_n])); break;
                        case STBI__F_avg_first: for (i = x - 1; i >= 1; --i, cur[img_n] = 255, raw += img_n, cur += out_n, prior += out_n) for (k = 0; k < img_n; ++k) cur[k] = (byte)(raw[k] + (cur[k - out_n] >> 1)); break;
                        case STBI__F_paeth_first: for (i = x - 1; i >= 1; --i, cur[img_n] = 255, raw += img_n, cur += out_n, prior += out_n) for (k = 0; k < img_n; ++k) cur[k] = (byte)(raw[k] + stbi__paeth(cur[k - out_n], 0, 0)); break;
                    }
                }
            }

            // we make a separate pass to expand bits to pixels; for performance,
            // this could run two scanlines behind the above code, so it won't
            // intefere with filtering but will still be in the cache.
            if (depth < 8)
            {
                for (j = 0; j < y; ++j)
                {
                    byte* cur = a->Out + stride * j;
                    byte* In = a->Out + stride * j + x * out_n - img_width_bytes;
                    // unpack 1/2/4-bit into a 8-bit buffer. allows us to keep the common 8-bit path optimal at minimal cost for 1/2/4-bit
                    // png guarante byte alignment, if width is not multiple of 8/4/2 we'll decode dummy trailing data that will be skipped in the later loop
                    byte scale = (color == 0) ? (byte)stbi__depth_scale_table[depth] : (byte)1; // scale grayscale values to 0..255 range

                    // note that the final byte might overshoot and write more data than desired.
                    // we can allocate enough data that this never writes out of memory, but it
                    // could also overwrite the next scanline. can it overwrite non-empty data
                    // on the next scanline? yes, consider 1-pixel-wide scanlines with 1-bit-per-pixel.
                    // so we need to explicitly clamp the final ones

                    if (depth == 4)
                    {
                        for (k = (int)x * img_n; k >= 2; k -= 2, ++In)
                        {
                            *cur++ = (byte)(scale * ((*In >> 4)));
                            *cur++ = (byte)(scale * ((*In) & 0x0f));
                        }
                        if (k > 0) *cur++ = (byte)(scale * ((*In >> 4)));
                    }
                    else if (depth == 2)
                    {
                        for (k = (int)x * img_n; k >= 4; k -= 4, ++In)
                        {
                            *cur++ = (byte)(scale * ((*In >> 6)));
                            *cur++ = (byte)(scale * ((*In >> 4) & 0x03));
                            *cur++ = (byte)(scale * ((*In >> 2) & 0x03));
                            *cur++ = (byte)(scale * ((*In) & 0x03));
                        }
                        if (k > 0) *cur++ = (byte)(scale * ((*In >> 6)));
                        if (k > 1) *cur++ = (byte)(scale * ((*In >> 4) & 0x03));
                        if (k > 2) *cur++ = (byte)(scale * ((*In >> 2) & 0x03));
                    }
                    else if (depth == 1)
                    {
                        for (k = (int)x * img_n; k >= 8; k -= 8, ++In)
                        {
                            *cur++ = (byte)(scale * ((*In >> 7)));
                            *cur++ = (byte)(scale * ((*In >> 6) & 0x01));
                            *cur++ = (byte)(scale * ((*In >> 5) & 0x01));
                            *cur++ = (byte)(scale * ((*In >> 4) & 0x01));
                            *cur++ = (byte)(scale * ((*In >> 3) & 0x01));
                            *cur++ = (byte)(scale * ((*In >> 2) & 0x01));
                            *cur++ = (byte)(scale * ((*In >> 1) & 0x01));
                            *cur++ = (byte)(scale * ((*In) & 0x01));
                        }
                        if (k > 0) *cur++ = (byte)(scale * ((*In >> 7)));
                        if (k > 1) *cur++ = (byte)(scale * ((*In >> 6) & 0x01));
                        if (k > 2) *cur++ = (byte)(scale * ((*In >> 5) & 0x01));
                        if (k > 3) *cur++ = (byte)(scale * ((*In >> 4) & 0x01));
                        if (k > 4) *cur++ = (byte)(scale * ((*In >> 3) & 0x01));
                        if (k > 5) *cur++ = (byte)(scale * ((*In >> 2) & 0x01));
                        if (k > 6) *cur++ = (byte)(scale * ((*In >> 1) & 0x01));
                    }
                    if (img_n != out_n)
                    {
                        // insert alpha = 255
                        byte* cur2 = a->Out + stride * j;
                        int ii;
                        if (img_n == 1)
                        {
                            for (ii = (int)x - 1; ii >= 0; --ii)
                            {
                                cur2[ii * 2 + 1] = 255;
                                cur2[ii * 2 + 0] = cur2[ii];
                            }
                        }
                        else
                        {
                            Debug.Assert(img_n == 3);
                            for (ii = (int)x - 1; ii >= 0; --ii)
                            {
                                cur2[ii * 4 + 3] = 255;
                                cur2[ii * 4 + 2] = cur2[ii * 3 + 2];
                                cur2[ii * 4 + 1] = cur2[ii * 3 + 1];
                                cur2[ii * 4 + 0] = cur2[ii * 3 + 0];
                            }
                        }
                    }
                }
            }

            return 1;
        }

        static int stbi__create_png_image(stbi__png* a, byte* image_data, uint image_data_len, int out_n, int depth, int color, int interlaced)
        {
            byte* final;
            int p;
            if (interlaced == 0)
                return stbi__create_png_image_raw(a, image_data, image_data_len, out_n, a->s->img_x, a->s->img_y, depth, color);

            // de-interlacing
            final = (byte*)CLib.CStdlib.malloc((int)(a->s->img_x * a->s->img_y * out_n));
            for (p = 0; p < 7; ++p)
            {
                int[] xorig = new int[] { 0, 4, 0, 2, 0, 1, 0 };
                int[] yorig = new int[] { 0, 0, 4, 0, 2, 0, 1 };
                int[] xspc = new int[] { 8, 8, 4, 4, 2, 2, 1 };
                int[] yspc = new int[] { 8, 8, 8, 4, 4, 2, 2 };
                int i, j, x, y;
                // pass1_x[4] = 0, pass1_x[5] = 1, pass1_x[12] = 1
                x = (int)(a->s->img_x - xorig[p] + xspc[p] - 1) / xspc[p];
                y = (int)(a->s->img_y - yorig[p] + yspc[p] - 1) / yspc[p];
                if (x != 0 && y != 0)
                {
                    uint img_len = (uint)(((((a->s->img_n * x * depth) + 7) >> 3) + 1) * y);
                    if (stbi__create_png_image_raw(a, image_data, image_data_len, out_n, (uint)x, (uint)y, depth, color) == 0)
                    {
                        CLib.CStdlib.free(final);
                        return 0;
                    }
                    for (j = 0; j < y; ++j)
                    {
                        for (i = 0; i < x; ++i)
                        {
                            int out_y = j * yspc[p] + yorig[p];
                            int out_x = i * xspc[p] + xorig[p];
                            CLib.CString.memcpy(final + out_y * a->s->img_x * out_n + out_x * out_n, a->Out + (j * x + i) * out_n, (uint)out_n);
                        }
                    }
                    CLib.CStdlib.free(a->Out);
                    image_data += img_len;
                    image_data_len -= img_len;
                }
            }
            a->Out = final;

            return 1;
        }

        static int stbi__compute_transparency(stbi__png* z, byte* tc, int out_n)
        {
            stbi__context* s = z->s;
            uint i, pixel_count = s->img_x * s->img_y;
            byte* p = z->Out;

            // compute color-based transparency, assuming we've
            // already got 255 as the alpha value in the output
            Debug.Assert(out_n == 2 || out_n == 4);

            if (out_n == 2)
            {
                for (i = 0; i < pixel_count; ++i)
                {
                    p[1] = (byte)(p[0] == tc[0] ? 0 : 255);
                    p += 2;
                }
            }
            else
            {
                for (i = 0; i < pixel_count; ++i)
                {
                    if (p[0] == tc[0] && p[1] == tc[1] && p[2] == tc[2])
                        p[3] = 0;
                    p += 4;
                }
            }
            return 1;
        }

        static int stbi__expand_png_palette(stbi__png* a, byte* palette, int len, int pal_img_n)
        {
            uint i, pixel_count = a->s->img_x * a->s->img_y;
            byte* p;
            byte* temp_out;
            byte* orig = a->Out;

            p = (byte*)CLib.CStdlib.malloc((int)(pixel_count * pal_img_n));
            if (p == null) throw new Exception("outofmem:Out of memory");

            // between here and free(out) below, exitting would leak
            temp_out = p;

            if (pal_img_n == 3)
            {
                for (i = 0; i < pixel_count; ++i)
                {
                    int n = orig[i] * 4;
                    p[0] = palette[n];
                    p[1] = palette[n + 1];
                    p[2] = palette[n + 2];
                    p += 3;
                }
            }
            else
            {
                for (i = 0; i < pixel_count; ++i)
                {
                    int n = orig[i] * 4;
                    p[0] = palette[n];
                    p[1] = palette[n + 1];
                    p[2] = palette[n + 2];
                    p[3] = palette[n + 3];
                    p += 4;
                }
            }
            CLib.CStdlib.free(a->Out);
            a->Out = temp_out;

            //STBI_NOTUSED(len);

            return 1;
        }

        static int stbi__unpremultiply_on_load = 0;
        static int stbi__de_iphone_flag = 0;

        public static void stbi_set_unpremultiply_on_load(int flag_true_if_should_unpremultiply)
        {
            stbi__unpremultiply_on_load = flag_true_if_should_unpremultiply;
        }

        public static void stbi_convert_iphone_png_to_rgb(int flag_true_if_should_convert)
        {
            stbi__de_iphone_flag = flag_true_if_should_convert;
        }

        static void stbi__de_iphone(stbi__png* z)
        {
            stbi__context* s = z->s;
            uint i, pixel_count = s->img_x * s->img_y;
            byte* p = z->Out;

            if (s->img_out_n == 3)
            {  // convert bgr to rgb
                for (i = 0; i < pixel_count; ++i)
                {
                    byte t = p[0];
                    p[0] = p[2];
                    p[2] = t;
                    p += 3;
                }
            }
            else
            {
                Debug.Assert(s->img_out_n == 4);
                if (stbi__unpremultiply_on_load != 0)
                {
                    // convert bgr to rgb and unpremultiply
                    for (i = 0; i < pixel_count; ++i)
                    {
                        byte a = p[3];
                        byte t = p[0];
                        if (a != 0)
                        {
                            p[0] = (byte)(p[2] * 255 / a);
                            p[1] = (byte)(p[1] * 255 / a);
                            p[2] = (byte)(t * 255 / a);
                        }
                        else
                        {
                            p[0] = p[2];
                            p[2] = t;
                        }
                        p += 4;
                    }
                }
                else
                {
                    // convert bgr to rgb
                    for (i = 0; i < pixel_count; ++i)
                    {
                        byte t = p[0];
                        p[0] = p[2];
                        p[2] = t;
                        p += 4;
                    }
                }
            }
        }

        static uint STBI__PNG_TYPE(int a, int b, int c, int d)
        {
            return (uint)(((a) << 24) + ((b) << 16) + ((c) << 8) + (d));
        }

        public static int stbi__parse_png_file(stbi__png* z, int scan, int req_comp)
        {
            byte* palette = stackalloc byte[1024]; byte pal_img_n = 0;
            byte has_trans = 0; byte* tc = stackalloc byte[3];
            uint ioff = 0, idata_limit = 0, i, pal_len = 0;
            int first = 1, k, interlace = 0, color = 0, depth = 0, is_iphone = 0;
            stbi__context* s = z->s;

            z->expanded = null;
            z->idata = null;
            z->Out = null;

            if (stbi__check_png_header(s) == 0) return 0;

            if (scan == STBI__SCAN_type) return 1;

            for (; ; )
            {
                stbi__pngchunk c = stbi__get_chunk_header(s);
                switch (c.type)
                {
                    case 0x43674249:    //STBI__PNG_TYPE('C', 'g', 'B', 'I'):
                        is_iphone = 1;
                        stbi__skip(s, (int)c.length);
                        break;
                    case 0x49484452:    //STBI__PNG_TYPE('I', 'H', 'D', 'R'):
                        {
                            int comp, filter;
                            if (first == 0) throw new Exception("multiple IHDR:Corrupt PNG");
                            first = 0;
                            if (c.length != 13) throw new Exception("bad IHDR len:Corrupt PNG");
                            s->img_x = stbi__get32be(s); if (s->img_x > (1 << 24)) throw new Exception("too large:Very large image (corrupt?)");
                            s->img_y = stbi__get32be(s); if (s->img_y > (1 << 24)) throw new Exception("too large:Very large image (corrupt?)");
                            depth = stbi__get8(s); if (depth != 1 && depth != 2 && depth != 4 && depth != 8) throw new Exception("1/2/4/8-bit only:PNG not supported: 1/2/4/8-bit only");
                            color = stbi__get8(s); if (color > 6) throw new Exception("bad ctype:Corrupt PNG");
                            if (color == 3) pal_img_n = 3; else if ((color & 1) != 0) throw new Exception("bad ctype:Corrupt PNG");
                            comp = stbi__get8(s); if (comp != 0) throw new Exception("bad comp method:Corrupt PNG");
                            filter = stbi__get8(s); if (filter != 0) throw new Exception("bad filter method:Corrupt PNG");
                            interlace = stbi__get8(s); if (interlace > 1) throw new Exception("bad interlace method:Corrupt PNG");
                            if (s->img_x == 0 || s->img_y == 0) throw new Exception("0-pixel image:Corrupt PNG");
                            if (pal_img_n == 0)
                            {
                                s->img_n = ((color & 2) != 0 ? 3 : 1) + ((color & 4) != 0 ? 1 : 0);
                                if ((1 << 30) / s->img_x / s->img_n < s->img_y) throw new Exception("too large:Image too large to decode");
                                if (scan == STBI__SCAN_header) return 1;
                            }
                            else
                            {
                                // if paletted, then pal_n is our final components, and
                                // img_n is # components to decompress/filter.
                                s->img_n = 1;
                                if ((1 << 30) / s->img_x / 4 < s->img_y) throw new Exception("too large:Corrupt PNG");
                                // if SCAN_header, have to scan to see if we have a tRNS
                            }
                            break;
                        }

                    case 0x504c5445:    //STBI__PNG_TYPE('P', 'L', 'T', 'E'):
                        {
                            if (first != 0) throw new Exception("first not IHDR:Corrupt PNG");
                            if (c.length > 256 * 3) throw new Exception("invalid PLTE:Corrupt PNG");
                            pal_len = c.length / 3;
                            if (pal_len * 3 != c.length) throw new Exception("invalid PLTE:Corrupt PNG");
                            for (i = 0; i < pal_len; ++i)
                            {
                                palette[i * 4 + 0] = stbi__get8(s);
                                palette[i * 4 + 1] = stbi__get8(s);
                                palette[i * 4 + 2] = stbi__get8(s);
                                palette[i * 4 + 3] = 255;
                            }
                            break;
                        }

                    case 0x74524e53:     //STBI__PNG_TYPE('t', 'R', 'N', 'S'):
                        {
                            if (first != 0) throw new Exception("first not IHDR:Corrupt PNG");
                            if (z->idata != null) throw new Exception("tRNS after IDAT:Corrupt PNG");
                            if (pal_img_n != 0)
                            {
                                if (scan == STBI__SCAN_header) { s->img_n = 4; return 1; }
                                if (pal_len == 0) throw new Exception("tRNS before PLTE:Corrupt PNG");
                                if (c.length > pal_len) throw new Exception("bad tRNS len:Corrupt PNG");
                                pal_img_n = 4;
                                for (i = 0; i < c.length; ++i)
                                    palette[i * 4 + 3] = stbi__get8(s);
                            }
                            else
                            {
                                if ((s->img_n & 1) == 0) throw new Exception("tRNS with alpha:Corrupt PNG");
                                if (c.length != (uint)s->img_n * 2) throw new Exception("bad tRNS len:Corrupt PNG");
                                has_trans = 1;
                                for (k = 0; k < s->img_n; ++k)
                                    tc[k] = (byte)((stbi__get16be(s) & 255) * stbi__depth_scale_table[depth]); // non 8-bit images will be larger
                            }
                            break;
                        }

                    case 0x49444154:    //STBI__PNG_TYPE('I', 'D', 'A', 'T'):
                        {
                            if (first != 0) throw new Exception("first not IHDR:Corrupt PNG");
                            if (pal_img_n != 0 && pal_len == 0) throw new Exception("no PLTE:Corrupt PNG");
                            if (scan == STBI__SCAN_header) { s->img_n = pal_img_n; return 1; }
                            if ((int)(ioff + c.length) < (int)ioff) return 0;
                            if (ioff + c.length > idata_limit)
                            {
                                byte* p;
                                if (idata_limit == 0) idata_limit = c.length > 4096 ? c.length : 4096;
                                while (ioff + c.length > idata_limit)
                                    idata_limit *= 2;
                                p = (byte*)CLib.CStdlib.realloc(z->idata, (int)idata_limit); if (p == null) throw new Exception("outofmem:Out of memory");
                                z->idata = p;
                            }
                            if (stbi__getn(s, z->idata + ioff, (int)c.length) == 0) throw new Exception("outofdata:Corrupt PNG");
                            ioff += c.length;
                            break;
                        }

                    case 0x49454e44:    //STBI__PNG_TYPE('I', 'E', 'N', 'D'):
                        {
                            uint raw_len, bpl;
                            if (first != 0) throw new Exception("first not IHDR:Corrupt PNG");
                            if (scan != STBI__SCAN_load) return 1;
                            if (z->idata == null) throw new Exception("no IDAT:Corrupt PNG");
                            // initial guess for decoded data size to avoid unnecessary reallocs
                            bpl = (s->img_x * (uint)depth + 7) / 8; // bytes per line, per component
                            raw_len = (uint)(bpl * s->img_y * s->img_n /* pixels */ + s->img_y /* filter mode per row */);
                            z->expanded = (byte*)stbi_zlib_decode_malloc_guesssize_headerflag((sbyte*)z->idata, (int)ioff, (int)raw_len, (int*)&raw_len, is_iphone);
                            if (z->expanded == null) return 0; // zlib should set error
                            CLib.CStdlib.free(z->idata); z->idata = null;
                            if ((req_comp == s->img_n + 1 && req_comp != 3 && pal_img_n == 0) || has_trans != 0)
                                s->img_out_n = s->img_n + 1;
                            else
                                s->img_out_n = s->img_n;
                            if (stbi__create_png_image(z, z->expanded, raw_len, s->img_out_n, depth, color, interlace) == 0) return 0;
                            if (has_trans != 0)
                                if (stbi__compute_transparency(z, tc, s->img_out_n) == 0) return 0;
                            if (is_iphone != 0 && stbi__de_iphone_flag != 0 && s->img_out_n > 2)
                                stbi__de_iphone(z);
                            if (pal_img_n != 0)
                            {
                                // pal_img_n == 3 or 4
                                s->img_n = pal_img_n; // record the actual colors we had
                                s->img_out_n = pal_img_n;
                                if (req_comp >= 3) s->img_out_n = req_comp;
                                if (stbi__expand_png_palette(z, palette, (int)pal_len, (int)s->img_out_n) == 0)
                                    return 0;
                            }
                            CLib.CStdlib.free(z->expanded); z->expanded = null;
                            return 1;
                        }

                    default:
                        // if critical, fail
                        if (first != 0) throw new Exception("first not IHDR:Corrupt PNG");
                        if ((c.type & (1 << 29)) == 0)
                        {
                            // not threadsafe
                            char[] invalid_chunk = "XXXX PNG chunk not known".ToCharArray();
                            invalid_chunk[0] = (char)(c.type >> 24);
                            invalid_chunk[1] = (char)(c.type >> 16);
                            invalid_chunk[2] = (char)(c.type >> 8);
                            invalid_chunk[3] = (char)(c.type >> 0);
                            throw new Exception(new string(invalid_chunk) + ":PNG not supported: unknown PNG chunk type");
                        }
                        stbi__skip(s, (int)c.length);
                        break;
                }
                // end of PNG chunk, read and skip CRC
                stbi__get32be(s);
            }
        }

        static byte* stbi__do_png(stbi__png* p, int* x, int* y, int* n, int req_comp)
        {
            byte* result = null;
            if (req_comp < 0 || req_comp > 4) throw new Exception("bad req_comp:Internal error");
            if (stbi__parse_png_file(p, STBI__SCAN_load, req_comp) != 0)
            {
                result = p->Out;
                p->Out = null;
                if (req_comp != 0 && req_comp != p->s->img_out_n)
                {
                    result = stbi__convert_format(result, p->s->img_out_n, req_comp, p->s->img_x, p->s->img_y);
                    p->s->img_out_n = req_comp;
                    if (result == null) return result;
                }
                *x = (int)p->s->img_x;
                *y = (int)p->s->img_y;
                if (n != null) *n = p->s->img_out_n;
            }
            CLib.CStdlib.free(p->Out); p->Out = null;
            CLib.CStdlib.free(p->expanded); p->expanded = null;
            CLib.CStdlib.free(p->idata); p->idata = null;

            return result;
        }

        static byte* stbi__png_load(stbi__context* s, int* x, int* y, int* comp, int req_comp)
        {
            stbi__png p;
            p.s = s;
            return stbi__do_png(&p, x, y, comp, req_comp);
        }

        static int stbi__png_test(stbi__context* s)
        {
            int r;
            r = stbi__check_png_header(s);
            stbi__rewind(s);
            return r;
        }

        static int stbi__png_info_raw(stbi__png* p, int* x, int* y, int* comp)
        {
            if (stbi__parse_png_file(p, STBI__SCAN_header, 0) == 0)
            {
                stbi__rewind(p->s);
                return 0;
            }
            if (x != null) *x = (int)p->s->img_x;
            if (y != null) *y = (int)p->s->img_y;
            if (comp != null) *comp = p->s->img_n;
            return 1;
        }

        static int stbi__png_info(stbi__context* s, int* x, int* y, int* comp)
        {
            stbi__png p;
            p.s = s;
            return stbi__png_info_raw(&p, x, y, comp);
        }

#endif //!STBI_NO_PNG

        #endregion png

    }
}