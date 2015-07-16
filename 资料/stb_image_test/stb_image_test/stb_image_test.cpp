// stb_image_test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


int _tmain(int argc, _TCHAR* argv[])
{
	int x; int y; int comp;
	stbi_uc* buffer = stbi_load("F:\\p.png", &x, &y, &comp, 0);

	return 0;
}

