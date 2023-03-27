#include "olcConsoleGameEngine.h"
#include <fstream>
#include <strstream>
#include <algorithm>

#define PI 3.1415925f

struct vector3
{
	float x, y, z;
	vector3(float X = 0.0f, float Y = 0.0f, float Z = 0.0f) : x(X), y(Y), z(Z) {};
	
	float magnitude() { return sqrtf(x * x + y * y + z * z); }
	vector3 normalized()
	{
		return { x / magnitude(), y / magnitude(), z / magnitude() };
	}

	vector3 operator+(const vector3& v2) const { return { x + v2.x, y + v2.y, z + v2.z, }; }
	vector3 operator-(const vector3& v2) const { return { x - v2.x, y - v2.y, z - v2.z, }; }
	float	operator*(const vector3& v2) const { return x * v2.x + y * v2.y + z * v2.z; }
	vector3 operator^(const vector3& v2) const {
		vector3 v;
		v.x = y * v2.z - z * v2.y;
		v.y = z * v2.x - x * v2.z;
		v.z = x * v2.y - y * v2.x;
		return v;
	}

};
struct Quaternion : vector3
{
	float w;
	Quaternion(float X = 0.0f, float Y = 0.0f, float Z = 0.0f, float W = 0.0f) : vector3(X, Y, Z), w(W) {};
	Quaternion(vector3 axis = {0, 1, 0}, float angle = 0.0f)
	{
		x = cosf(angle / 2);
		y = axis.x * sinf(angle / 2);
		z = axis.y * sinf(angle / 2);
		w = axis.z * sinf(angle / 2);
	};

	vector3 rotate(const vector3& p) const
	{
		vector3 r;
		r.x = ((-y * p.x - z * p.y - w * p.z) * (-y) + (x * p.x + z * p.z - w * p.y) * (x) + (x * p.y - y * p.z + w * p.x) * (-w) - (x * p.z + y * p.y - z * p.x) * (-z));
		r.y = ((-y * p.x - z * p.y - w * p.z) * (-z) - (x * p.x + z * p.z - w * p.y) * (-w) + (x * p.y - y * p.z + w * p.x) * (x) + (x * p.z + y * p.y - z * p.x) * (-y));
		r.z = ((-y * p.x - z * p.y - w * p.z) * (-w) + (x * p.x + z * p.z - w * p.y) * (-z) - (x * p.y - y * p.z + w * p.x) * (-y) + (x * p.z + y * p.y - z * p.x) * (x));
		return r;
	}
};

struct triangle
{
	vector3 p[3];

	wchar_t sym;
	short col;
};

struct mesh
{
	std::vector<triangle> tris;

	bool LoadObjectFromFile(std::string sFilename)
	{
		std::ifstream f(sFilename);
		if (!f.is_open())
			return false;

		// Local cache of verts
		std::vector<vector3> verts;

		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			std::strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				vector3 v;
				s >> junk >> v.x >> v.y >> v.z;
				verts.push_back(v);
			}

			if (line[0] == 'f')
			{
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];
				tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
			}
		}

		return true;
	}

};

struct mat4x4
{
	float m[4][4] = { 0 };
};

class olcEngine3D : public olcConsoleGameEngine
{
public:
	olcEngine3D()
	{
		m_sAppName = L"3D Demo";
	}


private:
	mesh meshCube;
	mat4x4 matProj;

	vector3 Camera;


	float theta;

	void MultiplyMatrixVector(vector3& i, vector3& o, mat4x4& m)
	{
		o.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
		o.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
		o.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
		float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];

		if (w != 0.0f)
		{
			o.x /= w; o.y /= w; o.z /= w;
		}
	}

	CHAR_INFO GetColour(float lum)
	{
		short bg_col, fg_col;
		wchar_t sym;
		int pixel_bw = (int)(13.0f * lum);
		switch (pixel_bw)
		{
		case 0: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID; break;

		case 1: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_QUARTER; break;
		case 2: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_HALF; break;
		case 3: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 4: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_SOLID; break;

		case 5: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_QUARTER; break;
		case 6: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_HALF; break;
		case 7: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 8: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_SOLID; break;

		case 9:  bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_QUARTER; break;
		case 10: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_HALF; break;
		case 11: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_THREEQUARTERS; break;
		case 12: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_SOLID; break;
		default:
			bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID;
		}

		CHAR_INFO c;
		c.Attributes = bg_col | fg_col;
		c.Char.UnicodeChar = sym;
		return c;
	}

public:
	bool OnUserCreate() override
	{

		meshCube.LoadObjectFromFile("meshes/bunny.obj");

		float Near = 0.1f;
		float Far = 1000.0f;
		float Fov = 90.0f;
		float AspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
		float FovRad = 1.0f / tanf(Fov * 0.5f / 180.0f * 3.14159f);

		matProj.m[0][0] = AspectRatio * FovRad;
		matProj.m[1][1] = FovRad;
		matProj.m[2][2] = Far / (Far - Near);
		matProj.m[3][2] = (-Far * Near) / (Far - Near);
		matProj.m[2][3] = 1.0f;
		matProj.m[3][3] = 0.0f;

		return true;
	}

	bool OnUserUpdate(float deltaTime) override
	{
		// Clear Screen
		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

		// Set up quaternions
		theta += 1.0f * deltaTime;

		Quaternion rotX( { 1, 0, 0 }, 0 );
		Quaternion rotY( { 0, 1, 0 }, theta );
		Quaternion rotZ( { 0, 0, 1 }, 0 );

		// Store triagles for rastering later
		std::vector<triangle> trianglesToRaster;

		// Draw Triangles
		for (triangle& tri : meshCube.tris)
		{
			triangle triProjected, triTranslated, triRotatedX, triRotatedXY, triRotatedXYZ;

			// Rotate in X-Axis
			triRotatedX.p[0] = rotX.rotate(tri.p[0]);
			triRotatedX.p[1] = rotX.rotate(tri.p[1]);
			triRotatedX.p[2] = rotX.rotate(tri.p[2]);

			// Rotate in Y-Axis
			triRotatedXY.p[0] = rotY.rotate(triRotatedX.p[0]);
			triRotatedXY.p[1] = rotY.rotate(triRotatedX.p[1]);
			triRotatedXY.p[2] = rotY.rotate(triRotatedX.p[2]);

			// Rotate in Z-Axis
			triRotatedXYZ.p[0] = rotZ.rotate(triRotatedXY.p[0]);
			triRotatedXYZ.p[1] = rotZ.rotate(triRotatedXY.p[1]);
			triRotatedXYZ.p[2] = rotZ.rotate(triRotatedXY.p[2]);


			// Offset into the screen
			triTranslated = triRotatedXYZ;
			triTranslated.p[0].z = triRotatedXYZ.p[0].z + 4.0f;
			triTranslated.p[1].z = triRotatedXYZ.p[1].z + 4.0f;
			triTranslated.p[2].z = triRotatedXYZ.p[2].z + 4.0f;

			// Use Cross-Product to get surface normal
			vector3 normal, line1, line2;
			
			line1 = triTranslated.p[1] - triTranslated.p[0];

			line2 = triTranslated.p[2] - triTranslated.p[0];
			
			normal = line1 ^ line2;

			// It's normally normal to normalise the normal
			normal = normal.normalized();

			vector3 cameraRay = triTranslated.p[0] - Camera;
			//if (normal.z < 0)
			if (normal * cameraRay < 0.0f)
			{
				// Illumination
				vector3 light_direction = { 0.0f, 0.0f, -1.0f };
				light_direction = light_direction.normalized();

				// How similar is normal to light direction
				float dp = normal * light_direction;

				// Choose console colours as required (much easier with RGB)
				CHAR_INFO c = GetColour(dp);
				triTranslated.col = c.Attributes;
				triTranslated.sym = c.Char.UnicodeChar;

				// Project triangles from 3D --> 2D
				MultiplyMatrixVector(triTranslated.p[0], triProjected.p[0], matProj);
				MultiplyMatrixVector(triTranslated.p[1], triProjected.p[1], matProj);
				MultiplyMatrixVector(triTranslated.p[2], triProjected.p[2], matProj);
				triProjected.col = triTranslated.col;
				triProjected.sym = triTranslated.sym;

				// Scale into view
				triProjected.p[0].x += 1.0f; triProjected.p[0].y += 1.0f;
				triProjected.p[1].x += 1.0f; triProjected.p[1].y += 1.0f;
				triProjected.p[2].x += 1.0f; triProjected.p[2].y += 1.0f;
				triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
				triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
				triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
				triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
				triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
				triProjected.p[2].y *= 0.5f * (float)ScreenHeight();

				// Store triangle for sorting
				trianglesToRaster.push_back(triProjected);
			}

		}

		// Sort triangles from back to front
		std::sort(trianglesToRaster.begin(), trianglesToRaster.end(), [](triangle& t1, triangle& t2)
			{
				float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
				float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
				return z1 > z2;
			});

		for (auto& triProjected : trianglesToRaster)
		{
			// Rasterize triangle
			FillTriangle(triProjected.p[0].x, triProjected.p[0].y,
				triProjected.p[1].x, triProjected.p[1].y,
				triProjected.p[2].x, triProjected.p[2].y,
				triProjected.sym, triProjected.col);

			/*DrawTriangle(triProjected.p[0].x, triProjected.p[0].y,
			triProjected.p[1].x, triProjected.p[1].y,
			triProjected.p[2].x, triProjected.p[2].y,
			PIXEL_SOLID, FG_BLACK);*/
		}


		return true;
	}

};




int main()
{
	olcEngine3D demo;
	if (demo.ConstructConsole(256, 240, 4, 4))
		demo.Start();
	return 0;
}