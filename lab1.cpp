#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

static std::default_random_engine engine (10) ; // random seed = 10
static std::uniform_real_distribution<double> uniform (0,1) ;

double sqr(double x) { return x*x;};

void boxMuller(double stdev, double &x, double &y) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(-2*log(r1))*cos(2.*M_PI*r2)*stdev;
	y = sqrt(-2*log(r1))*sin(2.*M_PI*r2)*stdev;
}
 
class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3]; 
};
 
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector&a, const Vector&b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Vector random_lightDir(Vector& N) {

    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
    double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
    double z = sqrt(r2); 

    Vector T1, T2;
    if (fabs(N[2]) <= std::min(fabs(N[0]), fabs(N[1]))) {
        T1 = Vector(-N[1], N[0], 0);
    } else if (fabs(N[1]) <= std::min(fabs(N[0]), fabs(N[2]))) {
        T1 = Vector(0, -N[2], N[1]);
    } else {
        T1 = Vector(-N[2], 0, N[0]);
    }

    T1.normalize();
    T2 = cross(N, T1);
    return x*T1 + y*T2 + z*N;
}


class Ray {
public:
    Ray(const Vector& O, const Vector &u) : O(O), u(u) {};
    // O is the origin, u is the directional vector
    Vector O,u;

};

class Geometry {
public:
    Geometry(const Vector& alb, bool mirror, bool transparent) : alb(alb),transparent(transparent), mirror(mirror) {};
    Vector alb;
    bool mirror, transparent;
    virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t) = 0;

};

class Sphere : public Geometry {
    public:
    Sphere(const Vector& C, double R, const Vector& alb, bool mirror = false, bool transparent = false)
    : Geometry(alb, mirror, transparent), C(C), R(R) {}
    Vector C; // center
    double R; // radius
    
    virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t) override { // P : point of intersection, N : normal vector (reflected ray)
        // determines whether the ray intersects with the sphere and updates P and N with relevant info
        double delta = sqr(dot(r.u, r.O - C)) - ((r.O - C).norm2() - R * R); // definition in lecture slides
        
        if (delta < 0) {
            return false;
        } else {
        double x = dot(r.u, C - r.O);
        double t1 = x - sqrt(delta);
        double t2 = x + sqrt(delta);
        if (t2 < 0) {
            return false;
        } // sphere is behind the ray
        if (t1 >= 0) { // if t1 negative then ray starts from in the middle of the sphere (direction negative)
            t = t1;
        } else {
            t = t2;
        }
        P = r.O + t * r.u; // extend ray from origin to surface of sphere (multiply unit vector by t)
        N = P - C; // directional vector from point of contact to center of sphere, gives normal
        N.normalize();
        return true;
        }

        return (delta >=0);
    }
};

class BoundingBox { // this is supposed to speed it up
    public:
    BoundingBox(const Vector& min = Vector(0,0,0), const Vector & max = Vector(0,0,0)) : min(min), max(max) {};

    bool intersect(const Ray& r) {
        // Checks if theres an intersection with a ray and the bounding box

        double tx1 = (min[0] - r.O[0])/r.u[0];
        double tx2 = (max[0] - r.O[0])/r.u[0];
        double txMin = std::min(tx1, tx2);
        double txMax = std::max(tx1, tx2);

        double ty1 = (min[1] - r.O[1])/r.u[1];
        double ty2 = (max[1] - r.O[1])/r.u[1];
        double tyMin = std::min(ty1, ty2);
        double tyMax = std::max(ty1, ty2);

        double tz1 = (min[2] - r.O[2])/r.u[2];
        double tz2 = (max[2] - r.O[2])/r.u[2];
        double tzMin = std::min(tz1, tz2);
        double tzMax = std::max(tz1, tz2);

        //Intersects if
        if (std::min(txMax, std::min(tyMax, tzMax)) > std::max(txMin, std::max(tyMin, tzMin))) {
            return true;
        }
        return false;
    }
    Vector min, max;
};

class TriangleIndices {
    public:
        TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
        };
        int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
        int uvi, uvj, uvk;  // indices within the uv coordinates array
        int ni, nj, nk;  // indices within the normals array
        int group;       // face group
};
     
class TriangleMesh : public Geometry {
    public:
        ~TriangleMesh() {};
        TriangleMesh() : Geometry(alb, mirror, transparent), smooth(smooth) {};
        bool smooth = true;

        void readOBJ(const char* obj) {
     
            char matfile[255];
            char grp[255];
     
            FILE* f;
            f = fopen(obj, "r");
            int curGroup = -1;
            while (!feof(f)) {
                char line[255];
                if (!fgets(line, 255, f)) break;
     
                std::string linetrim(line);
                linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
                strcpy(line, linetrim.c_str());
     
                if (line[0] == 'u' && line[1] == 's') {
                    sscanf(line, "usemtl %[^\n]\n", grp);
                    curGroup++;
                }
     
                if (line[0] == 'v' && line[1] == ' ') {
                    Vector vec;
     
                    Vector col;
                    if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                        col[0] = std::min(1., std::max(0., col[0]));
                        col[1] = std::min(1., std::max(0., col[1]));
                        col[2] = std::min(1., std::max(0., col[2]));
     
                        vertices.push_back(vec);
                        vertexcolors.push_back(col);
     
                    } else {
                        sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                        vertices.push_back(vec);
                    }
                }
                if (line[0] == 'v' && line[1] == 'n') {
                    Vector vec;
                    sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    normals.push_back(vec);
                }
                if (line[0] == 'v' && line[1] == 't') {
                    Vector vec;
                    sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                    uvs.push_back(vec);
                }
                if (line[0] == 'f') {
                    TriangleIndices t;
                    int i0, i1, i2, i3;
                    int j0, j1, j2, j3;
                    int k0, k1, k2, k3;
                    int nn;
                    t.group = curGroup;
     
                    char* consumedline = line + 1;
                    int offset;
     
                    nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                    if (nn == 9) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                        if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                        if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                        if (nn == 6) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                            if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                            if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                            if (nn == 3) {
                                if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                                if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                                if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                                indices.push_back(t);
                            } else {
                                nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                                if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                                if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                                if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                                if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                                if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                                if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                                indices.push_back(t);
                            }
                        }
                    }
     
                    consumedline = consumedline + offset;
     
                    while (true) {
                        if (consumedline[0] == '\n') break;
                        if (consumedline[0] == '\0') break;
                        nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                        TriangleIndices t2;
                        t2.group = curGroup;
                        if (nn == 3) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                            if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                            if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                            indices.push_back(t2);
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            k2 = k3;
                        } else {
                            nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                                if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                                if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                j2 = j3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                                if (nn == 2) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                    if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                    if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    k2 = k3;
                                    indices.push_back(t2);
                                } else {
                                    nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                    if (nn == 1) {
                                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                        consumedline = consumedline + offset;
                                        i2 = i3;
                                        indices.push_back(t2);
                                    } else {
                                        consumedline = consumedline + 1;
                                    }
                                }
                            }
                        }
                    }
     
                }
     
            }
            fclose(f);
     
        }

        void scale(double scale, const Vector& transpose) {
            for (int i = 0; i < vertices.size(); i++) {
                vertices[i][0] = vertices[i][0]*scale + transpose[0];
		        vertices[i][1] = vertices[i][1]*scale + transpose[1];
		        vertices[i][2] = vertices[i][2]*scale + transpose[2];
            }
        }

        virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t) override {

            bool intersection = false;

            if (!this->bounding_box.intersect(r)) {
                return false;
            }

            for (int i = 0; i < indices.size(); i++ ) {

                // Getting vertex coordinates
                Vector& A = vertices[indices[i].vtxi];
                Vector& B = vertices[indices[i].vtxj];
                Vector& C = vertices[indices[i].vtxk];

                // Solving intersection equation
                Vector e1 = B - A;
                Vector e2 = C - A;
                Vector z(0,0,0);
                Vector M = cross(e1, e2);
                Vector AO = A - r.O;

                double det = dot(r.u, M);
                if (std::abs(det) < 1e-6) continue;

                double beta = dot(e2, cross(AO, r.u)) / dot(r.u, M);
                double gamma = -dot(e1, cross(AO, r.u)) / dot(r.u, M);
                double alpha = 1 - beta - gamma;
                double localt = dot(AO, M) / dot(r.u, M);

                if (localt >= 0 && localt < t  && beta>=0 && beta <= 1 && alpha>=0 && alpha <= 1 && gamma>=0 && gamma <= 1) {
                    P = A + beta*e1 + gamma*e2;
					t = localt;
                    if (smooth) {
                        // Getting average of normals
					    N = alpha*normals[indices[i].ni] + beta*normals[indices[i].nj] + gamma*normals[indices[i].nk];
                    } else {
                        N = M;
                    }
					N.normalize();
                    intersection = true;
                }
            }
            return intersection;
        }

        void compute_bbox() {
            bounding_box.min = Vector(1E9, 1E9, 1E9); // Progressively try to improve these bounds (min and max)
            bounding_box.max = Vector(-1E9,-1E9, -1E9);

            for (int i = 0; i < vertices.size(); i++) {
                for (int j = 0; j < 3; j++) {
                    bounding_box.min[j] = std::min(bounding_box.min[j], vertices[i][j]);
                    bounding_box.max[j] = std::max(bounding_box.max[j], vertices[i][j]);
                }
            }
        }

        void no_bbox() {
            bounding_box.min = Vector(1E9, 1E9, 1E9);
            bounding_box.max = Vector(-1E9,-1E9, -1E9);
        }

        std::vector<TriangleIndices> indices;
        std::vector<Vector> vertices;
        std::vector<Vector> normals;
        std::vector<Vector> uvs;
        std::vector<Vector> vertexcolors;
        BoundingBox bounding_box;
};

 
class Scene {
    public:
        Scene() {};

        std::vector<Geometry*> objects;
        Vector L; // light source point coordinates
        double I;
        bool direct_lighting;
        int ray_count;
        bool antialiasing;

        bool intersect(const Ray& r, Vector& P, Vector& N, double& t, Geometry*& hit_obj) {
            t = 1E30;
            bool intersected = false;
            Vector localP, localN;
            double localt;

            for (int i = 0; i < objects.size(); i++) {
                if (objects[i]->intersect(r, localP, localN, localt)) {
                    if (localt < t) {
                        t = localt;
                        P = localP;
                        N = localN;
                        hit_obj = objects[i];

                        intersected = true;
                    }
                }
            }
            return intersected;
        }

        Vector getColor(const Ray &r, int bounce_number) {

            if (bounce_number <= 0) return Vector(0,0,0);

            Vector P,N;
            double t;
            Geometry *hit_obj = nullptr;
            
            if (!intersect(r, P, N, t, hit_obj)) {
                return Vector(0, 0, 0);
            }

            if(hit_obj->mirror) {

                // Handling mirror surfaces
                Vector reflection_dir = r.u - 2*dot(r.u, N)*N;
                Ray mirrorR(P + 0.001*N, reflection_dir);
                return getColor(mirrorR, bounce_number - 1);
            }

            if(hit_obj->transparent) {

                // Handling transparent surfaces without Fresnel's law
                double n1 = 1;
                double n2 = 1.4;
                Vector Nt = N;
                double t_n_coeff;
                Vector t_Tangent;

                if (dot(r.u, N)>0) {
                    std::swap(n1,n2);
                    Nt = -1* Nt;
                }
                t_Tangent = n1 / n2 * (r.u - dot(r.u, Nt) * Nt);
                //t_n_coeff = 1 - sqrt(1 - sqr(n1/n2) * ( 1- sqr(dot(r.u, Nt))));
                double rad = 1 - sqr(n1/n2)* (1-sqr(dot(r.u, Nt)));

                if (rad < 0) {
                    Vector ref_dir = r.u - 2*dot(r.u, N)*N;
                    ref_dir.normalize();
                    Ray mirrorR = Ray(P + 0.001*N, ref_dir);
                    return getColor(mirrorR, bounce_number-1);
                }

                t_n_coeff = -sqrt(rad);
                Vector ray_dir = t_n_coeff * Nt + t_Tangent;
                ray_dir.normalize();
                Ray refraction(P-0.001*Nt, ray_dir);
                return getColor(refraction, bounce_number - 1);
            }

            Vector lightDir = L - P;
            double d2 = lightDir.norm2();
            lightDir.normalize();
            
            Vector color = I/(4 * M_PI * d2) * hit_obj->alb / M_PI * std::max(0., dot(lightDir, N));

            // Checking for if the point is in the shadow
            Ray shadow_ray(P + 0.001*N, lightDir);
            double shadowt;
            Vector shadowP, shadowN;
            Geometry* empty = nullptr;
            bool in_shadow = intersect(shadow_ray, shadowP, shadowN, shadowt, empty);
            bool condition = ((shadowP - P).norm2() <= d2);
            bool both = (in_shadow && condition);
            if (both) {
                color = Vector(0,0,0);
            }

            if (!direct_lighting) {

                // Adding a random indirect light ray
                Vector random_dir = random_lightDir(N);
                random_dir.normalize();
                Ray bounce_ray = Ray(P + 0.0001*random_dir, random_dir);
                Vector new_color = hit_obj->alb * getColor(bounce_ray, bounce_number-1);
                color = (color + new_color) / 2;
            }

            return color;
            }

        ~Scene() {
            for (auto obj : objects) {
                delete obj;
            }
        }

};




int main() {
    int W = 512;
    int H = 512;

    // Defining camera properties
    Vector camera_origin(0, 0, 55); //camera placement
    double fov = 60 * M_PI / 180; //field of view
            
    Scene this_scene;

    //this_scene.objects.push_back(new Sphere(Vector(20,0,0), 10, Vector(1,1,1), false, true));
    //Sphere k(Vector(0,0,0), 10, Vector(0.4,0.7,0.2), false, false);
    //this_scene.add(k);
    //Sphere h(Vector(-20,0,0), 10, Vector(1,0.2,0.8), false, false);
    //this_scene.add(h);

    // Adding walls and floor/ceiling
    this_scene.objects.push_back(new Sphere(Vector(-1000, 0, 0), 940, Vector(0.5, 0.8, 0.1), false, false));
    this_scene.objects.push_back(new Sphere(Vector(1000, 0, 0), 940, Vector(0.9, 0.2, 0.3), false, false));
    this_scene.objects.push_back(new Sphere(Vector(0, 1000, 0), 940, Vector(0.3, 0.5, 0.3), false, false));
    this_scene.objects.push_back(new Sphere(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7), false, false));
    this_scene.objects.push_back(new Sphere(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.17), false, false));
    this_scene.objects.push_back(new Sphere(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9), false, false));
    
    // Adding cat mesh to scene
    TriangleMesh* mesh = new TriangleMesh();
	mesh->readOBJ("cat.obj");
	mesh->alb = Vector(1,1,1);
    mesh->mirror = false;
    mesh->smooth = true;
    mesh->transparent = false;


    // Bounding box
    mesh->scale(0.6, Vector(0,-10,0));
    mesh->compute_bbox();

    this_scene.objects.push_back(mesh);


    // Defining scene lighting
    this_scene.L = Vector(-10,20,40);
    this_scene.I = 8e10;
    this_scene.direct_lighting = false;
    bool antialiasing = true;
    int alias_rays = 16;
    int bounce_number = 5;

    std::vector<unsigned char> image(W * H * 3, 0);

    // iterating through every pixel (i,j)
    #pragma omp parallel for schedule(dynamic,1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            double d = -W/( 2. * tan(fov / 2.) );

            Vector r_dir(j - W/2 + 0.5, -i + 0.5 + H/2, d); // ray directions
            r_dir.normalize(); // normalize ray
            Ray r(camera_origin, r_dir); // ray for this pixel specifically

            Vector color = this_scene.getColor(r, bounce_number);           

            if (antialiasing) {
                color = Vector(0,0,0);
                for (int i = 0; i < alias_rays; i++) {
                    double dx, dy;
                    boxMuller(0.0005, dx, dy);
                    Vector muller_dir = r_dir + Vector(dx, dy, 0);
                    muller_dir.normalize();
                    Ray muller_ray(camera_origin, muller_dir);
                    
                    Vector random_color = this_scene.getColor(muller_ray, bounce_number);
                    color = color + random_color;
                }
                color = color / (this_scene.ray_count);
            }
            
            image[(i * W + j) * 3 + 0] =(std::max(0., std::min(255., std::pow(color.data[0], 1/2.2))));
            image[(i * W + j) * 3 + 1] = (std::max(0., std::min(255., std::pow(color.data[1], 1/2.2))));
            image[(i * W + j) * 3 + 2] =(std::max(0., std::min(255., std::pow(color.data[2], 1/2.2))));

        }

        }
    stbi_write_png("lab_image.png", W, H, 3, &image[0], 0);
 
    return 0;
}