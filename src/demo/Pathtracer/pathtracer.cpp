#include <stdlib.h>
#ifdef __APPLE__
#include <GLUT/GLUT.h>
#else
#include <GL/glut.h>
#endif

#include <GEL/CGLA/Vec2i.h>
#include <GEL/CGLA/Vec3uc.h>

#include "scene.h"
#include "camera.h"
#include "mesh.h"
#include "omni.h"

#include "matte.h"
#include "plastic.h"
#include "metal.h"
#include "glass.h"

#include "pathtracer.h"

using namespace CGLA;

//global pointer to the active scene
scene* current;

static pathtracer* renderer;

const int width = 64*8;
const int height = 64*8;

static Vec3f* film;
static Vec3uc* image;
static Vec2i pixel(0);
static bool done = false;

static float dgamma = 2.2f;
static float dexposure = 0.f;


pathtracer::pathtracer(int w, int h, bool explicit_direct, int subsamples)
: width_(w), height_(h), scene_(0)
{
    explicit_direct_ = explicit_direct;
    subsamples_ = subsamples;
}

void pathtracer::set_scene(scene* s)
{
    scene_ = s;
}

CGLA::Vec3f pathtracer::trace(const ray& r, bool include_emitted)
{
    //intersect ray with
    hit_info hi;
    bool hit = scene_->intersect(r, hi);

    if (!hit)
        return Vec3f(0.f, 0.f, 0.f);

    //Vec3f x, y, z = hi.shading_normal;
    //orthogonal(z, x, y);
    //return (Vec3f(hi.texcoords(0),hi.texcoords(1),0) + Vec3f(0.f)) / 1.f;

    //only include emitted light if requested
    CGLA::Vec3f Le(0.f);
    if (include_emitted)
        Le = hi.emitted;

    //compute reflectance for each bsdf component
    float rho_diffuse = intensity(hi.diffuse);
    float rho_glossy = intensity(hi.glossy);
    float rho_reflection = intensity(hi.reflection);
    float rho_refraction = intensity(hi.refraction);
    float rho_total = rho_diffuse+rho_glossy+rho_reflection+rho_refraction;
    assert(rho_total < 1.f);

    if (rho_total == 0.f)
        return Le;

    //compute direct lighting on hit point
    CGLA::Vec3f Ld(0.f);
    Vec3f wo = -r.direction;
    if (explicit_direct_ && rho_diffuse+rho_glossy>0.f)
    {
        size_t nlums = scene_->luminaires();
        for (size_t i=0; i<nlums; ++i)
        {
            const luminaire* lum = scene_->get_luminaire(i);
            int samples = lum->samples();

            Vec3f L(0.f);
            for (int j=0; j<samples; ++j)
            {
                Vec3f Li, wi;
                if (lum->sample(r, hi, Li, wi))
                {
                    float cost = std::max(dot(hi.shading_normal, wi),0.f);
                    L += Li * cost * bsdf_evaluate(hi, wi, wo);
                }

                Ld += L / float(samples);
            }
        }
    }

    //use russian roulette to terminate path
    float prussian = 1.f;
    float rr = mt_random();

    if (r.depth > 3)
    {
        prussian = rho_total;

        if (rr >= prussian)
            return Ld + Le;
    }

    //figure out which bXdf to sample
    float pdiffuse = rho_diffuse/rho_total;
    float pglossy = rho_glossy/rho_total;
    float preflection = rho_reflection/rho_total;
    float prefraction = rho_refraction/rho_total;

    Vec3f wi;
    float pwi;
    Vec3f fs;
    bool sample_emitted = !explicit_direct_;
    if (rr <= pdiffuse)
    {
        //sample diffuse
        float pwi = sample_lambertian(hi, wo, wi);
        fs = lambertian_brdf(hi, wi, wo) / (pwi * pdiffuse * prussian);
    }
    else if (rr <= pdiffuse+pglossy)
    {
        //sample glossy part
        float pwi = sample_phong(hi, wo, wi);

        if (dot(hi.shading_normal, wi) <= 0.f)
            return Le + Ld;

        fs = phong_brdf(hi, wi, wo) / (pwi * pglossy * prussian);
    }
    else if (rr <= pdiffuse+pglossy+preflection)
    {
        //sample perfect specular reflection
        wi = reflect(hi.shading_normal, wo);
        pwi = 1.f;
        float cost = dot(hi.shading_normal, wo);
        fs = hi.reflection / (pwi * preflection * cost);
        sample_emitted = true;
    }
    else if (rr <= pdiffuse+pglossy+preflection+prefraction)
    {
        //sample perfect specular refraction
        bool not_tir = refract(hi.shading_normal, wo, 1.f/hi.ior, wi);
        assert(not_tir);
        pwi = 1.f;
        float cost = std::abs(dot(hi.shading_normal, wi));
        fs = hi.refraction / (pwi * prefraction * cost);
        sample_emitted = true;
    }
    else
        assert(false);

    //create aux ray
    ray s;
    s.origin = hi.position + epsilon * wi;
    s.direction = wi;
    s.depth = r.depth + 1;
    s.distance = std::numeric_limits<float>::infinity();
    float cost = std::abs(dot(hi.shading_normal, wi));
    Vec3f Li = cost * fs * trace(s, sample_emitted);

    //returm sum of emitted + direct + indirect
    return Le + Ld + Li;
}

Vec3f pathtracer::compute_pixel(int w, int h)
{
    assert(scene_);

    //supersample
    Vec3f L(0.f);
    for (int j=0; j<subsamples_; ++j)
    {
        float y  = h + (j + 0.5f) / subsamples_;

        for (int i=0; i<subsamples_; ++i)
        {
            float x  = w + (i + 0.5f) / subsamples_;

            //ask camera for initial ray..
            Vec2f uv(x/width_, y/height_);
            ray r = scene_->get_camera()->generate(uv);

            //trace ray
            L += trace(r, true);
        }
    }

    return L / float(subsamples_ * subsamples_);
}


Vec3uc tonemap(const Vec3f& v)
{
    Vec3f u = v * std::pow(2.f, dexposure);
    float r = mt_random() - 0.5f;

    Vec3uc I;
    for (int i=0; i<3; ++i)
    {
        float c = std::pow(u[i], 1.f / dgamma);
        c = 255.f * c + 0.5f + r;
        I[i] = clamp(int(c), 0, 255);
    }
    return I;
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);

    for (int j=0; j<height; ++j)
        for (int i=0; i<width; ++i)
            image[i + j*width] = tonemap(film[i + j*width]);

    glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, image);
    glutSwapBuffers();

    GLenum e = glGetError();
    if (e != GL_NO_ERROR)
    {
        printf("OpenGL error: %s\n", gluErrorString(e));
        exit(EXIT_FAILURE);
    }
}

void reshape(GLint width, GLint height)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, width, 0, height);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void special(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_LEFT:
        dgamma = std::max(dgamma - 0.1f, 0.5f);
        break;

    case GLUT_KEY_RIGHT:
        dgamma = std::min(dgamma + 0.1f, 5.f);
        break;

    case GLUT_KEY_DOWN:
        dexposure = std::max(dexposure - 0.1f, -10.f);
        break;

    case GLUT_KEY_UP:
        dexposure = std::min(dexposure + 0.1f, 10.f);
        break;
    }

    printf("gamma: %.2f, exposure: %.2f\n", dgamma, dexposure);
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 27: //ESCAPE key
        exit(0);
        break;
    }
}

void idle(void)
{
    for (int i=0; i<width*8 && !done; ++i)
    {
        Vec3f L = renderer->compute_pixel(pixel[0], pixel[1]);
        film[pixel[0] + pixel[1] * width] = L;

        //advance
        ++pixel[0];

        if (pixel[0] == width)
        {
            pixel[0] = 0;
            ++pixel[1];

            if (pixel[1] == height)
            {
                done = true;
                glutIdleFunc(NULL);
            }
        }
    }

    glutPostRedisplay();
}

int main(int argc, char* argv[])
{
    //make the scene
    current = new scene;

    //setup some materials
    matte dull_white(Vec3f(0.5f, 0.5f, 0.5f));
    matte dull_gray(Vec3f(0.4f, 0.4f, 0.4f));
    matte dull_red(Vec3f(0.6f, 0.3f, 0.2f));
    matte dull_green(Vec3f(0.3f, 0.6f, 0.2f));
    matte dull_blue(Vec3f(0.2f, 0.2f, 0.6f));

    plastic glossy_white(Vec3f(0.4f, 0.4f, 0.4f), Vec3f(0.5f), 20.f);
    plastic glossy_purple(Vec3f(0.3f, 0.1f, 0.3f), Vec3f(0.65f), 8.f);
    plastic glossy_yellow(Vec3f(0.3f, 0.3f, 0.0f), Vec3f(0.65f), 32.f);

    //for exercise 2
    glass clear(Vec3f(0.99f, 0.99f, 0.99f), 1.5f);
    metal silver(Vec3f(0.9f, 0.9f, 0.9f), 22.f, 0.177f, 3.638f);

    //setup cornell box 1mx1mx1m
    mesh floor("../../data/cornell_box/floor.obj");
    floor.set_material(&glossy_white);
    current->insert(&floor);

    mesh ceiling("../../data/cornell_box/ceiling.obj");
    ceiling.set_material(&dull_gray);
    current->insert(&ceiling);

    mesh back("../../data/cornell_box/back.obj");
    back.set_material(&silver);
    current->insert(&back);

    mesh left("../../data/cornell_box/left.obj");
    left.set_material(&dull_red);
    current->insert(&left);

    mesh right("../../data/cornell_box/right.obj");
    right.set_material(&dull_blue);
    current->insert(&right);

    //add some objects to the box
    Mat4x4f tmp;

    mesh box1("../../data/cornell_box/box1.obj");
    tmp = rotation_Mat4x4f(YAXIS, -float(M_PI)/8.f);
    tmp = translation_Mat4x4f(Vec3f(-0.3f,0.f,-0.05f)) * tmp;
    box1.set_transform(tmp);
    box1.set_material(&dull_green);
    current->insert(&box1);

    mesh box2("../../data/cornell_box/box2.obj");
    tmp = translation_Mat4x4f(Vec3f(0.25f,0.f, -0.05f));
    box2.set_transform(tmp);
    box2.set_material(&glossy_purple);
    current->insert(&box2);

    mesh teapot("../../data/cornell_box/teapot1.obj");
    tmp = rotation_Mat4x4f(YAXIS, float(M_PI/4.f));
    tmp = translation_Mat4x4f(Vec3f(-0.25f,0.25f, -0.05)) * tmp;
    teapot.set_transform(tmp);
    teapot.set_material(&silver);
    current->insert(&teapot);

    mesh sphere1("../../data/cornell_box/sphere2.obj");
    tmp = translation_Mat4x4f(Vec3f(0.25f,0.45f,-0.05f));
    sphere1.set_transform(tmp);
    sphere1.set_material(&clear);
    current->insert(&sphere1);

    mesh sphere2("../../data/cornell_box/sphere2.obj");
    tmp = translation_Mat4x4f(Vec3f(0.22f,0.15f,0.25f));
    sphere2.set_transform(tmp);
    sphere2.set_material(&glossy_yellow);
    //current->insert(&sphere2);

    //light sources
    mesh quad_light("../../data/cornell_box/quad.obj");
    tmp = rotation_Mat4x4f(XAXIS, float(M_PI));
    tmp = translation_Mat4x4f(Vec3f(0.f,0.99f,0.2f)) * tmp;
    quad_light.set_transform(tmp);
    quad_light.set_exitance(Vec3f(300,300,300));
    current->insert(&quad_light);

    mesh quad_light1("../../data/cornell_box/quad.obj");
    tmp = rotation_Mat4x4f(XAXIS, float(M_PI));
    tmp = translation_Mat4x4f(Vec3f(0.f,0.99f,0.2f)) * tmp;
    quad_light1.set_transform(tmp);
    quad_light1.set_exitance(Vec3f(0,400,0));
//  current->insert(&quad_light1);

    mesh quad_light2("../../data/cornell_box/quad.obj");
    tmp = rotation_Mat4x4f(XAXIS, float(M_PI));
    tmp = translation_Mat4x4f(Vec3f(-0.2f,0.99f,0.2f)) * tmp;
    quad_light2.set_transform(tmp);
    quad_light2.set_exitance(Vec3f(400,0,0));
//  current->insert(&quad_light2);

    mesh quad_light3("../../data/cornell_box/quad.obj");
    tmp = rotation_Mat4x4f(XAXIS, float(M_PI));
    tmp = translation_Mat4x4f(Vec3f(0.2f,0.99f,0.2f)) * tmp;
    quad_light3.set_transform(tmp);
    quad_light3.set_exitance(Vec3f(0,0,400));
//  current->insert(&quad_light3);

    omni omni_light(Vec3f(30.f));
    omni_light.set_transform(translation_Mat4x4f(Vec3f(0.0f,0.95f,0.0f)));
    //current->insert(&omni_light);

    mesh sphere_light("../../data/cornell_box/small_sphere.obj");
    sphere_light.set_transform(translation_Mat4x4f(Vec3f(0.f,0.95f,0.f)));
    sphere_light.set_exitance(Vec3f(30.f/(4.f*float(M_PI)*0.01f*0.01f)));
    //current->insert(&sphere_light);

    //setup camera (eye, center, up, focal length)
    camera pentax(
        Vec3f(0.f,0.5f,2.0f),
        Vec3f(0.f,0.5f,0.5f),
        Vec3f(0,1,0),
        0.035f);

    current->set_camera(&pentax);

    //build acceleration structure
    current->initialize(8, 25);

    //create the renderer
    renderer = new pathtracer(width, height, true, 1);
    renderer->set_scene(current);

    //create the film
    film = new Vec3f[width * height];
    std::fill(film, film+width*height, Vec3f(0.3f));
    image = new Vec3uc[width * height];
    std::fill(image, image+width*height, Vec3uc(32,32,32));

    //init glut
    glutInit(&argc, argv);
    glutInitWindowSize(width, height);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Path tracer");

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);

    glutIdleFunc(idle);

    //Turn the flow of control over to GLUT
    glutMainLoop();

    //clean up
    delete current;
    delete renderer;

    return EXIT_SUCCESS;
}
