#ifndef __PATHTRACER_SCENE__H__
#define __PATHTRACER_SCENE__H__

#include "core.h"
#include "luminaire.h"
#include "camera.h"

namespace Geometry { class BSPTree; }

class scene
{
public:
    scene(void);
    ~scene(void);

    //insert object into scene (all objects are potential luminaires)
    size_t luminaires(void) const;
    void insert(const luminaire*);
    const luminaire* get_luminaire(size_t) const;

    //set the one and only camera
    const camera* get_camera(void) const;
    void set_camera(const camera*);

    //initialize acceleration structure
    void initialize(int max_objs=1<<20, int max_depth=20);

    //intersection testing
    bool intersect(const ray&);
    bool intersect(const ray&, hit_info&);

private:
    const camera* cam_;
    std::vector<const luminaire*> objs_;
    Geometry::BSPTree* accel_;
};

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
